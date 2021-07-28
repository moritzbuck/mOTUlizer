import tempfile
import sys, os
import gzip
from mOTUlizer.config import *
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from os.path import join as pjoin
import shutil
from mOTUlizer.errors import *
from mOTUlizer.utils import random_name
from mOTUlizer.classes.tools.Muscle import Muscle
import json

methods =  ['silixCOGs', 'mmseqsCluster']
class GeneClusters():
    def __repr__(self):
        return f"< a GeneClusters-ing of {len(self)} GeneClusters from {len(self.genomes)} genomes>"

    def __init__(self, motus_or_genome, gene_clusters = None, storage = None, quiet = False):
        if type(motus_or_genome) == list:
            self._from_motu = True
            if not all([hasattr(motu, "can_haz_gene_clusters") for motu in motus_or_genome]):
                raise CantGeneClusterError(f"you need to pass either an mOTU -list or a single genome,\n you passed a list of things where something don't work \n you passed {motus_or_genome}\n this is dirty, I am sorry, too lazy to implement interfaces")
            self.genomes = [g for motu in motus_or_genome for g in motu]
        else :
            self._from_motu = False
            if not hasattr(motus_or_genome, "can_haz_gene_clusters"):
                raise CantGeneClusterError(f"you need to pass either an mOTU-list or a genome as motus_or_genome,\n you passed a single thing that don't work\ the thing is {motus_or_genome}")
            self.genomes = [motus_or_genome]

        self.quiet = quiet
        self.storage = storage

        if self.storage:
            if not os.path.exists(self.storage):
                if not self.quiet:
                    print(f"Creating {self.storage} as storage folder", sys.stderr)
                    os.makedirs(storage, exist_ok=False)

        if gene_clusters:
            self.gene_clusters = [g for g in gene_clusters]
            if self._from_motu:
                cluster_sets = { g : [] for g in self.genomes}
                for clust in self.gene_clusters:
                    clust.set_storage(self.storage)
                    clust.quiet = self.quiet
                    for genom in clust.get_genomes():
                        cluster_sets[genom].append(clust)
                for g in self.genomes :
                    g.set_gene_clusters(cluster_sets[g])

        self._clust_set = None



    def get_genecluster_counts(self):
        return {g : len(g.get_genomes()) for g in self}

    def get_clust_set(self) :
        if not self._clust_set:
            self._clust_set =  set(self.gene_clusters)
        return self._clust_set

    def __getitem__(self, i):
        if type(i) == int and i < len(self):
            return self.gene_clusters[i]
        else:
            raise KeyError(str(i) + " is not a valid entry for your GeneCluster, try the name of the GC.")

    def __len__(self):
        return len(self.gene_clusters)


    def __iter__(self):
       return GeneClustersIterator(self)

    def item(self):
        return GeneClustersItemIterator(self)

    def export_gene_clusters(self, file, force = False):
        data = { c.name : c.export_data() for c in self}
        with open(file, "w") as handle:
            json.dump(data, handle, indent = 2, sort_keys = True)

    @classmethod
    def load(cls, file, motu, force = False, storage = None):
        if not os.path.exists(file):
            raise FileError("Your gene-cluster file does not exist or something")
        try :
            with open(file) as handle:
                data = json.load(handle)
                json_data = True
        except :
            try :
                with open(file) as handle:
                    data = {l.split("\t")[0] : l[:-1].split("\t")[1:] for l in handle}
                    json_data = False
            except :
                raise CantLoadGeneClusterError("your genes-clusters are neither json- or tsv- formated")

        if json_data and all(["genes2genome" in v for v in data.values()]):
            cluster_names = list(data.keys())
            genomes = {genome for k,v in data.items() for gene, genomes in v['genes2genome'].items() for genome in genomes}
            genes2genome = {k : {gene : [motu.get(genome) for genome in genomes] for gene, genomes in v['genes2genome'].items()} for k,v in data.items() }
            representatives = {k : v.get('representative') for k,v in genes2genome.items()}
        else :
            print("Running the vintage parser! (e.g. for gene-clusters save with mOTUpan v < 0.4.0)" , file = sys.stderr)
            cluster_names = list({vv for v in data.values() for vv in v})
            genomes = {k for k in data}
            clust2genome = {c : []  for c in cluster_names}
            for k,v in data.items():
                for vv in v:
                    clust2genome[vv].append(k)
            genes2genome = {k : {"MockGene_" + genome : [motu.get(genome)] for genome in v} for k,v in clust2genome.items() }
            representatives = {k : None for k in cluster_names}

            if not all([motu.get(genome) for genome in genomes]):
                if force:
                    genes2genome = {k : {gene : [genome for genome in genomes if genome] for gene, genomes in v.items()} for k,v in genes2genome.items() }
                else :
                    bads = [genome for genome in genomes if not motu.get(genome)]
                    raise CantLoadGeneClusterError(f"can't find some genome-IDs of your clusterings in your mOTU, I can filter them out but you need to explicitly say 'force=True'\n the missing genomes are {';'.join(bads)}")

        clusters = [GeneCluster(name = k, genes2genome = genes2genome[k], representative = representatives[k]) for k in cluster_names]
        return GeneClusters([motu], gene_clusters = clusters, storage = storage)

    @classmethod
    def compute_GeneClusters(cls, motus, name, precluster = False, threads = 4, storage = None, method =  "mmseqsCluster", **kwargs):
        name = name + method + "_"

        temp_folder = tempfile.mkdtemp(prefix = name)
        all_faas_file = pjoin(temp_folder, "concat.faa")
        gene_clusters2rep = None
        prot_ids = set()
        prot2genome ={}
        for motu in motus:
            for genome in motu:
                aminos = genome.get_amino_acids()
                for i in aminos:
                    if i in prot2genome:
                        prot2genome[i] += [genome]
                    else :
                        prot2genome[i] = [genome]
                with open(all_faas_file, "a") as handle:
                    SeqIO.write([SeqRecord(id = id, seq = seq, description = "") for id, seq in aminos.items()], handle, "fasta")

        if precluster:
            cdhit_file = tempfile.NamedTemporaryFile().name

            if not shutil.which('cd-hit'):
                print("You need cd-hit to run the preclustering, either install it or run mOTUpan without preclustering", file = sys.stderr)
                sys.exit(-1)

            exec = "cd-hit -i {input} -o {output} -c 0.95 -M 0 -T {threads} -d 0 -s 0.95 >> {log} 2>&1".format(input = all_faas_file, output = cdhit_file, threads = threads, log = "/dev/null")
            print("Running cd-hit preclustering", file = sys.stderr)
            os.system(exec)

            print("parsing cd-hit", file = sys.stderr)

            with open(cdhit_file + ".clstr") as handle:
                clusters = "\n".join(handle.readlines()).split("Cluster ")

            os.remove(cdhit_file + ".clstr")
            os.remove(cdhit_file)

            clusters = [c.split("\n\n") for c in clusters[1:] if "*" in c]
            clusters = [[cc.split(">")[1].split("... ") for cc in c if ">" in cc and cc != ">"] for c in clusters ]
            clusters = {[cc[0] for cc in c if cc[1] == "*" or cc[1] == "*\n"][0] : [cc[0] for cc in c] for c in clusters}

            print("For", len(prot2genome), "CDSes we got ", len(clusters), "preclusters", file = sys.stderr)
            seqs = [s for s in SeqIO.parse(all_faas_file, "fasta") if s.id in clusters]
            SeqIO.write(seqs, all_faas_file, "fasta")

        if method == "silixCOGs" :
            if not shutil.which("diamond") or not shutil.which("silix"):
                print("You need diamond and silix to run the silix gene-clustering, either install it or run mOTUpan with an other gene-clustering or your own traits", file = sys.stderr)
                sys.exit(-1)

            print("all v all diamond for silix", file = sys.stderr)

            os.system("diamond makedb --db {faas} --in {faas} > /dev/null 2> /dev/null".format(faas = all_faas_file))
            os.system("diamond blastp --more-sensitive -p {threads} -f 6 -q {faas} --db {faas} -o {out} 2> /dev/null > /dev/null".format(faas = all_faas_file, out = temp_out, threads = threads))

            print("running silix", file = sys.stderr)
            os.system("silix {faas} {out} > {clust_temp} #2> /dev/null".format(faas = all_faas_file, out = temp_out, clust_temp = temp_clust))

            print("parsing silix", file = sys.stderr)
            with open(temp_clust) as handle:
                if precluster:
                    recs = {g : l[:-1].split()[0]  for l in handle for g in clusters[l[:-1].split()[1]]}
                else :
                    recs = {l[:-1].split()[1] : l[:-1].split()[0]  for l in handle}


            #pretty formating names
            fill = max([len(v) for v in recs.values()])
            recs = {k : name + "_" + v.zfill(fill) for k, v in recs.items()}

            genome2gene_clusters = {k : set() for motu in motus for k in motu}
            for k,v in recs.items():
                for vv in prot2genome[k]:
                    genome2gene_clusters[vv].update([v])

        elif method == "mmseqsCluster" :
            "coverage = 80% with cov-mode = 0, minimal amino acid sequence identity = 80% and cluster-mode = 0"
            covmode = kwargs.get("mmseqsCluster_covmode",0)
            cov = kwargs.get("mmseqsCluster_cov",0.80)
            seqid = kwargs.get("mmseqsCluster_seqid",0.80)

            if not shutil.which("mmseqs"):
                print("You need mmseqs2 to run the silix gene-clustering, either install it or run mOTUpan with an other gene-clustering or your own traits", file = sys.stderr)
                sys.exit(-1)

            print("running mmseqs easy-cluster with params --min-seq-id {seqid} --cov-mode {covmode} -c {cov} in {}".format(temp_folder, covmode = covmode, cov = cov, seqid=seqid), file = sys.stderr)
            mmseqs_dat = pjoin(temp_folder, "mmseqs_")
            os.system("mmseqs easy-cluster --threads {threads} --min-seq-id {seqid} --cov-mode {covmode} -c {cov} {faas} {out} {tmp} 2> /dev/null > /dev/null".format(covmode = covmode, cov = cov, seqid=seqid, faas = all_faas_file, out = mmseqs_dat, tmp = temp_folder, threads = threads))

            with open(mmseqs_dat + "_cluster.tsv") as handle:
                if precluster:
                    recs = {g : l[:-1].split()[0]  for l in handle for g in clusters[l[:-1].split()[1]]}
                else :
                    recs = {l[:-1].split()[1] : l[:-1].split()[0]  for l in handle}

            #pretty formating names
            fill = len(str(len(set(recs.values()))))

            rep2clust = {k : name + str(i).zfill(fill) for i,k in enumerate(set(recs.values()))}
            gene_clusters2rep = {v: k for k,v in rep2clust.items()}

            print("For", len(recs), "CDSes we got ", len(gene_clusters2rep), " gene-clusters", file = sys.stderr)

            recs = {k : rep2clust[v] for k, v in recs.items()}
            clst2gene = { clst : [] for clst in gene_clusters2rep}
            for gene, clst in recs.items():
                clst2gene[clst].append(gene)
        else :
            print("The '{}' clustering method is not implemented yet".format(name) , file = sys.stderr)
            print("only allowed are :", methods, file = sys.stderr)
            sys.exit(-1)


        shutil.rmtree(temp_folder)

        clusts = [GeneCluster({g : prot2genome[g] for g in genes}, gene_clusters2rep[clust], name = clust, storage = None if not storage else (storage + "clusters/") ) for clust, genes in clst2gene.items()]
        return GeneClusters(motus, clusts, storage )
#        return { 'genome2gene_clusterss' : genome2gene_clusters, 'aa2gene_clusters' : recs, 'gene_clusters2rep' : gene_clusters2rep}

class GeneClustersIterator:
    def __init__(self, gene_clusters):
        self._gene_clusters = gene_clusters
        self.index = -1

    def __next__(self):
        if self.index < len(self._gene_clusters)-1:
            self.index += 1
            return self._gene_clusters[self.index]
        raise StopIteration

class GeneCluster():
    def __repr__(self):
        return f"< a GeneCluster of {len(self.genes)} genes from {len(self.get_genomes())} genomes>"

    def __init__(self, genes2genome, representative = None, name = None, storage = None, quiet = False):
        """write chekcs here"""
        if not name:
            self.name = "GeneCluster_" +random_name()
        else :
            self.name = name
        self.genes = list(genes2genome.keys())
        self.representative = representative
        self.gene2genome = genes2genome
        self._genomes = None
        self.set_storage(storage)
        self.quiet = quiet

    def set_storage(self, storage):
        if storage:
            if storage.endswith(self.name):
                self.storage = storage
            else :
                self.storage = pjoin(storage, self.name)
            os.makedirs(self.storage, exist_ok = True)
        else :
            self.storage = None

    def __len__(self):
        return len(self.genes)

    def export_data(self):
        return {'genes2genome' :
                    {gene : [g.name for g in genomes] for gene, genomes in self.gene2genome.items() },
                 'representative' : self.representative
               }

    def get_genomes(self):
        if not self._genomes:
            self._genomes = {gg for g in self.gene2genome.values() for gg in g}
        return self._genomes

    def get_seqs(self):
        seqs = dict()
        for gene, genomes in self.gene2genome.items():
            seqs[gene] = genomes[0].get_aa(gene)
        return seqs

    def get_alignment(self, method = "muscle"):
        if self.storage:
            out_file = pjoin(self.storage, "alignment.fasta")
        else :
            out_file = None
        if method == "muscle":
            tool = Muscle(seqsdict_or_infile = self.get_seqs(), out_file = out_file, quiet = True)
        else :
            raise CantMethodError(f"The alignment method ({method}) you asked for is not implemented yet.")

        if out_file and not os.path.exists(out_file):
            tool.run_command()
        elif not self.quiet:
            print("Alignment loaded from file", file = sys.stderr)

        ali = tool.parse_output()
        return ali

    def ali2codon(self):
        codoned = {gid for gid in self.genes}
        for gid, allied in self.get_alignment().items():
            seq = self.gene2genome[gid][0].get_gene(gid)
        return seq