import tempfile
import sys, os
import gzip
from mOTUlizer.config import *
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from os.path import join as pjoin
import shutil

methods =  ['silixCOGs', 'mmseqsCluster']
class GeneClusters():
    def __repr__(self):
        return f"< a GeneClusters-ing of {len(self)} GeneClusters from {len(self._motus)} mOTUs>"

    def __init__(self, motus, gene2gene_clust = None, gene2genome = None, gene_clust2rep = dict(), ):
        self._motus = motus
        self.genomes = [g for motu in motus for g in motu]
        self.gene_clusters = []
        if gene2gene_clust:
            all_gene_clust = {s : [] for s in set(gene2gene_clust.values())}
            for gene, clust in gene2gene_clust.items():
                all_gene_clust[clust].append(gene)
            self.gene_clusters = [GeneCluster(self, genes, gene_clust2rep.get(clust, None)) for clust, genes in all_gene_clust.items()]
        self.gene2genome = gene2genome

    def __getitem__(self, i):
        if type(i) == int and i < len(self):
            return self.gene_clusters[i]
        else:
            raise KeyError(str(i) + "is not a valid entry")

    def __len__(self):
        return len(self.gene_clusters)


    def __iter__(self):
       return GeneClustersIterator(self)


    @classmethod
    def compute_GeneClusters(cls, motus, name, precluster = False, threads = 4, method =  "mmseqsCluster", **kwargs):
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
            genome2gene_clusters = {k : set()  for motu in motus for k in motu}
            for k,v in recs.items():
                for vv in prot2genome[k]:
                    genome2gene_clusters[vv].update([v])

        else :
            print("The '{}' clustering method is not implemented yet".format(name) , file = sys.stderr)
            print("only allowed are :", methods, file = sys.stderr)
            sys.exit(-1)


        shutil.rmtree(temp_folder)

        for k,v in genome2gene_clusters.items():
            genome2gene_clusters[k] = set(v)


        return GeneClusters(motus, recs, prot2genome , gene_clusters2rep)
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

    def __init__(self, gene_clusters, genes, representative = None):
        self._gene_clusters = gene_clusters
        self.genes = genes
        self.representative = representative
