import tempfile
import sys, os
import gzip
from mOTUlizer.config import *
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import itertools
from os.path import join as pjoin
import shutil
import mOTUlizer
from mOTUlizer.errors import *
from mOTUlizer.utils import random_name, message
from mOTUlizer.classes.tools.Muscle import Muscle
from mOTUlizer.classes.tools.FastTree import FastTree
from mOTUlizer.db.SeqDb import SeqDb
import json
import re
from tqdm import tqdm
from mOTUlizer import get_quiet, get_threads

methods =  ['silixCOGs', 'mmseqsCluster', "ppanggolin"]
# class GeneClusters():
#     def __repr__(self):
#         return f"< a GeneClusters-ing of {len(self)} GeneClusters from {len(self.genomes)} genomes>"
#
#     def __init__(self, motus_or_genome, clustering_name = None, gene_clusters = None, quiet = False):
#         if not SeqDb.seq_db:
#             raise DataBaseNotInitialisedError("The database has not been initialised")
#         self.db = SeqDb.get_global_db()
#         if clustering_name and db.has_clustering():
#             self.genomes
#             self.gene_clusters = [GeneCluster(n,)]
#         if type(motus_or_genome) == list:
#             self._from_motu = True
#             if not all([hasattr(motu, "can_haz_gene_clusters") for motu in motus_or_genome]):
#                 raise CantGeneClusterError(f"you need to pass either an mOTU -list or a single genome,\n you passed a list of things where something don't work \n you passed {motus_or_genome}\n this is dirty, I am sorry, too lazy to implement interfaces")
#             self.genomes = [g for motu in motus_or_genome for g in motu]
#         else :
#             self._from_motu = False
#             if not hasattr(motus_or_genome, "can_haz_gene_clusters"):
#                 raise CantGeneClusterError(f"you need to pass either an mOTU-list or a genome as motus_or_genome,\n you passed a single thing that don't work\ the thing is {motus_or_genome}")
#             self.genomes = [motus_or_genome]
#
#         self.quiet = quiet
#
#         if gene_clusters:
#             self.gene_clusters = [g for g in gene_clusters]
#             if self._from_motu:
#                 cluster_sets = { g : [] for g in self.genomes}
#                 for clust in self.gene_clusters:
#                     clust.set_storage(self.storage)
#                     clust.quiet = self.quiet
#                     for genom in clust.get_genomes():
#                         cluster_sets[genom].append(clust)
#                 for g in self.genomes :
#                     g.set_gene_clusters(cluster_sets[g])
#
#         self._clust_set = None
#
#
#     def intersection(self, other_cluster):
#         return self.get_clust_set().intersection(other_cluster.get_clust_set())
#
#     def get_genecluster_counts(self):
#         return {g : len(g.get_genomes()) for g in self}
#
#     def get_clust_set(self) :
#         if not self._clust_set:
#             self._clust_set =  set(self.gene_clusters)
#         return self._clust_set
#
#     def __getitem__(self, i):
#         if type(i) == int and i < len(self):
#             return self.gene_clusters[i]
#         else:
#             raise KeyError(str(i) + " is not a valid entry for your GeneCluster, try the name of the GC.")
#
#     def __len__(self):
#         return len(self.gene_clusters)
#
#
#     def __iter__(self):
#        return GeneClustersIterator(self)
#
#     def item(self):
#         return GeneClustersItemIterator(self)
#
#     def export_gene_clusters(self, file, force = False):
#         data = { c.name : c.export_data() for c in self}
#         with open(file, "w") as handle:
#             json.dump(data, handle, indent = 2, sort_keys = True)
#
#     @classmethod
#     def load(cls, file, motu, force = False, storage = None):
#         if not os.path.exists(file):
#             raise FileError("Your gene-cluster file does not exist or something")
#         try :
#             with open(file) as handle:
#                 data = json.load(handle)
#                 json_data = True
#         except :
#             try :
#                 with open(file) as handle:
#                     data = {l.split("\t")[0] : l[:-1].split("\t")[1:] for l in handle}
#                     json_data = False
#             except :
#                 raise CantLoadGeneClusterError("your genes-clusters are neither json- or tsv- formated")
#
#         if json_data and all(["genes2genome" in v for v in data.values()]):
#             cluster_names = list(data.keys())
#             genomes = {genome for k,v in data.items() for gene, genomes in v['genes2genome'].items() for genome in genomes}
#             genes2genome = {k : {gene : [motu.get(genome) for genome in genomes] for gene, genomes in v['genes2genome'].items()} for k,v in data.items() }
#             representatives = {k : v.get('representative') for k,v in genes2genome.items()}
#         else :
#             print("Running the vintage parser! (e.g. for gene-clusters save with mOTUpan v < 0.4.0)" , file = sys.stderr)
#             cluster_names = list({vv for v in data.values() for vv in v})
#             genomes = {k for k in data}
#             clust2genome = {c : []  for c in cluster_names}
#             for k,v in data.items():
#                 for vv in v:
#                     clust2genome[vv].append(k)
#             genes2genome = {k : {"MockGene_" + genome : [motu.get(genome)] for genome in v} for k,v in clust2genome.items() }
#             representatives = {k : None for k in cluster_names}
#
#             if not all([motu.get(genome) for genome in genomes]):
#                 if force:
#                     genes2genome = {k : {gene : [genome for genome in genomes if genome] for gene, genomes in v.items()} for k,v in genes2genome.items() }
#                 else :
#                     bads = [genome for genome in genomes if not motu.get(genome)]
#                     raise CantLoadGeneClusterError(f"can't find some genome-IDs of your clusterings in your mOTU, I can filter them out but you need to explicitly say 'force=True'\n the missing genomes are {';'.join(bads)}")
#
#         clusters = [GeneCluster(name = k, genes2genome = genes2genome[k], representative = representatives[k]) for k in cluster_names]
#         return GeneClusters([motu], gene_clusters = clusters, storage = storage)



def compute_GeneClusters(motu_or_genome_list, name, precluster = False, method =  "mmseqsCluster", db = None,  **kwargs):
    name = name + method + "_"
    threads = get_threads()
    temp_folder = tempfile.mkdtemp(dir = mOTUlizer._temp_folder_ , prefix = name.replace(";", "_"))
    all_faas_file = pjoin(temp_folder, "concat.faa")
    if not db:
        db = SeqDb.get_global_db()



    gene_clusters2rep = None
    prot_ids = set()
    prot2genome ={}
    ppang_part = {}
    if hasattr(motu_or_genome_list, "_gene_clusters"):
        delattr(motu_or_genome_list, "_gene_clusters")
    for genome in motu_or_genome_list:
        if method != "ppanggolin":
            if not genome.db.has_features(genome.name):
                raise GenomeNotAnnotated("Your genome has no features, it probably wans't annotated")
            for k,seq in genome.amino_acids.items():
                prot2genome[str(k)] = genome
            with open(all_faas_file, "a") as handle:
                handle.writelines([f">{k}\n{seq}\n" for k,seq in genome.amino_acids.items()])


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
    elif method == "ppanggolin":
        from mOTUlizer.classes.tools.Ppanggolin import Ppanggolin

        pang = Ppanggolin(motu_or_genome_list)
        if not get_quiet():
            print("Running Ppanggolin  'annotate'")
        pang.run_command()
        data = pang.parse_output()
        clst2gene = data['clst2feature_id']
        gene_clusters2rep = data['clst2rep']
        clst2gene = {f"{name}{k}" : v for k,v in clst2gene.items()}
        gene_clusters2rep = {f"{name}{k}" : v for k,v in gene_clusters2rep.items()}
        if "get_ppangolin_partition" in kwargs and kwargs['get_ppangolin_partition']:
            if not get_quiet():
                print("Remove superpanf from ppangolin object")
            pang.remove_superpang()
            if not get_quiet():
                print("Running ppangolin partitioning")
            pang.run_just_partitioning()
            new_out = pang.parse_output()
            rep2id = { v : k for k,v in gene_clusters2rep.items()}
            new_id2id = { k : rep2id[v] for k,v in new_out['clst2rep'].items()}
            ppang_part = { new_id2id[k]  : v  for k,v in new_out['clst2partition'].items()}
            classes = {'C' : 'cloud', 'P' : 'persistent', 'S' : 'shell' }
            ppang_part = { k : classes[v] for k,v in ppang_part.items()}

    elif method == "mmseqsCluster" :
        covmode = kwargs.get("mmseqsCluster_covmode",0)
        cov = kwargs.get("mmseqsCluster_cov",0.80)
        seqid = kwargs.get("mmseqsCluster_seqid",0.80)
        f"coverage = 80% with cov-mode = 0, minimal amino acid sequence identity = 80% and cluster-mode = 0"

        if not shutil.which("mmseqs"):
            print("You need mmseqs2 to run the mmeseqs2 gene-clustering, either install it or run mOTUpan with an other gene-clustering or your own traits", file = sys.stderr)
            sys.exit(-1)

        print("running mmseqs easy-cluster with params --min-seq-id {seqid} --cov-mode {covmode} -c {cov} in {}".format(temp_folder, covmode = covmode, cov = cov, seqid=seqid), file = sys.stderr)
        mmseqs_dat = pjoin(temp_folder, "mmseqs_")
        logs = "2> /dev/null > /dev/null" if mOTUlizer._quiet_ else ""
        os.system("mmseqs easy-cluster --threads {threads}  --min-seq-id {seqid} --cov-mode {covmode} -c {cov} {faas} {out} {tmp}  {logs}".format(covmode = covmode, cov = cov, seqid=seqid, faas = all_faas_file, out = mmseqs_dat, tmp = temp_folder, threads = threads, logs = logs))

        with open(mmseqs_dat + "_cluster.tsv") as handle:
            if precluster:
                recs = {g : l[:-1].split()[0]  for l in handle for g in clusters[l[:-1].split()[1]]}
            else :
                recs = {l[:-1].split()[1] : l[:-1].split()[0]  for l in handle}

        #pretty formating names
        fill = len(str(len(set(recs.values()))))

        rep2clust = {k : name + str(i).zfill(fill) for i,k in enumerate(set(recs.values()))}
        gene_clusters2rep = {v: k for k,v in rep2clust.items()}
        if not mOTUlizer._quiet_:
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

    gc2genome = {clust : {prot2genome[g] for g in genes} for clust, genes in clst2gene.items()}

    clusts = [GeneCluster(representative = gene_clusters2rep[clust], annotations = {"ppanggolin_partition" : ppang_part.get(clust, "NA")}, name = clust, bunched = True, db = db) for clust, genes in clst2gene.items()]

    feat_tuples = [(c.name, gene) for c in clusts for g in clst2gene[c.name]]
    genome_tuples = [(c.name, g.name) for c in clusts for g in gc2genome[c.name]]

    db.add_feature_tuple2gc(feat_tuples, bunched = True)
    db.add_genome_tuple2gc(genome_tuples, bunched = True)

    db.commit()
    if "get_ppangolin_partition" in kwargs and kwargs['get_ppangolin_partition']:
        clusts = { 'clusts' : clusts, 'ppangolin_partitioning' : ppang_part}
    return clusts


class GeneCluster():
    def __eq__(self, gc):
        return self.name == gc.name

    def __hash__(self):
        return hash(self.name)

    def __repr__(self):
        return f"< GeneCluster '{self.name} 'of {0 if not self.features else len(self.features)} genes from {len(self.genomes)} genomes>"

    def __init__(self, name = None, features = None,  representative = None, annotations = None, genomes = None, bunched = False, db = None):

        if not db:
            if not SeqDb.seq_db:
                raise DataBaseNotInitialisedError(f"The database has not been initialised")
            self.db = SeqDb.get_global_db()
        else :
            self.db = db
        """write chekcs here"""


        if not name:
            self.name = "GeneCluster_" + random_name()
        else :
            self.name = name

        if name and self.db.has_gc(name):
            self._data = self.db.get_gc(name)
            self.features = self._data['features']
            self.annotations = self._data['annotations']
            self._representative_id = self._data['representative']
        else :
            self.db.add_gene_cluster(name = name, features = features, representative = representative, genomes = genomes, annotations = annotations, bunched = bunched)
            self.features = features
            self.annotations = annotations
            self._representative_id = representative

        self._mutation_dict = None

    def __len__(self):
        return len(self.features)

    @property
    def representative(self):
        if not hasattr(self, "_representative"):
            dd = self.db.get_features_seq([self._representative_id])
            self._representative = dd[int(self._representative_id)]
        return self._representative

    @property
    def feature2genome(self):
        if not hasattr(self, "_feature2genome"):
            self._feature2genome = {f : self.db.get_genomesFromFeat(f) for f in self.features}
        return self._feature2genome

    @property
    def genome2feature(self):
        if not hasattr(self, "_genome2feature"):
            self._genome2feature = {f : [] for f in self.genomes}
            for k,v in self.feature2genome.items():
                self._genome2feature[v] += [k]
        return self._genome2feature

    @property
    def genomes(self):
        return self.db.get_genome_from_GC(self.name)


    @property
    def genes(self):
        return self.db.get_genes_from_GC(self.name)


    def get_seqs(self):
        seqs = dict()
        for gene, genomes in self.gene2genome.items():
            seqs[gene] = genomes[0].get_aa(gene)
        return seqs

    def get_alignment(self, method = "muscle"):
        if len(self) < 2 :
            raise CantRunError("Making alignments for clusters with only one gene makes little sense...")
        if self.storage:
            out_file = pjoin(self.storage, "alignment.fasta")
        else :
            out_file = None
        if method == "muscle":
            tool = Muscle(seqsdict_or_infile = self.get_seqs(), out_file = out_file, quiet = True)
        else :
            raise CantMethodError(f"The alignment method ({method}) you asked for is not implemented yet.")

        if out_file is None or (out_file and not os.path.exists(out_file)):
            tool.run_command()
        elif not self.quiet:
            print("Alignment loaded from file", file = sys.stderr)

        ali = tool.parse_output()
        return ali

    def alignment_with_codons(self, write = False):
        codoned = {gid : dict() for gid in self.genes}
        if write and not self.storage:
            raise FileError("You want to write the codon alignment but you have no storage folder for your gene-cluster")
        for gid, allied in self.get_alignment().items():
            seq = self.gene2genome[gid][0].get_gene(gid)
            # a little dance to replace trailing and leadin "-" by X, I do not want to count them in the ANI also, add a stop
            allied = allied.rstrip("*")
            cleaned = allied.rstrip("-")
            trailing_gap = (len(allied) - len(cleaned))
            paded = cleaned + "X" * trailing_gap
            cleaned = paded.lstrip("-")
            leading_gap = (len(paded) - len(cleaned))
            allied = "X" * leading_gap  + cleaned + "*"
            codons = re.findall('...', str(seq))
#            assert len(allied.ungap("-")[:-1].strip("X")) == len(codons[:-1])
            codons = ["XXX"]  * leading_gap + codons[:-1] + ["XXX"] * trailing_gap + [codons[-1]]
            counter = 0
            for aa in allied:
                if aa == "-":
                    codons = codons[0:counter] + ["---"] + codons[counter:]
                counter += 1
            assert len(codons) == len(allied)
            codoned[gid]['AAs'] = allied
            codoned[gid]['codons'] = codons
        if write :
            SeqIO.write([SeqRecord(id = gid, description = "", seq = Seq("".join(v['codons']))) for gid, v in codoned.items()], pjoin(self.storage, "alignment_codons.fasta"), "fasta")
        return codoned

    @classmethod
    def _compare_two_seqs(cls,ali1, ali2):
        aa1, codon1 = (ali1['AAs'],ali1['codons'])
        aa2, codon2 = (ali2['AAs'],ali2['codons'])
        synonymous_muts = 0
        non_synonymous_muts = 0
        counted_bases = 0
        codon_diff = lambda a,b : sum([m != n for m,n in zip(a,b)])
        transitions = ""
        for i in range(len(aa1)):
            if codon1[i] != "XXX" or codon2[i] != "XXX":
                if codon1[i] != "---" and codon2[i] != "---":
                    counted_bases += 3
                    if codon1[i] != codon2[i]:
                        if aa1[i] == aa2[i]:
                            synonymous_muts += codon_diff(codon1[i],codon2[i])
                        else :
                            non_synonymous_muts += codon_diff(codon1[i],codon2[i])
                            transitions += "".join(sorted([aa1[i], aa2[i]]))
        return {'synonymous_muts' : synonymous_muts, 'non_synonymous_muts' : non_synonymous_muts, 'counted_bases' : counted_bases, 'transitions' : transitions}


    def within_cluster_mutations(self, write = False):
        if not self._mutation_dict:
            alignment = self.alignment_with_codons(write = write)
            ali_len = len(next(iter(alignment.values()))['AAs'])
            assert all([len(v['AAs']) == ali_len for v in alignment.values()])
            pairs = { frozenset(k) for k in itertools.permutations(alignment.keys(), 2)}
            pairs = [tuple(p) for p in pairs]
            self._mutation_dict = {p : GeneCluster._compare_two_seqs(alignment[p[0]], alignment[p[1]]) for p in pairs}
        return self._mutation_dict

    def synonymous_rate(self):
        return sum([v['synonymous_muts'] for v in self.within_cluster_mutations().values()])/sum([v['counted_bases'] for v in self.within_cluster_mutations().values()])

    def non_synonymous_rate(self):
        return sum([v['non_synonymous_muts'] for v in self.within_cluster_mutations().values()])/sum([v['counted_bases'] for v in self.within_cluster_mutations().values()])

    def tree(self, seq_type = "nucl", method = "fasttree"):
        ali = self.alignment_with_codons()
        nucl_seq = lambda data : Seq("".join([ codon.replace("X","-") for codon in data['codons'] ]))
        aa_seq = lambda data : Seq("".join([ "-" if aa == "X" else aa for aa in str(data['AAs']) ]))

        if self.storage:
            out_file = pjoin(self.storage, "tree.nwk")
        else :
            out_file = None
        if method == "fasttree":
            if seq_type == "nucl":
                seqs = {gid :  nucl_seq(data) for gid, data in ali.items()}
            elif seq_type == "aas":
                seqs = {gid :  aa_seq(data) for gid, data in ali.items()}
            else :
                raise CantMethodError(f"You can only align 'nucl' or 'aas', you tried {seq_type}")
            tool = FastTree(seqsdict_or_infile = seqs, out_file = out_file, quiet = self.quiet)
        else :
            raise CantMethodError(f"The treemaking method ({method}) you asked for is not implemented yet.")

        if out_file and not os.path.exists(out_file):
            tool.run_command()
        elif not self.quiet:
            print("Alignment loaded from file", file = sys.stderr)

        tree = tool.parse_output()
        return tree
