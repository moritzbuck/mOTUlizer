import sys
import os

import subprocess
from subprocess import call
import multiprocessing

import tempfile
import json

from random import shuffle, choice, choices
from math import log10

from mOTUlizer import __version__
from mOTUlizer.classes.GFF import GFF
from mOTUlizer.classes.MetaBin import MetaBin
from mOTUlizer.classes.GeneClusters import *
from mOTUlizer.config import *

mean = lambda x : sum(x)/len(x)

class mOTU:
    def __len__(self):
        return len(self.members)

    def __repr__(self):
        return "< {tax} mOTU {name}, of {len} members >".format(name = self.name, len = len(self), tax =  None ) #self.consensus_tax()[0].split(";")[-1])

    def __init__(self, genomes, name = "mOTU", gene_clusters = None, make_gene_clustering = False, compute_core = False, quiet = False, storage = None, **kwargs):
        self.can_haz_gene_clusters = True
        self.quiet = quiet
        if not self.quiet:
            print("Initializing mOTU", file = sys.stderr)
        self._genecluster_poolsize = None
        self.name = name
        self.likelies = None
        self.mock = []
        self.members = genomes
        self.gene_clusters = gene_clusters
        self.core = None
        self.method = None
        self.storage = storage
        if gene_clusters:
            if not all([g in self.gene_clusters.genomes for g in self]):
                GenomeMismatchError("the genomes in your genome clustering are not the same as the ones of your mOTU")

        self.genome_id2genome = {g.name : g for g in self}
        self._core_computed = False
        if  1 < sum([g.original_complet is None for g in self]) < len(self)-1 :
            self.merens_trick()
        if any([not g.original_complet for g in self]):
            self.estimate_complete_from_length()

        if make_gene_clustering:
            if self.gene_clusters and not kwargs.get('force', False):
                raise CantGeneClusterError("You want to recompute the genome clustering but you already have one\n, if you sure you want to, pass 'True' to the mOTU contructor")
            precluster = kwargs.get('precluster', False)
            threads = kwargs.get('threads', multiprocessing.cpu_count())
            self.gene_clusters = GeneClusters.compute_GeneClusters([self], name = name, precluster = precluster, threads = threads, storage = None if not self.storage else pjoin(self.storage, "gene_clusters"))

        if compute_core:
            max_it = kwargs.get('motupan_maxit', 100)
            method = kwargs.get('motupan_method', 'motupan_v0_3_2')
            self.compute_core(method, max_it)

    def load_gene_clusters(self, file, force = False):
        self.gene_clusters = GeneClusters.load(file, self, force = force, storage = None if not self.storage else pjoin(self.storage, "gene_clusters"))

    def export_gene_clusters(self, file, force = False):
        self.gene_clusters.export_gene_clusters(file = file , force = force)


    def compute_core(self, method = "motupan_v0_3_2", max_it = 100):
        self.method = method
        if self.method == "motupan_v0_3_2":
            likelies = self._core_likelyhood(max_it = max_it)
            self._core_computed = True
        self.likelies = likelies

    def __getitem__(self, i):
        if type(i) == int and i < len(self):
            return self.members[i]
        elif i in self.genome_id2genome:
            return self.genome_id2genome[i]
        else :
            raise KeyError(str(i) + "is not a valid entry")

    def get(self, key, default = None):
        if key in self.genome_id2genome:
            return self.genome_id2genome[key]
        else :
            return default

    def roc_values(self, boots):
        if not self._core_computed:
            raise CoreNotComputedError("Self explanatory error")
        if boots > 0 or len(self.mock) >0 :
            from mOTUlizer.classes.MockData import MockmOTU
            mean = lambda data: float(sum(data)/len(data))
            variance 	= lambda data, avg: sum([x**2 for x in [i-avg for i in data]])/float(len(data))
            std_dev = lambda data: variance(data, mean(data))**0.5

            while len(self.mock) < boots:
                print("Running bootstrap {}/{}".format(len(self.mock)+1, boots), file = sys.stderr)
                completnesses = {"Genome_{}".format(i) : c.new_completness for i,c in enumerate(self)}

                accessory = sorted([v for g,v in  self.gene_clusters.get_genecluster_counts().items() if g not in self.core])
                missing = int(sum(accessory)*(1-mean(list(completnesses.values()))/100))
                if len(accessory) > 0:
                    addeds = choices(list(range(len(accessory))), weights = accessory, k = missing)
                    for k in addeds :
                        accessory[k] += 1

                self.mock += [MockmOTU(self.name + "_mock", len(self.core), len(self), lambda g : completnesses[g], accessory = accessory, method = self.method)]
            return { 'mean_recall' : mean([m.recall for m in self.mock]),
                     'sd_recall' : std_dev([m.recall for m in self.mock]),
                     'mean_fpr' : mean([m.fpr for m in self.mock]),
                     'sd_fpr' : std_dev([m.fpr for m in self.mock]),
                     'mean_lowest_false' : mean([m.lowest_false for m in self.mock]),
                     'sd_lowest_false' : std_dev([m.lowest_false for m in self.mock]),
                     'nb_bootstraps' : len(self.mock)
                     }
        else:
            return { 'mean_recall' : "NA",
                     'sd_recall' : "NA",
                     'mean_fpr' : "NA",
                     'sd_fpr' : "NA",
                     'mean_lowest_false' : "NA",
                     'sd_lowest_false' : "NA",
                     'nb_bootstraps' : 0}


    def avg_gene_clusters_content(self):
        return sum([len(m.gene_clusters) for m in self.members])/len(self)

    def consensus_tax(self):
        taxos = [bin.taxonomy for bin in self.members]
        temp = {}
        for i in range(0,7):
            leveled = [t[i] for t in taxos]
            counts = {l : leveled.count(l) for l in set(leveled)}
            win = max(counts.items(), key = lambda x : x[1])
            temp[i] = win
        output = ["", []]
        for i in range(7):
            if temp[i][0] == "":
                break
            output[0] += ";" + temp[i][0]
            output[1] += [temp[i][1] / len(self)]
        output[0] = output[0][1:]
        return tuple(output)

    def overlap_matrix(self):
        if not hasattr(self, 'overlap_dict'):
            self.overlap_dict = {(i,j) : len(i.overlap(j))/len(i.gene_clusters) for i in self.members for j in self.members if i != j}
        return self.overlap_dict

    def mean_overlap(self) :
        return mean(list(self.overlap_matrix().values()))

    def get_stats(self):
        out = {self.name : {
                            "nb_genomes" : len(self),
                            }
              }
        if self.gene_clusters:
            out[self.name]["gene_clusters"] = {g.name : g.genes for g in self.gene_clusters}
        if self._core_computed:
            out[self.name].update({
            "core" : [g.name for g in self.core],
            "accessory_genome" :  [g.name for g in self.gene_clusters if g not in self.core],
            "singleton_gene_clusters" : [g.name for g in self.gene_clusters if len(g.get_genomes()) == 1],
            "genomes" : [v.name for v in self],
            "likelies" : {g.name : self.likelies[g] for g in self.gene_clusters}
            })
        if all([g.nucleotide_file for g in self]):
            out[self.name].update({
                "mean_ANI" : self.get_mean_ani(),
                "ANIs" : self.get_anis()
                })
        return out

    def get_representative(self, method = "complex", max_redund = 5, min_complete = 95):
            tt = [v.get_data() for v in self]
            data = { t['name'] : t for t in tt}
            if all([v['original_contamin'] > max_redund for v in data.values()]):
                return max( [ (k , v['original_complet']) for k,v in data.items()], key = lambda x: x[1])[0]
            data = {k : v for k,v in data.items() if v['original_contamin'] < max_redund}
            if any([v['original_complet'] > min_complete for v in data.values()])  :
                data = {k : v for k,v in data.items() if v['original_complet'] > min_complete}
                best_redund = min(data.items(), key = lambda x : x[1]['original_contamin'])[1]['original_contamin']
                return max( [ (k , v['original_complet']) for k,v in data.items() if v['original_contamin'] == best_redund], key = lambda x: x[1])[0]
            else:
                return max( [ (k , v['original_complet']) for k,v in data.items()], key = lambda x: x[1])[0]

    # def get_representative(tt, max_redund = 5, min_complete = 95):
    #         data = { t['name'] : t for t in tt}
    #
    #         data = {k : v for k,v in data.items() if v['original_contamin'] < max_redund}
    #         if any([v['original_complet'] > min_complete for v in data.values()])  :
    #             data = {k : v for k,v in data.items() if v['original_complet'] > min_complete}
    #             best_redund = min(data.items(), key = lambda x : x[1]['original_contamin'])[1]['original_contamin']
    #             return min( [ (k , v['original_contamin']) for k,v in data.items() if v['original_contamin'] == best_redund], key = lambda x: x[1])[0]
    #         elif len(data) >0 :
    #             return max( [ (k , v['original_complet']) for k,v in data.items()], key = lambda x: x[1])[0]
    #         else :
    #             return None

    def get_mean_ani(self):
        dist_dict = self.get_anis()
        dists = [dist_dict.get((a.name,b.name)) for a in self for b in self if a != b ]
        missing_edges = sum([d is None for d in dists])
        found_edges = [d is None for d in dists]

        return {'mean_ANI' : sum([d for d in dists if d])/len(found_edges) if len(found_edges) > 0 else None, 'missing_edges' : missing_edges, 'total_edges' : len(found_edges) + missing_edges}

    def _core_likelyhood(self, max_it = 20, likeli_cutof = 0 ):
        likelies = {gene_cluster : self._core_likely(gene_cluster) for gene_cluster in self.gene_clusters}
        core = set([c for c, v in likelies.items() if v > likeli_cutof])
        core_len = len(core)
        i = 1
        if not self.quiet:
            print("iteration 1 : ", core_len, "sum_abs_LLHR:" , sum([l if l > 0 else -l for l in likelies.values()]), file = sys.stderr)
        for mag in self:
            if len(core) > 0:
                mag.new_completness = 100*len(mag.get_clust_set().intersection(core))/len(core)
            else :
                mag.new_completness = 0
            mag.new_completness = mag.new_completness if mag.new_completness < 99.9 else 99.9
            mag.new_completness = mag.new_completness if mag.new_completness > 0 else 0.01
        for i in range(2,max_it):
            likelies = {gene_cluster : self._core_likely(gene_cluster, complet = "new", core_size = core_len) for gene_cluster in self.gene_clusters}
            old_core = core
            core = set([c for c, v in likelies.items() if v > likeli_cutof])
            new_core_len = len(core)
            for mag in self:
                if len(core) > 0:
                    mag.new_completness = 100*len(mag.get_clust_set().intersection(core))/len(core)
                else :
                    mag.new_completness = 0
                mag.new_completness = mag.new_completness if mag.new_completness < 99.9 else 99.9
                mag.new_completness = mag.new_completness if mag.new_completness > 0 else 0.01
            if not self.quiet:
                print("iteration",i, ": ", new_core_len, "sum_abs_LLHR:" , sum([l if l > 0 else -l for l in likelies.values()]), file = sys.stderr)
            if core == old_core:
               break
            else :
                core_len = new_core_len

        if not self.quiet:
            pp =  "\nYour {name}-run for {nb_mags} genomes (with mean initial completeness {mean_start:.2f}) resulted\n"
            pp += "in a core of {core_len} traits with a total sum of loglikelihood-ratios {llhr:.2f} and a corrected \n"
            pp += "mean completness of {mean_new:.2f}, resulting to a estimated mean traits per genome count of {trait_count:.2f}\n"
            pp = pp.format(name = self.name, nb_mags = len(self), core_len = core_len, mean_start = mean([b.original_complet for b in self]),
                        mean_new =  mean([b.new_completness for b in self]), llhr =  sum([l if l > 0 else -l for l in likelies.values()]),
                        trait_count = mean([100*len(b.gene_clusters)/b.new_completness for b in self]))
            print(pp, file = sys.stderr)
        self.iterations = i -1
        self.core = core
        self._core_computed = True
        return likelies

    def _core_prob(self, gene_cluster, complet = "checkm"):
        comp = lambda mag : (mag.original_complet if complet =="checkm" else mag.new_completness)/100
        presence = [log10(comp(mag)) for mag in self if gene_cluster in mag.get_clust_set()]
        abscence = [log10(1 - comp(mag)) for mag in self if gene_cluster not in mag.get_clust_set()]
        return sum(presence + abscence)

    def get_genecluster_poolsize(self):
        if not self._genecluster_poolsize:
            self._genecluster_poolsize = sum([len(c.get_genomes()) for c in  self.gene_clusters])
        return self._genecluster_poolsize

    def _pange_prob(self, gene_cluster, core_size, complet = "checkm"):
#        pool_size = sum(self.gene_clustersCounts.values())
        pool_size = self.get_genecluster_poolsize()
        comp = lambda mag : (mag.original_complet if complet =="checkm" else mag.new_completness)/100
        #presence = [1 - (1-self.gene_clustersCounts[gene_clusters]/pool_size)**(len(mag.gene_clusterss)-(core_size*comp(mag))) for mag in self if gene_clusters in mag.gene_clusterss]
        #abscence = [ (1-self.gene_clustersCounts[gene_clusters]/pool_size)**(len(mag.gene_clusterss)-(core_size*comp(mag))) for mag in self if gene_clusters not in mag.gene_clusterss]

#        mag_prob = {mag : ( 1-1/pool_size )**len(mag.gene_clusterss) for mag in self}
#        mag_prob = {mag : ( 1-1/pool_size )**(len(mag.gene_clusterss)-(core_size*comp(mag))) for mag in self}
#        mag_prob = {mag : ( 1-self.gene_clustersCounts[gene_clusters]/pool_size )**(len(mag.gene_clusterss)-(core_size*comp(mag))) for mag in self}

#        mag_prob = {mag : ( 1-self.gene_clustersCounts[gene_clusters]/pool_size )**(len(mag.gene_clusterss)-(core_size*comp(mag))) for mag in self}
        mag_prob = {mag : ( 1-len(gene_cluster.get_genomes())/pool_size)**len(mag.gene_clusters) for mag in self}

        presence = [ log10(1 -   mag_prob[mag]) if mag_prob[mag] < 1 else MIN_PROB                for mag in self if mag.name in gene_cluster.get_genomes()]
        abscence = [ log10(      mag_prob[mag]) if mag_prob[mag] > 0 else log10(1-(10**MIN_PROB)) for mag in self if mag.name not in gene_cluster.get_genomes()]

        #abscence = [ 1-self.gene_clustersCounts[gene_clusters]/len(self)*comp(mag) for mag in self if gene_clusters not in mag.gene_clusterss]
        #presence = [ self.gene_clustersCounts[gene_clusters]/len(self)*comp(mag) for mag in self if gene_clusters not in mag.gene_clusterss]

        return sum(presence + abscence)

    def _core_likely(self, gene_cluster, complet = "checkm", core_size = 0):
        pange_prob = self._pange_prob(gene_cluster, core_size, complet)
        return self._core_prob(gene_cluster, complet) - pange_prob

    def nb_gene_clusters(self):
        return mean([b.estimate_nb_gene_clusters() for b in self if b.new_completness > 40 ])

    def get_pangenome_size(self, singletons = False):
        return len([k for k,v in self.gene_clustersCounts.items() if k not in self.core and v > (0 if singletons else 1)])

    def pretty_pan_table(self):

        out_dict = {}
        stats = self.get_stats()
        stats = list(stats.values())[0]
        stats.update(self.roc_values(boots = len(self.mock)))
        for g in self.gene_clusters:
            k = g.name
            out_dict[k] = {}
            out_dict[k]['type'] = 'core' if k in stats['core'] else 'accessory'
            out_dict[k]['genome_occurences'] = 0
            out_dict[k]['log_likelihood_to_be_core'] = stats['likelies'][k]
            out_dict[k]['genomes'] = []
            out_dict[k]['trait_name'] = k
            out_dict[k]['genes'] = g.genes if g.genes else "NA"
            out_dict[k]['representative'] = g.representative if g.representative else "NA"
            out_dict[k]['genomes'] += [gg.name for gg in g.get_genomes()]
            out_dict[k]['genome_occurences'] = len(g.get_genomes())

        for k,v in out_dict.items():
            v['mean_copy_per_genome'] = "NA" if not v['genes'] else len(v['genes'])/len(v['genomes'])
            v['genes'] = ";".join(v['genes'])
            v['genomes'] = ";".join(v['genomes'])

        header = ['trait_name','type', 'genome_occurences', 'log_likelihood_to_be_core', 'mean_copy_per_genome','genomes', 'genes']
        genome_line = "genomes=" + ";".join(["{}:prior_complete={}:posterior_complete={}".format(k.name, k.original_complet, k.new_completness) for k in self])
        mean = lambda l : sum([ll for ll in l])/len(l)

        if stats['mean_recall'] != "NA":
            bootsy = """
#bootstrapped_mean_false_positive_rate={fpr:.2f};bootstrapped_sd_false_positive_rate={sd_fpr:.2f}
#bootstrapped_mean_recall={recall:.2f};bootstrapped_sd_recall={sd_recall:.2f}
#bootstrapped_mean_lowest_false_positive={lowest:.2f};bootstrapped_sd_lowest_false_positive={sd_lowest:.2f}
#bootstrapped_nb_reps={boots}
#""".format( boots=len(self.mock),
            fpr=stats['mean_fpr'],
            recall = stats['mean_recall'], lowest = stats['mean_lowest_false'],sd_fpr=stats['sd_fpr'],
            sd_recall = stats['sd_recall'], sd_lowest = stats['sd_lowest_false'])
        else :
            bootsy=""


        outformat ="""#mOTUlizer:mOTUpan:{version}
#run_name={name}
#
#genome_count={nb}
#core_length={core_len}
#mean_prior_completeness={prior_complete:.2f}
#mean_posterior_completeness={post_complete:.2f}
#sum_abs_loglikelihood_ratios={SALLHR:.2f}
#mean_est_genome_size={size:.2f};traits_per_genome
#{genomes}
#{boostrap}
{header}
{data}
"""

        return outformat.format(version = __version__ , nb = stats['nb_genomes'],
            name = self.name.strip("_"),
            core_len = len(self.core),
            genomes=genome_line,
            prior_complete=mean([b.original_complet for b in self]),
            post_complete=mean([b.new_completness for b in self]),
            SALLHR= -1 if not self.likelies else sum([l if l > 0 else -l for l in self.likelies.values()]),
            size=mean([100*len(b.gene_clusters)/b.new_completness for b in self]),
            boostrap=bootsy,
            header = "\t".join(header), data = "\n".join(["\t".join([str(v[hh]) for hh in header]) for v in out_dict.values()]))

    def load_anis(self, ani_dict_or_file, force = False):
        if type(ani_dict_or_file) == dict:
            ani_dict = ani_dict_or_file
            if not self.quiet:
                print(f"checking if the {len(ani_dict)} entries of your imported ani correspond to genomes")
            if not force:
                for k,v in ani_dict.items():
                    if "ani" not in v:
                        raise ValueError("there should be a dict with a keys 'ani' in there")
                    for kk in k:
                        if kk not in self.genome_id2genome:
                            raise GenomeIdError(f"Some genome IDs of your imported similarity file are not in the genome set,\n use 'force = True' if you just want to subset,\n the genome id that crashed it is {kk}")
            self.anis = {(self.genome_id2genome[k[0]], self.genome_id2genome[k[1]]) : v for k,v in ani_dict.items()}
        elif os.path.exists(ani_dict_or_file):
            self.anis = dict()
            with open(ani_dict_or_file) as handle:
                for l in handle:
                    if "query" not in l:
                        ll = l.split("\t")
                        if "." in ll[0] :
                            g1 = ".".join(os.path.basename(ll[0]).split(".")[:-1]) if any([ll[0].endswith(ext) for ext in FASTA_EXTS]) else ll[0]
                        else :
                            g1 = ll[0]
                        if "." in ll[1] :
                            g2 = ".".join(os.path.basename(ll[1]).split(".")[:-1]) if any([ll[1].endswith(ext) for ext in FASTA_EXTS]) else ll[1]
                        else :
                            g2 = ll[1]
                        dist = float(ll[2])
                        if g1 not in self.genome_id2genome and not force:
                            raise GenomeIdError(f"Some genome IDs of your imported similarity file are not in the genome set,\n use 'force = True' if you just want to subset,\n the genome id that crashed it is {g1}")
                        if g1 not in self.genome_id2genome and not force:
                            raise GenomeIdError(f"Some genome IDs of your imported similarity file are not in the genome set,\n use 'force = True' if you just want to subset,\n the genome id that crashed it is {g2}")
                        self.anis[(self.genome_id2genome[g1],self.genome_id2genome[g2])] = {'ani' : dist, 'query_chunks' : None if len(ll) < 3 else float(ll[3]) , 'reference_chunks' : None if len(ll) < 4 else float(ll[4]) }

    def export_anis(self, file_path):
        head = ["query" , "reference" , "ani" , "query_chunks" , "reference_chunks"]
        with open(file_path,"w") as handle:
            handle.writelines(["\t".join(head) + "\n"] + [f"{k[0].name}\t{k[1].name}\t{v['ani']}\t{v['query_chunks']}\t{v['reference_chunks']}\n" for k,v in self.get_anis().items() ] )
    def get_anis(self, method = "fastANI", block_size = 500, threads=1):
        if not hasattr(self, 'anis'):
            if method == "fastANI":
                if not shutil.which('fastANI'):
                    print("You need fastANI if you do not provide a file with pairwise similarities, either install it or provide pairwise similarities (see doc...)", file = sys.stderr)
                    sys.exit(-1)
                fastani_file = tempfile.NamedTemporaryFile().name

                mags = [b.nucleotide_file for b in self]
                if any([m is None for m in mags]):
                    raise CantNucleotideError("All the genomes you want to have ANI for have to have a a nucleotide fasta file")

                mag_blocks = [mags[i:(i+block_size)] for i in list(range(0,len(mags), block_size))]

                if len(mag_blocks) > 1:
                    print("You have more then {bsize} bins, so we will run fastANI in blocks, if it crashes due to memory, make smaller blocks".format(bsize = block_size), file=sys.stderr)

                with open(fastani_file, "w") as handle:
                    handle.writelines(["query\tsubject\tani\tsize_q\tsize_s\n"])

                for i,bloc1 in enumerate(mag_blocks):
                    b1_tfile = tempfile.NamedTemporaryFile().name

                    with open(b1_tfile, "w") as handle:
                        handle.writelines([l +"\n" for l in bloc1])

                    for j,bloc2 in enumerate(mag_blocks):
                            print("doing bloc {i} and {j}".format(i = i, j=j), file = sys.stderr)
                            b2_tfile = tempfile.NamedTemporaryFile().name
                            with open(b2_tfile, "w") as handle:
                                handle.writelines([l  +"\n" for l in bloc2])

                            out_tfile = tempfile.NamedTemporaryFile().name
                            call("fastANI --ql {b1} --rl {b2} -o {out} -t {threads} 2> /dev/null".format(b1 = b1_tfile, b2 = b2_tfile, out = out_tfile, threads = threads), shell = True)
                            with open(out_tfile) as handle:
                                new_dat = ["\t".join([ll for ll in l.split()]) +"\n" for l in handle.readlines()]
                            with open(fastani_file, "a") as handle:
                                handle.writelines(new_dat)

                            os.remove(out_tfile)
                            os.remove(b2_tfile)

                os.remove(b1_tfile)
                with open(fastani_file) as handle:
                    handle.readline()
                    out_dists = {(os.path.basename(l.split()[0]), os.path.basename(l.strip().split()[1])) : {'ani' : float(l.split()[2]), 'query_chunks' : float(l.split()[3]), 'reference_chunks' : float(l.split()[2]) } for l in handle}
                    out_dists = {( ".".join(k[0].split(".")[:-1]) if any([k[0].endswith(ext) for ext in FASTA_EXTS]) else k[0],
                                   ".".join(k[1].split(".")[:-1]) if any([k[1].endswith(ext) for ext in FASTA_EXTS]) else k[1] ): v
                                   for k,v in out_dists.items() }
                os.remove(fastani_file)
            else :
                print("No other method for ani computation implemented yet")
                sys.exit()
            self.anis = {(self.genome_id2genome[k[0]], self.genome_id2genome[k[1]]) : v for k,v in out_dists.items()}
        return self.anis

    def cluster_MetaBins(self, ani_cutoff = 95, prefix = "mOTU_", mag_complete = 40, mag_contamin = 5, sub_complete = 0, sub_contamin = 100, threads = 1):
        import igraph

        dist_dict = self.get_anis(threads = threads)

        if not self.quiet:
            print("seeding bin-graph", file = sys.stderr )

        all_bins = {a.name : a for a in self}

        good_mag = lambda b : b.original_complet > mag_complete and b.original_contamin < mag_contamin
        decent_sub = lambda b : b.original_complet > sub_complete and b.original_contamin < sub_contamin and not good_mag(b)
        good_pairs = [k for k,v  in dist_dict.items() if v['ani'] > ani_cutoff and dist_dict.get((k[1],k[0]), 0)['ani'] > ani_cutoff and good_mag(k[0]) and good_mag(k[1])]
        species_graph = igraph.Graph()
        vertexDeict = { v : i for i,v in enumerate(set([x for k in good_pairs for x in k]))}
        rev_vertexDeict = { v : i for i,v in vertexDeict.items()}
        species_graph.add_vertices(len(vertexDeict))
        species_graph.add_edges([(vertexDeict[k[0]], vertexDeict[k[1]]) for k in good_pairs])

        print("getting clusters", file = sys.stderr)

        genome_clusters = [[rev_vertexDeict[cc] for cc in c ] for c in species_graph.components(mode=igraph.STRONG)]

        mean = lambda l : sum([len(ll) for ll in l])/len(l)

        print("recruiting to graph of the", len(genome_clusters) ," mOTUs of mean length", mean(genome_clusters), file = sys.stderr)


        left_pairs = {k : v['ani'] for k, v in dist_dict.items() if v['ani'] > ani_cutoff and k[0] != k[1] and ((decent_sub(k[0]) and good_mag(k[1])) or (decent_sub(k[1]) and good_mag(k[0])))}
        print("looking for good_left pairs", file = sys.stderr)
#        print(left_pairs)

        subs = {l : (None,0) for ll in left_pairs.keys() for l in ll if not good_mag(l)}
#        print(subs)
        print("looking for best mOTU match", file = sys.stderr)
        for p,ani in left_pairs.items():
            if p[0] in subs and subs[p[0]][1] < ani:
                subs[p[0]] = (p[1], ani)
            if p[1] in subs and subs[p[1]][1] < ani:
                subs[p[1]] = (p[0], ani)

        genome_clusters = [set(gg) for gg in genome_clusters]

        print("append to the", len(genome_clusters) ,"mOTUs of mean length", mean(genome_clusters), file = sys.stderr)
        for k, v in subs.items():
            for g in genome_clusters:
                if v[0] in g :
                    g.add(k)

        genome_clusters = [list(gg) for gg in genome_clusters]

        print("processing the", len(genome_clusters) ,"mOTUs of mean length", mean(genome_clusters), file = sys.stderr)
        #print(genome_clusters)

        zeros = len(str(len(genome_clusters)))

        genome2clust = {gg : i for i, gs in enumerate(genome_clusters) for gg in gs}
        dd_dicts = [dict() for i in range(len(genome_clusters))]
        for k,v in dist_dict.items():
            c1 = genome2clust.get(k[0], "z1")
            c2 = genome2clust.get(k[1], "z2")
            if c1 == c2:
                dd_dicts[c1][(k[0].name, k[1].name)] = v

        motus = [ mOTU(genomes = gs, name = prefix + str(i).zfill(zeros), quiet = True) for i, gs in enumerate(genome_clusters)]
        for motu, dists in zip(motus, dd_dicts):
            motu.load_anis(dists)

        return motus

    def __iter__(self):
       return mOTUIterator(self)

    def items(self):
        return mOTUItemIterator(self)

class mOTUIterator:
    def __init__(self, motu):
        self._motu = motu
        self.index = -1

    def __next__(self):
        if self.index < len(self._motu)-1:
            self.index += 1
            return self._motu[self.index]
        raise StopIteration

class mOTUItemIterator:
    def __init__(self, motu):
        self._motu = motu
        self.index = -1

    def __iter__(self):
        return self

    def __next__(self):
        if self.index < len(self._motu)-1:
            self.index += 1
            return (self._motu[self.index].name, self._motu[self.index])
        raise StopIteration
