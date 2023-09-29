import sys
import os

import subprocess
from subprocess import call
import multiprocessing

import tempfile
import json

from random import shuffle, choice, choices
from math import log10
from statistics import mean, median

from mOTUlizer import __version__
import mOTUlizer
from mOTUlizer.classes.MetaBin import MetaBin
from mOTUlizer.classes.GeneClusters import *
from mOTUlizer.config import *
from mOTUlizer.db.SeqDb import SeqDb
from mOTUlizer import get_quiet, get_threads
from mOTUlizer.classes.tools.SuperPang import SuperPang
from mOTUlizer.classes.GeneClusters import compute_GeneClusters
import multiprocessing as mp

mean = lambda x : "NA"  if "NA" in x else sum(x)/len(x)
def _annot_single(metabin, temp_dir, method, tool_args):
    print(f"Doing {metabin}")
    source, features = MetaBin(metabin, temp_dir = temp_dir).annotate(method, threads=1, tool_args = tool_args)
    return (metabin, features, source)

def _roc_single(name, len_core, nb_genomes, completnesses, accessory, method):
    from mOTUlizer.classes.MockData import MockmOTU
    mock = MockmOTU(name, len_core, nb_genomes, lambda g : completnesses[g], accessory = accessory, method = method)
    return { 'recall' : mock.recall,
             'fpr' : mock.recall,
             'lowest_false' : mock.lowest_false
             }

class mOTU:
    def __eq__(self, motu):
        return self.name == motu.name

    def __len__(self):
        return len(self.members)

    def __repr__(self):
        return "< {tax} mOTU {name}, of {len} members >".format(name = self.name, len = len(self), tax =  None ) #self.consensus_tax()[0].split(";")[-1])

    @classmethod
    def from_datastructure(cls, datastructure):
        #make gene_clusters:

        #make genomes :
        pass

    def __init__(self, name, genomes = None, make_gene_clustering = False, compute_core = False, quiet = None, storage = None, db = None, **kwargs):
        if not db:
            self.db = SeqDb.get_global_db()
            if not SeqDb.seq_db:
                raise DataBaseNotInitialisedError("The database has not been initialised")
        else:
            self.db = db
        if self.db.has_mOTU(name):
            self.members = [MetaBin(g, db = self.db) for g in self.db.get_mOTU(name)]
        else :
            self.members = genomes
            self.db.add_mOTU(name, [g.name for g in genomes])

        self.can_haz_gene_clusters = True
        self.quiet = quiet
        if not mOTUlizer._quiet_:
            print("Initializing mOTU", file = sys.stderr)
        self._genecluster_poolsize = None
        self.name = name
        self.likelies = None
        self.mock = []
        self.core = self.db.get_core(self.name)
        if self.core :
            self.core = [GeneCluster(c, db = self.db) for c in self.core]
            self.likelies = self.db.get_likelies(self.name)

        self.method = None

        self.genome_id2genome = {g.name : g for g in self}

        self._core_computed = self.core is not None

        if make_gene_clustering:
            if self.gene_clusters and not kwargs.get('force', False):
                raise CantGeneClusterError("You want to recompute the genome clustering but you already have one\n, if you sure you want to, pass 'True' to the mOTU contructor")
            precluster = kwargs.get('precluster', False)
            threads = kwargs.get('threads', multiprocessing.cpu_count())
            compute_GeneClusters([self], name = name, precluster = precluster, threads = threads, db = self.db)

        if  1 < sum([g.completeness is None for g in self]) < len(self)-1 :
            self.merens_trick()
        if self.gene_clusters and any([not g.completeness for g in self]):
            self.estimate_complete_from_length()

        if compute_core:
            max_it = kwargs.get('motupan_maxit', 100)
            method = kwargs.get('motupan_method', 'motupan_v0_3_2')
            self.compute_core(method, max_it)

    @property
    def superpang_pangenome(self):
        if len(self) <2:
            raise CantPangenome("You need more than one genome to run SuperPang")
        pange_name = (self.name + "_SuperPang").replace(";","_")
        if not hasattr(self, "_superpang_pangenome"):
            if self.db.has_genome(pange_name):
                self._superpang_pangenome = MetaBin(pange_name)
            else :
                sp =  SuperPang(self, quiet = False, threads = get_threads())
                sp.run_command()
                seqs = sp.parse_output()
                seqs = [{'contig_name' : self.name + "_SuperPang_" + s_id.split("_length")[0], 'sequence' : seq, 'annotations' : {'core' : s_id.endswith("_core")}} for s_id, seq in seqs]
                self._superpang_pangenome = MetaBin(name = pange_name)
                self.db.add_contigs(pange_name, contigs = seqs)
                self.db.commit()
        return self._superpang_pangenome

    def minipipe(self):
        pange = self.superpang_pangenome
        if not self.db.has_features(pange.name):
            print("gene-call SuperPang")
            pange.annotate(commit = True, tool_args = {'annotate' : False})
        if not self.core:
            print("Computing panggolin GCs and partition")
            ppang_partitions = compute_GeneClusters(self.members + [pange], name = self.name.replace(";","_") + "_w_pange_gcs", method = "ppanggolin", get_ppangolin_partition = True)['ppangolin_partitioning']

            print("Running mOTUpan")
            self.compute_core()

        print("loading cores")
        superpang_core = {gc for gc in pange.gene_clusters if any([gc.db.get_features_data(feat)['contig_name'].endswith('-core') for feat in gc.genome2feature[pange.name]])}
        motupan_core = self.core
        ppang_core = {k for k in self.gene_clusters if "ppanggolin_partition" in k.annotations and k.annotations['ppanggolin_partition'] == "persistent" }
        return {'superpang_core' : superpang_core, 'motupan_core' : motupan_core, 'ppang_core': ppang_core}


    def estimate_complete_from_length(self):
        max_len = max([len(self.genome2gcs[g.name]) for g in self])
        for g in self:
            g.completeness = 100*len(self.genome2gcs[g.name])/max_len

    @property
    def genome2gcs(self):
        if not hasattr(self, "_genome2gcs"):
            if not mOTUlizer._quiet_:
                print(f"loading genome GCs for {self.name} for the first time")
            tt = self.db.get_mOTU_GCs(self.name)
            if len(tt) == 0:
                return None
            objs = {gg for g in tt.values() for gg in g}
            objs = {gg : GeneCluster(gg, db = self.db) for gg in tqdm(objs)}
            self._genome2gcs = {k :  {objs[vv] for vv in v} for k,v in tt.items()}
        return self._genome2gcs

    @property
    def pangenome(self):
        if not hasattr(self, "_pangenome"):
            self._pangenome = set.union(*self.genome2gcs.values())
        return self._pangenome


    @property
    def gene_clusters(self):
        if not hasattr(self, "_gene_clusters"):
            gcs = self.genome2gcs
            if gcs is None : 
                return gcs
            if len(gcs) > 0 : 
                gcs = set.union(*list(gcs.values()))
            if len(gcs) == 0 :
                return None
            self._gene_clusters = gcs
        return self._gene_clusters

    def stats(self):
        if not mOTUlizer.get_quiet() :
            print("Getting pangenome")
        pangenome = {gc.name for gcs in self.genome2gcs.values() for gc in gcs} 
        if not mOTUlizer.get_quiet() :
            print("Getting estimate genome lengths")
        est_len = [100*len(set(self.genome2gcs[g.name]))/g.completeness for g in self]
        if not mOTUlizer.get_quiet() :
            print("Getting completnesses and redundancies")
        cs = [g.completeness for g in self]
        rs = [g.redundancy for g in self]
        new_cs = [-1 if not self.core else g.motupan_completeness(self.core, set(self.genome2gcs[g.name])) for g in self]
        return { 'mean_completeness' : mean(cs) ,
                 'median_completeness' : median(cs) , 
                 'min_completeness'  : min(cs),
                 'max_completeness'  : max(cs),
                 'mean_redundancy' : mean(rs) ,
                 'median_redundancy' : median(rs) , 
                 'min_redundancy'  : min(rs),
                 'max_redundancy'  : max(rs),
                 'mean_new_completeness' : mean(new_cs) ,
                 'median_new_completeness' : median(new_cs) , 
                 'min_new_completeness'  : min(new_cs),
                 'max_new_completeness'  : max(new_cs),
                 'mean_est_len' : mean(est_len),
                 'median_est_len' : median(est_len),
                 'pangenome_len' : len(pangenome),
                 'core_len' : len(self.core),
                 'number_ags' : len(self)
        }

    def compute_core(self, method = "motupan_v0_3_2", max_it = 100):
        self.method = method
        if self.method == "motupan_v0_3_2":
            likelies = self._core_likelyhood(max_it = max_it)
            self._core_computed = True

        self.db.add_core(self.name, [c.name for c in self.core], {k.name : v for k,v in likelies.items()})
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

    def roc_values(self, boots, threads = 20):
        if not self._core_computed:
            raise CoreNotComputedError("Self explanatory error")

        if boots > 0 :
            mean = lambda data: float(sum(data)/len(data)) if 'NA' not in data else "NA"
            variance 	= lambda data, avg: sum([x**2 for x in [i-avg for i in data]])/float(len(data))  if 'NA' not in data else "NA"
            std_dev = lambda data: variance(data, mean(data))**0.5  if 'NA' not in data else "NA"

            to_do = [] 
            print(f"Running {len(to_do)} bootstrapped FPRs")
            pool = mp.Pool(threads)
            for x in range(boots):
                completnesses = {"Genome_{}".format(i) : c.motupan_completeness(self.core) for i,c in enumerate(self)}
                accessory = sorted([len(gc.genomes) for gc in  self.gene_clusters if gc not in self.core])
                missing = int(sum(accessory)*(1-mean(list(completnesses.values()))/100))
                if len(accessory) > 0:
                    addeds = choices(list(range(len(accessory))), weights = accessory, k = missing)
                    for k in addeds :
                        accessory[k] += 1
                to_do += [ (self.name, len(self.core), len(self), completnesses, accessory, self.method ) ]


            to_stat =  pool.starmap_async(_roc_single, to_do).get()
            return { 'mean_recall' : mean([m['recall'] for m in to_stat]),
             'sd_recall' : std_dev([m['recall'] for m in to_stat]),
             'mean_fpr' : mean([m['fpr'] for m in to_stat]),
             'sd_fpr' : std_dev([m['fpr'] for m in to_stat]),
             'mean_lowest_false' : mean([m['lowest_false'] for m in to_stat]),
             'sd_lowest_false' : std_dev([m['lowest_false'] for m in to_stat]),
             'nb_bootstraps' : len(to_stat)
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
        output[1] = ";".join([f"{o:02}" for o in output[1]])
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
            "singleton_gene_clusters" : [g.name for g in self.gene_clusters if len(g.genomes) == 1],
            "genomes" : [v.name for v in self],
            "likelies" : {g.name : self.likelies[g] for g in self.gene_clusters}
            })
        if all([g.nucleotide_file for g in self]):
            out[self.name].update({
                "mean_ANI" : self.get_mean_ani(),
                "ANIs" : self.get_anis()
                })
        return out

    def get_gc_proteoms(self):
        return {gc.name : gc.representative for gc in tqdm(self.gene_clusters)}
    
    def get_representative(self, method = "complex", max_redund = 5, min_complete = 95):
        if method == "complex":
            tt = [v.get_data() for v in self]
            data = { t['name'] : t for t in tt}
            if all([v['original_redundancy'] > max_redund for v in data.values()]):
                return max( [ (k , v['_original_complet']) for k,v in data.items()], key = lambda x: x[1])[0]
            data = {k : v for k,v in data.items() if v['original_redundancy'] < max_redund}
            if any([v['_original_complet'] > min_complete for v in data.values()])  :
                data = {k : v for k,v in data.items() if v['_original_complet'] > min_complete}
                best_redund = min(data.items(), key = lambda x : x[1]['original_redundancy'])[1]['original_redundancy']
                return max( [ (k , v['_original_complet']) for k,v in data.items() if v['original_redundancy'] == best_redund], key = lambda x: x[1])[0]
            else:
                return max( [ (k , v['_original_complet']) for k,v in data.items()], key = lambda x: x[1])[0]
        if method == "good_centroid":
            goods = [g for g in self if g.completeness > min_complete and g.redundancy < max_redund]
            anis = {g.name : [] for g in goods}

            for k,v in self.get_anis().items():
                for gg, nn in [k, tuple(reversed(k))]:
                    if gg in anis:
                        anis[gg] +=  [v['ani']*self[nn].completeness/100]
            return self[max(anis.items(), key = lambda p: sum(p[1]))[0]]

    # def get_representative(tt, max_redund = 5, min_complete = 95):
    #         data = { t['name'] : t for t in tt}
    #
    #         data = {k : v for k,v in data.items() if v['original_contamin'] < max_redund}
    #         if any([v['_original_complet'] > min_complete for v in data.values()])  :
    #             data = {k : v for k,v in data.items() if v['_original_complet'] > min_complete}
    #             best_redund = min(data.items(), key = lambda x : x[1]['original_contamin'])[1]['original_contamin']
    #             return min( [ (k , v['original_contamin']) for k,v in data.items() if v['original_contamin'] == best_redund], key = lambda x: x[1])[0]
    #         elif len(data) >0 :
    #             return max( [ (k , v['_original_complet']) for k,v in data.items()], key = lambda x: x[1])[0]
    #         else :
    #             return None

    def get_mean_ani(self):
        dist_dict = self.get_anis()
        dists = [dist_dict.get((a.name,b.name))['ani'] for a in self for b in self if a != b ]
        missing_edges = sum([d is None for d in dists])
        found_edges = sum([not d is None for d in dists])

        return {'mean_ANI' : sum([d for d in dists if d])/found_edges if found_edges > 0 else None, 'missing_edges' : missing_edges, 'total_edges' : found_edges + missing_edges}

    def _core_likelyhood(self, max_it = 20, likeli_cutof = 0 ):


        likelies = {gene_cluster : self._core_likely(gene_cluster) for gene_cluster in self.gene_clusters}
        core = set([c for c, v in likelies.items() if v > likeli_cutof])
        core_len = len(core)
        i = 1
        if not self.quiet:
            print("iteration 1 : ", core_len, "sum_abs_LLHR:" , sum([l if l > 0 else -l for l in likelies.values()]), file = sys.stderr)
        for mag in self:
            if len(core) > 0:
                mag.new_completness = 100*len(self.genome2gcs[mag.name].intersection(core))/len(core)
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
                    mag.new_completness = 100*len(self.genome2gcs[mag.name].intersection(core))/len(core)
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
            pp += "in a fcore of {core_len} traits with a total sum of loglikelihood-ratios {llhr:.2f} and a corrected \n"
            pp += "mean completness of {mean_new:.2f}, resulting to a estimated mean traits per genome count of {trait_count:.2f}\n"
            pp = pp.format(name = self.name, nb_mags = len(self), core_len = core_len, mean_start = mean([b.completeness for b in self]),
                        mean_new =  mean([b.new_completness for b in self]), llhr =  sum([l if l > 0 else -l for l in likelies.values()]),
                        trait_count = mean([100*len(self.genome2gcs[b.name])/b.new_completness for b in self]))
            print(pp, file = sys.stderr)
        self.iterations = i -1
        self.core = core
        self._core_computed = True
        return likelies

    def _core_prob(self, gene_cluster, complet = "checkm"):
        comp = lambda mag : (mag.completeness if complet =="checkm" else mag.new_completness)/100
        presence = [log10(comp(mag)) for mag in self if gene_cluster in self.genome2gcs[mag.name]]
        abscence = [log10(1 - comp(mag)) for mag in self if gene_cluster not in self.genome2gcs[mag.name]]
        return sum(presence + abscence)

    def get_genecluster_poolsize(self):
        if not self._genecluster_poolsize:
            self._genecluster_poolsize = sum([len(c.genomes) for c in  self.gene_clusters])
        return self._genecluster_poolsize

    def _pange_prob(self, gene_cluster, core_size, complet = "checkm"):
#        pool_size = sum(self.gene_clustersCounts.values())
        pool_size = self.get_genecluster_poolsize()
        comp = lambda mag : (mag.completeness if complet =="checkm" else mag.new_completness)/100
        current_genomes = gene_cluster.genomes

        #presence = [1 - (1-self.gene_clustersCounts[gene_clusters]/pool_size)**(len(mag.gene_clusterss)-(core_size*comp(mag))) for mag in self if gene_clusters in mag.gene_clusterss]
        #abscence = [ (1-self.gene_clustersCounts[gene_clusters]/pool_size)**(len(mag.gene_clusterss)-(core_size*comp(mag))) for mag in self if gene_clusters not in mag.gene_clusterss]

#        mag_prob = {mag : ( 1-1/pool_size )**len(mag.gene_clusterss) for mag in self}
#        mag_prob = {mag : ( 1-1/pool_size )**(len(mag.gene_clusterss)-(core_size*comp(mag))) for mag in self}
#        mag_prob = {mag : ( 1-self.gene_clustersCounts[gene_clusters]/pool_size )**(len(mag.gene_clusterss)-(core_size*comp(mag))) for mag in self}

#        mag_prob = {mag : ( 1-self.gene_clustersCounts[gene_clusters]/pool_size )**(len(mag.gene_clusterss)-(core_size*comp(mag))) for mag in self}
        mag_prob = {mag : ( 1-len(gene_cluster.genomes)/pool_size)**len(self.genome2gcs[mag.name]) for mag in self}

        presence = [ log10(1 -   mag_prob[mag]) if mag_prob[mag] < 1 else MIN_PROB                for mag in self if mag.name in current_genomes]
        abscence = [ log10(      mag_prob[mag]) if mag_prob[mag] > 0 else log10(1-(10**MIN_PROB)) for mag in self if mag.name not in current_genomes]

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
        if not stats['core']:
            stats['core'] = []
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
            out_dict[k]['genomes'] += [gg.name for gg in g.genomes]
            out_dict[k]['genome_occurences'] = len(g.genomes)

        for k,v in out_dict.items():
            v['mean_copy_per_genome'] = "NA" if not v['genes'] else len(v['genes'])/len(v['genomes'])
            v['genes'] = ";".join(v['genes'])
            v['genomes'] = ";".join(v['genomes'])

        header = ['trait_name','type', 'genome_occurences', 'log_likelihood_to_be_core', 'mean_copy_per_genome','genomes', 'genes']
        genome_line = "genomes=" + ";".join(["{}:prior_complete={}:posterior_complete={}".format(k.name, k.completeness, k.new_completness) for k in self])
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
            prior_complete=mean([b.completeness for b in self]),
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

    def write_checkm_file(self, path, for_temp_files = False):
        with open(path, "w") as handle:
            handle.writelines([ "Bin Id\tCompleteness\tContamination\n"])
            if for_temp_files:
                handle.writelines([f"{os.path.basename(t.nucleotide_file)[:-4]}\t{t.completeness}\t{t.redundancy}\n" for t in self])
            else :
                handle.writelines([f"{t.name}\t{t.completeness}\t{t.redundancy}\n" for t in self])


    def get_anis(self, method = "fastANI", block_size = 500, threads=20, recompute = False):
        if not hasattr(self, 'anis') or recompute:
            if not recompute:
                anis = self.db.get_anis(self)
            else :
                anis = dict()
            not_in = { (g1.name, g2.name)  for g1 in self for g2 in self if (g1.name, g2.name) not in anis}
            to_do = {gg for g in not_in for gg in g}
            if len(to_do) > 0 :
                if method == "fastANI":
                    if not shutil.which('fastANI'):
                        print("You need fastANI if you do not provide a file with pairwise similarities, either install it or provide pairwise similarities (see doc...)", file = sys.stderr)
                        sys.exit(-1)
                    fastani_file = tempfile.NamedTemporaryFile().name

                    mags2file = {b : self[b].nucleotide_file for b in to_do}
                    file2mags = {v.split("/")[-1][:-4] : k for k,v in mags2file.items()}
                    mags = list(mags2file.values())

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
                        out_dists = {(os.path.basename(l.split()[0]), os.path.basename(l.strip().split()[1])) : {'ani' : float(l.split()[2]), 'query_chunks' : float(l.split()[3]), 'reference_chunks' : float(l.split()[4]) } for l in handle}
                        out_dists = {( ".".join(k[0].split(".")[:-1]) if any([k[0].endswith(ext) for ext in FASTA_EXTS]) else k[0],
                                       ".".join(k[1].split(".")[:-1]) if any([k[1].endswith(ext) for ext in FASTA_EXTS]) else k[1] ): v
                                       for k,v in out_dists.items() }
                    os.remove(fastani_file)
                    new_anis = {(file2mags[k[0]], file2mags[k[1]]) : v for k,v in out_dists.items() if (file2mags[k[0]], file2mags[k[1]]) not in anis}
                elif method == "sourmash" :
                    if not self.quiet:
                        print("Computing ANIs using sourmash containment as a proxy:")
                    new_anis = { (g1.name,g2.name) : {'ani' : 100*(1 - g1.signature.containment_ani(g2.signature).dist), 'query_chunks' : -2, 'reference_chunks' : -2 } for g1 in tqdm(self) for g2 in self if (g1.name, g2.name) in not_in}
                else :
                    print("No other method for ani computation implemented yet")
                    sys.exit()
                for pair in not_in:
                        if pair not in new_anis:
                            new_anis[(a,b)] = {'ani' : -1, 'query_chunks' : -1, 'reference_chunks' : -1}
                tt = list(new_anis.items())
                if recompute:
                    if not mOTUlizer._quiet_:
                        print("updating anis")
                    for i in tqdm(list(range(0, len(tt), 1_000_000))):
                        self.db.update_anis(dict(tt[i:(i+1_000_000)]))
                else:
                    if not mOTUlizer._quiet_:
                        print("inserting anis")
                    for i in tqdm(list(range(0, len(tt), 1_000_000))):
                        self.db.add_anis(dict(tt[i:(i+1_000_000)]))
                anis.update(new_anis)
            self.anis = anis
        return self.anis

    def cluster_MetaBins(self, ani_cutoff = 95, prefix = "mOTU_", mag_complete = 40, mag_redundancy = 5, sub_complete = 0, sub_redundancy = 100, threads = 1, method = "fastANI"):
        import igraph

        dist_dict = self.get_anis(threads = threads, method = method)

        if not self.quiet:
            print("seeding bin-graph", file = sys.stderr )

        all_bins = {a.name : a for a in self}

        good_mag = lambda b : self[b].completeness > mag_complete and self[b].redundancy < mag_redundancy
        decent_sub = lambda b : self[b].completeness > sub_complete and self[b].redundancy < sub_redundancy and not good_mag(b)
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
                dd_dicts[c1][(k[0], k[1])] = v

        motus = [ mOTU(genomes = [self[g] for g in gs], name = prefix + str(i).zfill(zeros), quiet = True) for i, gs in enumerate(genome_clusters)]
        for motu, dists in zip(motus, dd_dicts):
            motu.load_anis(dists)

        return motus

    def __iter__(self):
       return mOTUIterator(self)

    def items(self):
        return mOTUItemIterator(self)

    def parallel_annotate(self, method = "prokka", threads = 24, temp_dir = None, tool_args = dict()):
        to_annots = [(k.name,temp_dir, method, tool_args) for k in self if not k.has_features()]
        print(f"Annotating {len(to_annots)} genomes of {self.name} with {method}")
        pool = mp.Pool(threads)
        to_commit =  pool.starmap_async(_annot_single, to_annots).get()
        if not mOTUlizer._quiet_:
            print("pushing features to db ")
        for name, features, source  in tqdm(to_commit):
            self.db.add_features(name, features = features, source = source)
        self.db.commit()

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
