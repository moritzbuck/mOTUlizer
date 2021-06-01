from mOTUlizer.classes.MetaBin import MetaBin
from mOTUlizer.classes.COGs import *
import subprocess
import tempfile
import os
from mOTUlizer.config import *
from random import shuffle, choice, choices
from math import log10
import sys
import json


mean = lambda x : sum(x)/len(x)

class mOTU:
    def __len__(self):
        return len(self.members)

    def __repr__(self):
        return "< {tax} mOTU {name}, of {len} members >".format(name = self.name, len = len(self), tax =  None ) #self.consensus_tax()[0].split(";")[-1])

    def __init__(self, **kwargs):
        self.quiet = kwargs['quiet'] if "quiet" in kwargs else False
        self.likelies = None
        self.mock = []
        if "cog_dict" in kwargs:
            self.__for_mOTUpan(**kwargs)

        if "bins" in kwargs:
            self.__from_bins(**kwargs)


    def __for_mOTUpan(self, name, faas, cog_dict, checkm_dict, threads = 4, precluster = False, max_it = 20, method = None, quiet=False):
        self.name = name
        self.faas = faas
        if not quiet:
            print("Creating mOTU for mOTUpan", file = sys.stderr)
        if  not cog_dict :
            tt = compute_COGs(self.faas, name = name, precluster = precluster, threads = threads)
            self.cog_dict = tt['genome2cogs']
            self.aa2cog = tt['aa2cog']
        else :
            self.cog_dict = cog_dict
            self.aa2cog = {}

        if checkm_dict == "length_seed" :
            max_len = max([len(cogs) for cogs in self.cog_dict.values()])
            checkm_dict = {}
            for f in self.cog_dict:
                checkm_dict[f] = 100*len(self.cog_dict[f])/max_len

        self.members = [MetaBin(bin_name, cogs = self.cog_dict[bin_name], faas = self.faas.get(bin_name), fnas = None, complet = checkm_dict.get(bin_name)) for bin_name in self.cog_dict.keys()]
        self.cogCounts = {c : 0 for c in set.union(set([cog for mag in self.members for cog in mag.cogs]))}
        for mag in self.members:
            for cog in mag.cogs:
                    self.cogCounts[cog] += 1

        self.core = {cog for cog, counts in self.cogCounts.items() if (100*counts/len(self)) > mean(checkm_dict.values())}

        self.method = method
        if self.method :
            self.likelies = self.__core_likelyhood(max_it = max_it)


    def __getitem__(self, i):
        return self.members[i]


    def roc_values(self, boots):
        if boots > 0:
            from mOTUlizer.classes.MockData import MockmOTU
            mean = lambda data: float(sum(data)/len(data))
            variance 	= lambda data, avg: sum([x**2 for x in [i-avg for i in data]])/float(len(data))
            std_dev = lambda data: variance(data, mean(data))**0.5

            while len(self.mock) < boots:
                print("Running bootstrap {}/{}".format(len(self.mock)+1, boots), file = sys.stderr)
                completnesses = {"Genome_{}".format(i) : c.new_completness for i,c in enumerate(self)}

                accessory = sorted([v for k,v in self.cogCounts.items() if k not in self.core])
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
                     'sd_lowest_false' : std_dev([m.lowest_false for m in self.mock])}
        else:
            return { 'mean_recall' : "NA",
                     'sd_recall' : "NA",
                     'mean_fpr' : "NA",
                     'sd_fpr' : "NA",
                     'mean_lowest_false' : "NA",
                     'sd_lowest_false' : "NA"}


    def avg_cog_content(self):
        return sum([len(m.cogs) for m in self.members])/len(self)

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
            self.overlap_dict = {(i,j) : len(i.overlap(j))/len(i.cogs) for i in self.members for j in self.members if i != j}
        return self.overlap_dict

    def mean_overlap(self) :
        return mean(list(self.overlap_matrix().values()))

    def fastani_matrix(self):
        if not hasattr(self, 'fastani_dict'):
            if not shutil.which('fastANI'):
                print("You need fastANI if you do not provide a file with pairwise similarities, either install it or provide pairwise similarities (see doc...)", file = sys.stderr)
                sys.exit(-1)
            cmd = "fastANI --ql {temp_file}  --rl {temp_file}  -o {temp_out} -t {threads} 2> /dev/null"
            temp_inp = tempfile.NamedTemporaryFile()
            temp_out = tempfile.NamedTemporaryFile()

            with open(temp_inp.name, "w") as handle:
                handle.writelines([t.genome + "\n" for t in self.members ])

            cmd = cmd.format(temp_file = temp_inp.name, temp_out = temp_out.name, threads = THREADS)

            subprocess.call(cmd, shell=True)

            with open(temp_out.name ) as handle:
                self.fastani_dict = {}
                for l in handle:
                    q = os.path.basename(l.split()[0][:-4])
                    s = os.path.basename(l.split()[1][:-4])
                    if q  != s:
                        value = float(l.split()[2])
                        self.fastani_dict[(q,s)] = value
        return self.fastani_dict


    def get_stats(self):
        out = {}

        out[self.name] = { "nb_genomes" : len(self),
        "core" : list(self.core) if self.core else None,
        "aux_genome" : [k for k,v in self.cogCounts.items() if k not in self.core] if self.core else None ,
        "singleton_cogs" : [k for k,v in self.cogCounts.items() if k not in self.core if v == 1] if self.core else None,
        "cogs" : None if self.cog_dict is None else {'genome' : {k : list(v) for k,v in self.cog_dict.items()}, 'aa' : self.aa2cog} if self.aa2cog else ({k : list(v) for k,v  in self.cog_dict.items()} if self.cog_dict else None),
        "mean_ANI" : self.get_mean_ani() if (hasattr(self, 'fastani_dict') or all([hasattr(g, "genome") for g in self])) else None,
        "ANIs" : [[k[0], k[1], v] for k, v in self.fastani_matrix().items()] if (hasattr(self, 'fastani_dict')  or all([hasattr(g, "genome") for g in self])) else None,
        "genomes" : [v.get_data() for v in self],
        "likelies" : self.likelies
        }
        return out

    def get_representative(self, method = "complex", max_redund = 5, min_complete = 95):
            tt = [v.get_data() for v in self]
            data = { t['name'] : t for t in tt}
            if all([v['checkm_contamin'] > max_redund for v in data.values()]):
                return max( [ (k , v['checkm_complet']) for k,v in data.items()], key = lambda x: x[1])[0]
            data = {k : v for k,v in data.items() if v['checkm_contamin'] < max_redund}
            if any([v['checkm_complet'] > min_complete for v in data.values()])  :
                data = {k : v for k,v in data.items() if v['checkm_complet'] > min_complete}
                best_redund = min(data.items(), key = lambda x : x[1]['checkm_contamin'])[1]['checkm_contamin']
                return max( [ (k , v['checkm_complet']) for k,v in data.items() if v['checkm_contamin'] == best_redund], key = lambda x: x[1])[0]
            else:
                return max( [ (k , v['checkm_complet']) for k,v in data.items()], key = lambda x: x[1])[0]

    # def get_representative(tt, max_redund = 5, min_complete = 95):
    #         data = { t['name'] : t for t in tt}
    #
    #         data = {k : v for k,v in data.items() if v['checkm_contamin'] < max_redund}
    #         if any([v['checkm_complet'] > min_complete for v in data.values()])  :
    #             data = {k : v for k,v in data.items() if v['checkm_complet'] > min_complete}
    #             best_redund = min(data.items(), key = lambda x : x[1]['checkm_contamin'])[1]['checkm_contamin']
    #             return min( [ (k , v['checkm_contamin']) for k,v in data.items() if v['checkm_contamin'] == best_redund], key = lambda x: x[1])[0]
    #         elif len(data) >0 :
    #             return max( [ (k , v['checkm_complet']) for k,v in data.items()], key = lambda x: x[1])[0]
    #         else :
    #             return None

    def get_mean_ani(self):
        dist_dict = self.fastani_matrix()
        dists = [dist_dict.get((a.name,b.name)) for a in self for b in self if a != b ]
        missing_edges = sum([d is None for d in dists])
        found_edges = [d is None for d in dists]

        return {'mean_ANI' : sum([d for d in dists if d])/len(found_edges) if len(found_edges) > 0 else None, 'missing_edges' : missing_edges, 'total_edges' : len(found_edges) + missing_edges}

    def __core_likelyhood(self, max_it = 20, likeli_cutof = 0 ):
        likelies = {cog : self.__core_likely(cog) for cog in self.cogCounts}
        self.core = set([c for c, v in likelies.items() if v > likeli_cutof])
        core_len = len(self.core)
        i = 1
        if not self.quiet:
            print("iteration 1 : ", core_len, "sum_abs_LLHR:" , sum([l if l > 0 else -l for l in likelies.values()]), file = sys.stderr)
        for mag in self:
            if len(self.core) > 0:
                mag.new_completness = 100*len(mag.cogs.intersection(self.core))/len(self.core)
            else :
                mag.new_completness = 0
            mag.new_completness = mag.new_completness if mag.new_completness < 99.9 else 99.9
            mag.new_completness = mag.new_completness if mag.new_completness > 0 else 0.01
        for i in range(2,max_it):
            likelies = {cog : self.__core_likely(cog, complet = "new", core_size = core_len) for cog in self.cogCounts}
            old_core = self.core
            self.core = set([c for c, v in likelies.items() if v > likeli_cutof])
            new_core_len = len(self.core)
            for mag in self:
                if len(self.core) > 0:
                    mag.new_completness = 100*len(mag.cogs.intersection(self.core))/len(self.core)
                else :
                    mag.new_completness = 0
                mag.new_completness = mag.new_completness if mag.new_completness < 99.9 else 99.9
                mag.new_completness = mag.new_completness if mag.new_completness > 0 else 0.01
            if not self.quiet:
                print("iteration",i, ": ", new_core_len, "sum_abs_LLHR:" , sum([l if l > 0 else -l for l in likelies.values()]), file = sys.stderr)
            if self.core == old_core:
               break
            else :
                core_len = new_core_len

        if not self.quiet:
            json.dump({ self.name : {"nb_mags" : len(self), "core_len" : core_len, "mean_starting_completeness" :  mean([b.checkm_complet for b in self]), "mean_new_completness" : mean([b.new_completness for b in self]), "LHR"  :  sum([l if l > 0 else -l for l in likelies.values()]), "mean_est_binsize" : mean([100*len(b.cogs)/b.new_completness for b in self])}}, sys.stderr )
            print(file = sys.stderr)
        self.iterations = i -1
        return likelies

    def __core_prob(self, cog, complet = "checkm"):
        comp = lambda mag : (mag.checkm_complet if complet =="checkm" else mag.new_completness)/100
        presence = [log10(comp(mag)) for mag in self if cog in mag.cogs]
        abscence = [log10(1 - comp(mag)) for mag in self if cog not in mag.cogs]
        return sum(presence + abscence)

    def __pange_prob(self, cog, core_size, complet = "checkm"):
#        pool_size = sum(self.cogCounts.values())
        pool_size = sum([c for k,c in  self.cogCounts.items()])
        comp = lambda mag : (mag.checkm_complet if complet =="checkm" else mag.new_completness)/100
        #presence = [1 - (1-self.cogCounts[cog]/pool_size)**(len(mag.cogs)-(core_size*comp(mag))) for mag in self if cog in mag.cogs]
        #abscence = [ (1-self.cogCounts[cog]/pool_size)**(len(mag.cogs)-(core_size*comp(mag))) for mag in self if cog not in mag.cogs]

#        mag_prob = {mag : ( 1-1/pool_size )**len(mag.cogs) for mag in self}
#        mag_prob = {mag : ( 1-1/pool_size )**(len(mag.cogs)-(core_size*comp(mag))) for mag in self}
#        mag_prob = {mag : ( 1-self.cogCounts[cog]/pool_size )**(len(mag.cogs)-(core_size*comp(mag))) for mag in self}

#        mag_prob = {mag : ( 1-self.cogCounts[cog]/pool_size )**(len(mag.cogs)-(core_size*comp(mag))) for mag in self}
        mag_prob = {mag : ( 1-self.cogCounts[cog]/pool_size)**len(mag.cogs) for mag in self}

        presence = [ log10(1 -   mag_prob[mag]) if mag_prob[mag] < 1 else MIN_PROB                for mag in self if cog in mag.cogs]
        abscence = [ log10(      mag_prob[mag]) if mag_prob[mag] > 0 else log10(1-(10**MIN_PROB)) for mag in self if cog not in mag.cogs]

        #abscence = [ 1-self.cogCounts[cog]/len(self)*comp(mag) for mag in self if cog not in mag.cogs]
        #presence = [ self.cogCounts[cog]/len(self)*comp(mag) for mag in self if cog not in mag.cogs]

        return sum(presence + abscence)

    def __core_likely(self, cog, complet = "checkm", core_size = 0):
        pange_prob = self.__pange_prob(cog, core_size, complet)
        return self.__core_prob(cog, complet) - pange_prob

    def nb_cogs(self):
        return mean([b.estimate_nb_cogs() for b in self if b.new_completness > 40 ])

    def get_pangenome_size(self, singletons = False):
        return len([k for k,v in self.cogCounts.items() if k not in self.core and v > (0 if singletons else 1)])


    def __from_bins(self, bins, name,dist_dict = None ):
        self.name = name

        self.members = bins
        self.core = None
        self.fastani_dict = dist_dict
        self.aa2cog = None
        self.likelies = None
        self.cog_dict = None

    @classmethod
    def cluster_MetaBins(cls , all_bins, dist_dict, ani_cutoff = 95, prefix = "mOTU_", mag_complete = 40, mag_contamin = 5, sub_complete = 0, sub_contamin = 100):
        import igraph

        print("seeding bin-graph", file = sys.stderr )

        all_bins = {a.name : a for a in all_bins}

        tt = [(k, v.checkm_complet, v.checkm_contamin) for k, v in all_bins.items() if v.checkm_complet > 0]

        good_mag = lambda b : all_bins[b].checkm_complet > mag_complete and all_bins[b].checkm_contamin < mag_contamin
        decent_sub = lambda b : all_bins[b].checkm_complet > sub_complete and all_bins[b].checkm_contamin < sub_contamin and not good_mag(b)
        good_pairs = [k for k,v  in dist_dict.items() if v > ani_cutoff and dist_dict.get((k[1],k[0]), 0) > ani_cutoff and good_mag(k[0]) and good_mag(k[1])]
        species_graph = igraph.Graph()
        vertexDeict = { v : i for i,v in enumerate(set([x for k in good_pairs for x in k]))}
        rev_vertexDeict = { v : i for i,v in vertexDeict.items()}
        species_graph.add_vertices(len(vertexDeict))
        species_graph.add_edges([(vertexDeict[k[0]], vertexDeict[k[1]]) for k in good_pairs])

        print("getting clusters", file = sys.stderr)

        genome_clusters = [[rev_vertexDeict[cc] for cc in c ] for c in species_graph.components(mode=igraph.STRONG)]

        mean = lambda l : sum([len(ll) for ll in l])/len(l)

        print("recruiting to graph of the", len(genome_clusters) ," mOTUs of mean length", mean(genome_clusters), file = sys.stderr)


        left_pairs = {k : v for k, v in dist_dict.items() if v > ani_cutoff and k[0] != k[1] and ((decent_sub(k[0]) and good_mag(k[1])) or (decent_sub(k[1]) and good_mag(k[0])))}
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
        motus = [ mOTU(bins = [all_bins[gg] for gg in gs], name = prefix + str(i).zfill(zeros), dist_dict = {(k,l) : dist_dict[(k,l)] for k in gs for l in gs if (k,l) in dist_dict}) for i, gs in enumerate(genome_clusters)]


        return motus
