from mOTUlizer.classes.MetaBin import MetaBin
from numpy import mean, prod, std, ceil
import subprocess
import tempfile
import os
from mOTUlizer.config import *
from numpy import inf
from random import shuffle, choice
import pandas
from tqdm import tqdm

class mOTU:
    def __len__(self):
        return len(self.members)

    def __repr__(self):
        return "< {tax} mOTU {name}, of {len} members >".format(name = self.name, len = len(self), tax = self.consensus_tax()[0].split(";")[-1])

    def __init__(self, name, members, data_pack, precomp_ani = None, funct_derep = 0.95, min_complet = 70):
        self.name = name
        self.members = [MetaBin(m, **data_pack) for m in members]
        self.core = None
        self.cogCounts = {c : 0 for c in set.union(*[mag.cogs for mag in self.members])}
        print("Counting cogs")
        for mag in tqdm(self.members):
            for cog in mag.cogs:
                    self.cogCounts[cog] += 1
        non_singles = {c for c,v in self.cogCounts.items() if v > 1}
        self.likelies = self.__core_likelyhood()
        if precomp_ani:
            memb_dict = {m.name : m for m in self}
            self.fastani_dict = { (a,b) : precomp_ani[(a,b)] for a in memb_dict for b in memb_dict if a != b  and(a,b) in precomp_ani }


    def __getitem__(self, i):
        return self.members[i]

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

    def mean_fastani(self) :
        return mean(list(self.fastani_matrix().values()))

    def get_otu_stats(self):
        tax = self.consensus_tax()
        tax_cert = ";".join([str(t) for t in tax[1]])

        stats = { "otu" : self.name,
        "nb_genomes" : len(self),
        "taxonomy" : tax[0],
        "certainty" : tax_cert,
        "mean_fastani" : str(self.mean_fastani()),
        "mean_COG_overlap" : str(self.mean_overlap()),
        "core_size" : len(self.core),
        "full_pan" : self.get_pangenome_size(singletons = True),
        "nosingle_pan" : self.get_pangenome_size(singletons = False),
        "mean_cog_count" : self.nb_cogs(),
        }
        return stats

    def get_count_stats(self):
        tax = self.consensus_tax()
        tax_cert = ";".join([str(t) for t in tax[1]])
        data = self.count_distribution()

        stats = { "otu" : self.name,
        "nb_genomes" : len(self),
        "taxonomy" : tax[0],
        "certainty" : tax_cert,
        "mean_fastani" : str(self.mean_fastani()),
        "mean_COG_overlap" : str(self.mean_overlap()),
        "core_size" : len(self.core),
        "full_pan" : self.get_pangenome_size(singletons = True),
        "nosingle_pan" : self.get_pangenome_size(singletons = False),
        "mean_cog_count" : self.nb_cogs(),
        }

        out_dat = []
        for k,v in data.items():
            tt = stats.copy()
            tt.update({'fract' : k/len(self), 'cog_count' : v})
            out_dat += [tt]
        return out_dat

    def __core_likelyhood(self, max_it = 10):
        print("Computing lieklihoods:")
        likelies = {cog : self.__core_likely(cog) for cog in tqdm(self.cogCounts)}
        self.core = set([c for c, v in likelies.items() if v > 1])
        core_len = len(self.core)
        i = 1
        if VERBOSE:
            print("iteration 1 : ", core_len )
        for mag in self:
            mag.new_completness = 100*len(mag.cogs.intersection(self.core))/len(self.core)
            mag.new_completness = mag.new_completness if mag.new_completness < 95 else 95
        for i in range(2,max_it):
            likelies = {cog : self.__core_likely(cog, complet = "new", core_size = core_len) for cog in self.cogCounts}
            self.core = set([c for c, v in likelies.items() if v > 1])
            new_core_len = len(self.core)
            for mag in self:
                mag.new_completness = 100*len(mag.cogs.intersection(self.core))/len(self.core)
            if VERBOSE:
                print("iteration",i, ": ", new_core_len)
            if new_core_len == core_len:
               break
            else :
                core_len =new_core_len
        if VERBOSE:
            print("Mean checkM completeness", mean([b.checkm_complet for b in self]), ", mean new_completness",  mean([b.new_completness for b in self]))
        self.iterations = i -1
        return likelies

    def __core_prob(self, cog, complet = "checkm"):
        comp = lambda mag : (mag.checkm_complet if complet =="checkm" else mag.new_completness)/100
        presence = [comp(mag) for mag in self if cog in mag.cogs]
        abscence = [1 - comp(mag) for mag in self if cog not in mag.cogs]
        return prod(presence + abscence)

    def __pange_prob(self, cog, core_size, complet = "checkm"):
        pool_size = sum(self.cogCounts.values())
        comp = lambda mag : (mag.checkm_complet if complet =="checkm" else mag.new_completness)/100
        presence = [1 - (1-self.cogCounts[cog]/pool_size)**(len(mag.cogs)-(core_size*comp(mag))) for mag in self if cog in mag.cogs]
        abscence = [ (1-self.cogCounts[cog]/pool_size)**(len(mag.cogs)-(core_size*comp(mag))) for mag in self if cog not in mag.cogs]
        return prod(presence + abscence)

    def __core_likely(self, cog, complet = "checkm", core_size = 0):
        pange_prob = self.__pange_prob(cog, core_size, complet)
        if pange_prob != 0:
            return self.__core_prob(cog, complet)/pange_prob
        else :
            return inf

    def nb_cogs(self):
        return mean([b.estimate_nb_cogs() for b in self if b.new_completness > 40 ])

    def get_pangenome_size(self, singletons = False):
        return len([k for k,v in self.cogCounts.items() if k not in self.core and v > (0 if singletons else 1)])

    def count_distribution(self):
        counts = { i : 0 for i in range(1,len(self)+1) }
        for k,v in self.cogCounts.items():
            if k in self.core:
                counts[len(self)] += 1
            else :
                counts[v] +=1
        return counts


    def simulate_rarefaction(self, reps = 100, genomes = 1000):
        cog_pool = [k for k, v in self.cogCounts.items() for i in range(v) if k not in self.core]
        nosings_cog_pool = [k for k, v in self.cogCounts.items() for i in range(v) if k not in self.core and v > 1]
        draws = self.nb_cogs()-len(self.core)
        series = []
        for rep in tqdm(range(reps)):
            random_auxs = [{choice(cog_pool) for i in range(int(ceil(draws)))} for m in range(genomes)]
            tt = self.rarefy_pangenome(custom_cogs = random_auxs)
            for k, v in tt.items():
                v['rep_sim'] = rep
                v['singletons'] = False
                v['sim'] = True
            series += list(tt.values())

            random_auxs = [{choice(nosings_cog_pool) for i in range(int(ceil(draws)))} for m in range(genomes)]
            tt = self.rarefy_pangenome(custom_cogs = random_auxs)
            for k, v in tt.items():
                v['rep_sim'] = rep
                v['singletons'] = True
                v['sim'] = True
            series += list(tt.values())

        tt = self.rarefy_pangenome()
        for k, v in tt.items():
                v['rep_sim'] = -1
                v['singletons'] = False
                v['sim'] = False
        series += list(tt.values())

        tt = self.rarefy_pangenome(singletons = True)
        for k, v in tt.items():
                v['rep_sim'] = -1
                v['singletons'] = True
                v['sim'] = False
        series += list(tt.values())
        return series


    def rarefy_pangenome(self, reps = 100, singletons = False, custom_cogs = None):
        def __min_95(ll):
            ll.sort_values()
            return list(ll.sort_values())[round(len(ll)*0.05)]

        def __max_95(ll):
            ll.sort_values()
            return list(ll.sort_values())[round(len(ll)*0.95)]

        def __genome_count(ll):
            return ll[0]

        __genome_count.__name__ = "genome_count"
        __max_95.__name__ = "max_95"
        __min_95.__name__ = "min_95"

        pange = set.union(*custom_cogs) if custom_cogs else {k for k,v in self.cogCounts.items() if k not in self.core and v > (0 if singletons else 1)}
        series = []
        for i in range(reps):
            series += [{ 'rep' : i , 'genome_count' : 0, 'pangenome_size' : 0}]
            m = custom_cogs if custom_cogs else [m.cogs for m in self.members.copy()]
            shuffle(m)
            founds = set()
            for j,mm in enumerate(m):
                founds = founds.union(mm.intersection(pange))
                series += [{ 'rep' : i , 'genome_count' : j+1, 'pangenome_size' : len(founds)}]

        t_pandas = pandas.DataFrame.from_records(series)[['genome_count', 'pangenome_size']]
        t_pandas = t_pandas.groupby('genome_count').agg({'genome_count' : [__genome_count] ,'pangenome_size' : [mean, std, __min_95, __max_95]} )
        t_pandas.columns = ["rr_" + p[1] for p in t_pandas.columns]
        t_pandas = t_pandas.to_dict(orient="index")
        tt = self.get_otu_stats()
        for v in t_pandas.values():
            v.update(tt)
        return t_pandas
