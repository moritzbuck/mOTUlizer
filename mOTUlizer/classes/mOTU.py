from mOTUlizer.classes.MetaBin import MetaBin
from mOTUlizer.classes.COGs import *
import subprocess
import tempfile
import os
from mOTUlizer.config import *
from random import shuffle, choice
from math import log10
import sys

mean = lambda x : sum(x)/len(x)
class mOTU:
    def __len__(self):
        return len(self.members)

    def __repr__(self):
        return "< {tax} mOTU {name}, of {len} members >".format(name = self.name, len = len(self), tax = self.consensus_tax()[0].split(";")[-1])

    def __init__(self, name, faas, cog_dict, checkm_dict = {}):
        self.name = name
        self.faas = faas
        if  not cog_dict :
            tt = compute_COGs(self.faas, name = name + "COG")
            self.cog_dict = tt['genome2cogs']
            self.aa2cog = tt['aa2cog']
        else :
            self.cog_dict = cog_dict
            self.aa2cog = {}
        self.members = [MetaBin(bin_name, self.cog_dict[bin_name], self.faas.get(bin_name), checkm_dict.get(bin_name)) for bin_name in self.faas.keys()]
        self.core = None
        self.cogCounts = {c : 0 for c in set.union(*[mag.cogs for mag in self.members])}
        for mag in self.members:
            for cog in mag.cogs:
                    self.cogCounts[cog] += 1

        self.likelies = self.__core_likelyhood()

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


    def get_stats(self):
        out = {}
        out[self.name] = { "nb_genomes" : len(self),
        "core" : list(self.core),
        "aux_genome" : [k for k,v in self.cogCounts.items() if k not in self.core],
        "singleton_cogs" : [k for k,v in self.cogCounts.items() if k not in self.core if v == 1],
        "cogs" : {'genome' : {k : list(v) for k,v in self.cog_dict.items()}, 'aa' : self.aa2cog}
        }
        return out

    def __core_likelyhood(self, max_it = 20):
        likelies = {cog : self.__core_likely(cog) for cog in self.cogCounts}
        self.core = set([c for c, v in likelies.items() if v > 0])
        core_len = len(self.core)
        i = 1
        print("iteration 1 : ", core_len, "LHR:" , sum(likelies.values()), file = sys.stderr)
        for mag in self:
            mag.new_completness = 100*len(mag.cogs.intersection(self.core))/len(self.core)
            mag.new_completness = mag.new_completness if mag.new_completness < 99.9 else 99.9
            mag.new_completness = mag.new_completness if mag.new_completness > 0 else 0.01
        for i in range(2,max_it):
            likelies = {cog : self.__core_likely(cog, complet = "new", core_size = core_len) for cog in self.cogCounts}
            self.core = set([c for c, v in likelies.items() if v > 0])
            new_core_len = len(self.core)
            for mag in self:
                mag.new_completness = 100*len(mag.cogs.intersection(self.core))/len(self.core)
                mag.new_completness = mag.new_completness if mag.new_completness < 99.9 else 99.9
                mag.new_completness = mag.new_completness if mag.new_completness > 0 else 0.01

            print("iteration",i, ": ", new_core_len, "LHR:" , sum(likelies.values()), file = sys.stderr)
            if new_core_len == core_len:
               break
            else :
                core_len =new_core_len

        print("starting completeness", mean([b.checkm_complet for b in self]), ", mean new_completness",  mean([b.new_completness for b in self]), file = sys.stderr)
        self.iterations = i -1
        return likelies

    def __core_prob(self, cog, complet = "checkm"):
        comp = lambda mag : (mag.checkm_complet if complet =="checkm" else mag.new_completness)/100
        presence = [log10(comp(mag)) for mag in self if cog in mag.cogs]
        abscence = [log10(1 - comp(mag)) for mag in self if cog not in mag.cogs]
        return sum(presence + abscence)

    def __pange_prob(self, cog, core_size, complet = "checkm"):
        pool_size = sum(self.cogCounts.values())
        comp = lambda mag : (mag.checkm_complet if complet =="checkm" else mag.new_completness)/100
        #presence = [1 - (1-self.cogCounts[cog]/pool_size)**(len(mag.cogs)-(core_size*comp(mag))) for mag in self if cog in mag.cogs]
        #abscence = [ (1-self.cogCounts[cog]/pool_size)**(len(mag.cogs)-(core_size*comp(mag))) for mag in self if cog not in mag.cogs]

#        presence = [ log10(1 -   ( 1 - 1/len(self.cogCounts))**(len(mag.cogs)-(core_size*comp(mag)))) for mag in self if cog in mag.cogs]
#        abscence = [       log10(( 1 - 1/len(self.cogCounts))**(len(mag.cogs)-(core_size*comp(mag)))) for mag in self if cog not in mag.cogs]

        presence = [ log10(1 -   ( 1-self.cogCounts[cog]/pool_size )**(len(mag.cogs)-(core_size*comp(mag)))) for mag in self if cog in mag.cogs]
        abscence = [       log10(( 1-self.cogCounts[cog]/pool_size )**(len(mag.cogs)-(core_size*comp(mag)))) for mag in self if cog not in mag.cogs]


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
