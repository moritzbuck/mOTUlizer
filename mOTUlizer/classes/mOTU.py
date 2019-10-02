from mOTUlizer.classes.MetaBin import MetaBin
from numpy import mean, prod
import subprocess
import tempfile
import os
from mOTUlizer.config import *
from numpy import inf

class mOTU:
    def __len__(self):
        return len(self.members)

    def __repr__(self):
        return "< {tax} mOTU {name}, of {len} members >".format(name = self.name, len = len(self), tax = self.consensus_tax()[0].split(";")[-1])

    def __init__(self, name, members, data_pack):
        self.name = name
        self.members = [MetaBin(m, **data_pack) for m in members]
        self.core = None

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
        "mean_COG_overlap" : str(self.mean_overlap())
        }
        return stats

    def core_likelyhood(self):
        cogCounts = {c : 0 for c in set.union(*[mag.cogs for mag in self])}
        for mag in self:
            for cog in mag.cogs:
                    cogCounts[cog] += 1
        likelies = {cog : self.__core_likely(cog, cogCounts) for cog in cogCounts}
        self.core = set([c for c, v in likelies.items() if v > 1])
        core_len = len(self.core)
        print("iteration 1 : ", core_len )
        for mag in self:
            mag.new_completness = 100*len(mag.cogs.intersection(self.core))/len(self.core)
        for i in range(10):
            likelies = {cog : self.__core_likely(cog, cogCounts, complet = "new", core_size = core_len) for cog in cogCounts}
            self.core = set([c for c, v in likelies.items() if v > 1])
            new_core_len = len(self.core)
            for mag in self:
                mag.new_completness = 100*len(mag.cogs.intersection(self.core))/len(self.core)
            print("iteration",i, ": ", new_core_len)
            if new_core_len == core_len:
                break
            else :
                core_len =new_core_len
        return likelies

    def __core_prob(self, cog, complet = "checkm"):
        comp = lambda mag : (mag.checkm_complet if complet =="checkm" else mag.new_completness)/100
        presence = [comp(mag) for mag in self if cog in mag.cogs]
        abscence = [1 - comp(mag) for mag in self if cog not in mag.cogs]
        return prod(presence + abscence)

    def __pange_prob(self, cog, pool, core_size, complet = "checkm"):
        pool_size = sum(pool.values())
        comp = lambda mag : (mag.checkm_complet if complet =="checkm" else mag.new_completness)/100
        presence = [1 - (1-pool[cog]/pool_size)**(len(mag.cogs)-(core_size*comp(mag))) for mag in self if cog in mag.cogs]
        abscence = [ (1-pool[cog]/pool_size)**(len(mag.cogs)-(core_size*comp(mag))) for mag in self if cog not in mag.cogs]
        return prod(presence + abscence)

    def __core_likely(self, cog, pool, complet = "checkm", core_size = 0):
        pange_prob = self.__pange_prob(cog, pool, core_size, complet)
        if pange_prob != 0:
            return self.__core_prob(cog, complet)/pange_prob
        else :
            return inf
