from random import random, sample, choice
from mOTUlizer.classes.mOTU import mOTU, mean
from collections import defaultdict
from numpy.random import normal
from numpy import floor, mean

genome2guass = {}

class MockmOTU(mOTU):
    def __repr__(self) :
        return "< MockmOTU with {n} genomes, of average {c}% completness, with core/genome_len of {r} >".format(c = 100*self.mean_completeness, n = len(self), r = self.ratio)

    def __init__(self, name, core_len, nb_genomes, completeness, max_it = 20, accessory = None, method = None):

        core = {"CoreTrait_{}".format(i) for i in range(core_len)}

        if accessory is None:
            sub_dist = [int(nb_genomes/i) for i in range(2,1000) if int(nb_genomes/i) > 0] + [1]*100
            sub_dist = list(range(nb_genomes-1, 1,-1))
        else :
            sub_dist = accessory

        self.size_accessory = sum(sub_dist)
        self.mean_size_accessory = sum(sub_dist)/nb_genomes

        mock_genomes = dict()
        for k in range(nb_genomes):
            mock_genomes["Genome_{}".format(k)] = list(core)

        for i,v in enumerate(sub_dist):
            genomes = sample(list(mock_genomes.keys()), v if v < len(mock_genomes) else len(mock_genomes) )
            for g in genomes:
                mock_genomes[g] += ["AccessoryTrait_{}".format(i)]

        self.incompletes = {g : {vv for vv in v if random() < (completeness(g)/100)} for g, v in mock_genomes.items()}

        to_rm = []
        for k, v in self.incompletes.items():
            if len(v) == 0:
                if len(core) > 0:
                    self.incompletes[k] = choice(list(core))

        if core_len == 0:
            self.mean_completeness = "NA"
            self.completenesses = {c : 0 for c,v in self.incompletes.items()}
        else :
            self.mean_completeness = mean([len({vv for vv in v if vv.startswith("CoreTrait_")})/core_len for c,v in self.incompletes.items()])
            self.completenesses = {c : 100*len({vv for vv in v if vv.startswith("CoreTrait_")})/core_len for c,v in self.incompletes.items()}
    #        self.accessory = accessory
        self.mean_size = mean([len(m) for m in mock_genomes.values()])
        self.real_core_len = core_len

        zerifneg = lambda g: 0.001 if g < 0 else g
        super().__init__(name = name, faas = {}, cog_dict = self.incompletes, checkm_dict = { k : zerifneg(normal(v, 10)) for k,v in self.completenesses.items()}, max_it = max_it, method = method, quiet=True)
        if core_len == 0:
            self.recall = "NA"
            self.fpr = "NA"
        else :
            self.recall = len(core.intersection(self.core))/core_len
            self.fpr = sum([not c.startswith("CoreTrait_") for c in self.core])/len([not c.startswith("CoreTrait_") for c in self.core])

        self.lowest_false = {k : v for k,v in self.cogCounts.items() if k in self.core and k not in core}
        self.lowest_false = 1 if(len(self.lowest_false) ==0) else min(self.lowest_false.items(), key = lambda x : x[1])[1]/len(self)


    def mock_cog_stats(self):
        all_genes = set.union(*self.incompletes.values())
        outp = {t : {} for t in all_genes}
        for t,dd  in outp.items():
            dd['freq'] = sum([t in zz for zz in self.incompletes.values()])/len(self.incompletes)
            dd['core'] = t in self.core
            dd['type'] = "core" if t.startswith("CoreTrait_") else "accessory"
            dd['nb_genomes'] = len(self)
            dd['core_len'] = len(self.core)
            dd['real_core_len'] = self.read_core_len
            dd['llikelihood'] = self.likelies[t]
            dd['len_accessory_genome'] = len(all_genes) - dd['real_core_len']
        return outp

    @classmethod
    def guauss_completes(cls, g, mean_completeness = 60, stdev = 10):
        if g in genome2guass:
            return genome2guass[g]
        else :
            out_prob = 1000
            while(not (20 < out_prob < 99) ):
                out_prob = normal(mean_completeness, stdev)
            genome2guass[g] = out_prob
            return out_prob

    @classmethod
    def run_boots(cls):
        out = {}

        for i in range(1000, 2500, 500):
            for nb in range(10, 250, 15):
                for c in range(30, 100, 5):
                    MockData.genome2guass = {}
                    mockmotu = MockmOTU("complete_{}_core_size_{}_nbgenomes_{}".format(c,i,nb).format(c), i, nb, lambda g : MockmOTU.guauss_completes(g, mean_completeness = c, stdev = 10), max_it = 100)
                    out[mockmotu.name] = {}
                    out[mockmotu.name]['core_size'] = i
                    out[mockmotu.name]['nb_genomes'] = nb
                    out[mockmotu.name]['mean_completeness'] = mockmotu.mean_completeness
                    out[mockmotu.name]['mean_genome_size'] = mockmotu.mean_size
                    out[mockmotu.name]['recall'] = mockmotu.recall
                    out[mockmotu.name]['lowest_false'] = mockmotu.lowest_false
                    out[mockmotu.name]['accessory_genepool'] = mockmotu.size_accessory
                    out[mockmotu.name]['mean_new_completness'] = mean([b.new_completness for b in mockmotu])
        return out
