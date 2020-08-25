from random import random, sample
from mOTUlizer.classes.mOTU import mOTU, mean
from collections import defaultdict
from tqdm import tqdm
from numpy.random import normal
from numpy import floor

genome2guass = {}

with open("mOTUlizer/data/Efeacalis.csv") as handle:
    emp_disp = [int(l[:-1].split(",")[-1])/691 for l in handle if "accessory" in l]
gene_count = sum(emp_disp)

class MockmOTU(mOTU):
    def __repr__(self) :
        return "< MockmOTU with {n} genomes, of average {c}% completness, with core/genome_len of {r} >".format(c = 100*self.mean_completeness, n = len(self), r = self.ratio)

    def __init__(self, name, ratio, genome_len, nb_genomes, completeness):

        self.ratio = ratio
        core_len = int(genome_len*self.ratio)
        core = {"CoreTrait_{}".format(i) for i in range(core_len)}


        sub_dist = list(range(nb_genomes-1, 1,-1))

        mock_genomes = dict()
        for k in range(nb_genomes):
            mock_genomes["Genome_{}".format(k)] = list(core)

        for i,v in enumerate(sub_dist):
            genomes = sample(list(mock_genomes.keys()), v)
            for g in genomes:
                mock_genomes[g] += ["AccessoryTrait_{}".format(i)]

        self.incompletes = {g : {vv for vv in v if random() < (completeness(g)/100)} for g, v in mock_genomes.items()}
        self.mean_completeness = mean([len({vv for vv in v if vv.startswith("CoreTrait_")})/core_len for c,v in self.incompletes.items()])
        self.completenesses = {c : 100*len({vv for vv in v if vv.startswith("CoreTrait_")})/core_len for c,v in self.incompletes.items()}
#        self.accessory = accessory
        self.mean_size = mean([len(m) for m in mock_genomes.values()])
        self.read_core_len = core_len
        super().__init__(name = name, faas = {}, cog_dict = self.incompletes, checkm_dict = self.completenesses, max_it = 50)


    @classmethod
    def guauss_completes(cls, g, mean_completeness = 60, stdev = 10):
        if g in genome2guass:
            return genome2guass[g]
        else :
            out_prob = 1000
            while(not (20 < out_prob < 95) ):
                out_prob = normal(mean_completeness, stdev)
            genome2guass[g] = out_prob
            return out_prob

    @classmethod
    def run_boots(cls, name, nb_genomes, genome_len, ratio, accessory_len):
        mockmotu = MockmOTU(name, ratio, genome_len, nb_genomes, accessory_len, MockmOTU.guauss_completes)
