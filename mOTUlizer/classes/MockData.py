from random import random
from mOTUlizer.classes.mOTU import mOTU, mean
from collections import defaultdict
from tqdm import tqdm
from numpy.random import normal

genome2guass = {}


class MockmOTU(mOTU):
    def __repr__(self) :
        return "< MockmOTU with {n} genomes, of average {c}% completness, with core/genome_len of {r} >".format(c = 100*self.mean_completeness, n = len(self), r = self.ratio)

    def __init__(self, name, ratio, genome_len, nb_genomes, accessory_len, completeness):

        self.ratio = ratio
        core_len = int(genome_len*self.ratio)
        core = {"CoreTrait_{}".format(i) for i in range(core_len)}

        rank2prob = lambda i: 1/(i+1)+0.01
        accessory = [rank2prob(i) for i in range(accessory_len)]

        mock_genomes = dict()


        for k in tqdm(range(nb_genomes)):
            variable = []
            i = 1
            while(len(variable) < (genome_len*(1-ratio))):
                if (len(accessory)-1) < i :
                    variable += ["AccessoryTrait_{}".format(i)]
                    accessory += [rank2prob(i)]
                else :
                    if random() < accessory[i]:
                        variable += ["AccessoryTrait_{}".format(i)]
                i += 1
            mock_genomes["Genome_{}".format(k)] = core.union(variable)

            self.incompletes = {g : {vv for vv in v if random() < (completeness(g)/100)} for g, v in mock_genomes.items()}
            self.mean_completeness = mean([len({vv for vv in v if vv.startswith("CoreTrait_")})/core_len for c,v in self.incompletes.items()])
        self.accessory = accessory
        super().__init__(name = name, faas = {}, cog_dict = self.incompletes, checkm_dict = {g : 95 for g in self.incompletes}, max_it = 50)


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
