from os.path import join as pjoin
import sys

class MetaBin:
    def __repr__(self) :
        return "< bin {name} with {n} cogs>".format(n = len(self.cogs) if cogs else "NA", name = self.name)

    def __init__(self, name, cogs,fnas, faas, complet, contamin = 0):
        self.name = name
        self.cogs = cogs
        self.faas = faas
        self.checkm_complet = complet
        self.checkm_contamin = contamin
        self.checkm_complet = self.checkm_complet if self.checkm_complet < 95 else 95
        self.new_completness = None

    def overlap(self, target):
        return self.cogs.intersection(target.cogs)

    def estimate_nb_cogs(self):
        assert self.new_completness != None, "new_completness not computed, please do"
        return 100*len(self.cogs)/self.new_completness
