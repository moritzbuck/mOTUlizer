from os.path import join as pjoin

class MetaBin:
    def __repr__(self) :
        return "< {tax} bin {name} >".format(tax = self.tail(), name = self.name)

    def __init__(self, name, taxonomy, stats, mag2cog, base_folder):
        self.name = name
        self.taxonomy = taxonomy[self.name]
        self.cogs = set(mag2cog[self.name])
        self.gc = stats.loc[self.name, "GC"]
        self.coding_density = stats.loc[self.name, "coding_density"]
        self.genome_length = stats.loc[self.name, "length"]
        self.nb_contig = stats.loc[self.name, "nb_contigs"]
        self.nb_prots = stats.loc[self.name, "nb_proteins"]
        self.checkm_complet = stats.loc[self.name, "completeness"]
        self.checkm_complet = self.checkm_complet if self.checkm_complet < 95 else 95
        self.checkm_contam = stats.loc[self.name, "contamination"]
        self.checkm_hetero = stats.loc[self.name, "strain_heterogeneity"]
        self.proteom = pjoin(base_folder, "proteoms", name + ".faa")
        self.genome = pjoin(base_folder, "genomes", name + ".fna")
        self.gff = pjoin(base_folder, "gffs", name + ".gff")
        self.new_completness = None

    def tail(self) :
        return [f for f in self.taxonomy if f != ''][-1]

    def overlap(self, target):
        return self.cogs.intersection(target.cogs)

    def estimate_nb_cogs(self):

        assert self.new_completness != None, "new_completness not computed, please do"
        return 100*len(self.cogs)/self.new_completness
