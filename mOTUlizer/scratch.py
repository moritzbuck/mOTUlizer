import os
import pandas
import shutil
from numpy import mean
from Bio import SeqIO



stats = pandas.read_csv("test_data/magstats.csv", index_col = 0)

with open("test_data/taxonomy.tax") as handle:
    taxonomy = {l.split(",")[0] : l[:-1].split(",")[1:] for l in handle}

with open("test_data/mag2cog.tsv") as handle:
    mag2cog = {l.split()[0] : l[:-1].split()[1:] for l in handle}


    
class MetaBin:
    def __repr__(self) :
        return "< {tax} bin {name} >".format(tax = self.tail(), name = self.name)
    
    
    def __init__(self, name):
        self.name = name
        self.taxonomy = taxonomy[self.name]
        self.cogs = set(mag2cog[self.name])
        self.gc = stats.loc[self.name, "GC"]
        self.coding_density = stats.loc[self.name, "coding_density"]
        self.genome_length = stats.loc[self.name, "length"]
        self.nb_contig = stats.loc[self.name, "nb_contigs"]
        self.nb_prots = stats.loc[self.name, "nb_proteins"]
        self.checkm_complet = stats.loc[self.name, "completeness"]
        self.checkm_contam = stats.loc[self.name, "contamination"]
        self.checkm_hetero = stats.loc[self.name, "strain_heterogeneity"]
        
    def tail(self) :
        return [f for f in self.taxonomy if f != ''][-1]

    def overlap(self, target):
        return self.cogs.intersection(target.cogs)
    
class mOTU:
    def __len__(self):
        return len(self.members)

    def __repr__(self):
        return "< {tax} mOTU {name}, of {len} members >".format(name = self.name, len = len(self), tax = self.consensus_tax()[0].split(";")[-1])
    
    def __init__(self, name, members):
        self.name = name
        self.members = [MetaBin(m) for m in members]

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
        return mean(list(self.overlap_dict.values()))
    
otu_list = []
with open("test_data/mOTUs.txt") as handle:
    for l in handle:
        name = l.split()[0]
        bins = l.split()[1].split(";")
        otu_list += [ mOTU( name = name, members = bins ) ]
        
