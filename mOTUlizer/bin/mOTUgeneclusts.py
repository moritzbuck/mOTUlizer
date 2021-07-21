from mOTUlizer.classes.mOTU import mOTU
import os
fna_folder = "example_files/fnas/"
gff_folder = "example_files/gffs/"

fnas = { f[:-4] : fna_folder + f for f in os.listdir(fna_folder)}
gffs = { f[:-4] : gff_folder + f for f in os.listdir(gff_folder)}


for g in fnas:
    genomes = [MetaBin(g, fnas = fnas[g], gffs = gffs[g]) for bin_name in self.gene_clusters_dict.keys()]

tt = mOTU(name = "test", genomes = genomes)
