from mOTUlizer.classes.mOTU import mOTU
from mOTUlizer.classes.MetaBin import MetaBin

import os
fna_folder = "example_files/fnas_smol/"
gff_folder = "example_files/gffs/"

fnas = { f[:-4] : fna_folder + f for f in os.listdir(fna_folder)}
gffs = { f[:-4] : gff_folder + f for f in os.listdir(gff_folder)}

genomes = []
i = 0
for g in fnas:
    test = MetaBin(g, amino_acid_file = fnas[g], gff_file = gffs[g])
    _ = test.get_amino_acids()
    genomes.append(test)

tt = mOTU(name = "test", genomes = genomes)
print(tt.gene_clustering)
