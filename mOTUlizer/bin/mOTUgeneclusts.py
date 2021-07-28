from mOTUlizer.classes.mOTU import mOTU
from mOTUlizer.classes.MetaBin import MetaBin
from mOTUlizer.classes.tools.Muscle import Muscle
from tqdm import tqdm
import os
fna_folder = "example_files/fnas/"
gff_folder = "example_files/gffs/"

fnas = { f[:-4] : fna_folder + f for f in os.listdir(fna_folder)}
gffs = { f[:-4] : gff_folder + f for f in os.listdir(gff_folder)}

genomes = []
i = 0
for g in fnas:
    test = MetaBin(g, nucleotide_file = fnas[g], gff_file = gffs[g])
    _ = test.get_amino_acids()
    genomes.append(test)

tt = mOTU(name = "test", genomes = genomes, make_gene_clustering = True, thread=4)
tt.export_gene_clusters(file = "test_stuff.json")
#tt = mOTU(name = "test", genomes = genomes, storage = "/home/moritz/temp/test_motusuite/")
#tt.load_gene_clusters("test_stuff.json")

for t in tt.gene_clusters:
    t.quiet = True
gene_clusters = [c for c in tqdm(tt.gene_clusters)]
print(gene_clusters[0].get_alignment())
print(tt['ERR599236.12'].get_gene('ERR599236.12_00324').translate())
print(tt['ERR599236.12'].nucleotide_file)
