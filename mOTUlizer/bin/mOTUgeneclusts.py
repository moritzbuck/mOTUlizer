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

#tt = mOTU(name = "test", genomes = genomes, make_gene_clustering = True, thread=4, storage = "/home/moritz/temp/test_motusuite/")
#tt.export_gene_clusters(file = "test_stuff.json")
tt = mOTU(name = "test", genomes = genomes, storage = "/home/moritz/temp/test_motusuite/")
tt.load_gene_clusters("test_stuff.json")

for t in tt.gene_clusters:
    t.quiet = True
    _ = t.get_genes()
lazy_div = lambda a,b : None if b == 0 else a/b
gene_clusters = [c for c in tqdm(tt.gene_clusters) if len(c) > 3 and max([len(l[1]) for l in c._genome2genes.items()]) > 1 ]
# test_id = "ERR2098377.12_00548"
# print(tt[test_id.split("_")[0]].get_gene(test_id).translate())
# print(tt[test_id.split("_")[0]].get_gene(test_id))
# print(tt[test_id.split("_")[0]].get_aa(test_id))
print(gene_clusters[1].get_alignment())
[print(c.name, len(c)/len(c.get_genomes())) for c in gene_clusters]
