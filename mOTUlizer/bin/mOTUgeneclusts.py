from mOTUlizer.classes.mOTU import mOTU
from mOTUlizer.classes.MetaBin import MetaBin
from mOTUlizer.classes.tools.Muscle import Muscle
from tqdm import tqdm
import os
import time
from numpy import mean
import pandas
import multiprocessing as mp
print("Number of processors: ", mp.cpu_count())


folder = "example_files/aquadb_mOTU_00113/"
ffolder = "/home/moritz/temp/pnecs/fnas/"
gfolder = "/home/moritz/temp/pnecs/gffs/"

fnas = { f[:-4] : ffolder + f for f in os.listdir(ffolder) if f.endswith(".fna")}
gffs = { f[:-4] : gfolder + f for f in os.listdir(gfolder) if f.endswith(".gff")}

genomes = []
i = 0
for g in tqdm(fnas):
    test = MetaBin(g, nucleotide_file = fnas[g], gff_file = gffs[g])
#    _ = test.get_amino_acids()
    genomes.append(test)

#tt = mOTU(name = "pnecs", genomes = genomes, make_gene_clustering = True, thread=20, storage = "/home/moritz/temp/pnecs/motusuite_data/")
#tt.export_gene_clusters(file = "/home/moritz/temp/pnecs/motusuite_data/gene_clusters.json")
tt = mOTU(name = "test", genomes = genomes, storage = "/home/moritz/temp/test_motusuite/")
tt.load_gene_clusters("/home/moritz/temp/pnecs/motusuite_data/gene_clusters.json")

for t in tqdm(tt.gene_clusters):
    if len(t) > 1:
        _ = t.within_cluster_mutations()
lazy_div = lambda a,b : None if b == 0 else a/b

non_single_genes = [t for t in tt.gene_clusters if len(t) > 1]
genomes = [g for g in tt]

dat = {frozenset({genome1.name, genome2.name}) : genome1.gene_cluster_ani(genome2) for genome1 in tqdm(genomes) for genome2 in tqdm(genomes) if genome1 != genome2}
tt.compute_core()

gc_data = {}
for g in tqdm(non_single_genes):
    gc_data[g.name] = {}
    gc_data[g.name]['prevalence'] = len(g.get_genomes())/len(tt)
    gc_data[g.name]['core'] = g in tt.core
    mut_data = g.within_cluster_mutations()
    gc_data[g.name]['mean_dNdSp1'] = mean([ v['synonymous_muts']/(v['non_synonymous_muts']+1) for v in mut_data.values()])
    gc_data[g.name]['mean_ANI'] = mean([(v['synonymous_muts']+v['non_synonymous_muts'])/v['counted_bases'] for v in mut_data.values()])


anis = tt.get_anis(threads = 24)

anis = {frozenset((k[0].name, k[1].name)) : v['ani'] for k, v in anis.items() if k[1] != k[0]}

genome2mOTU = {g.name : motu.name if len(motu) > 1 else "singleton" for motu in clstrs for g in motu}

for k,v in dat.items():
    v['fastANI'] = anis.get(k, -1)/100
    k = tuple(k)
    v['estimated_nb_pairs'] = tt[k[0]].new_completness*tt[k[0]].new_completness*(len(tt[k[0]].gene_clusters)+len(tt[k[1]].gene_clusters))/10000/2
    if genome2mOTU[k[0]] == genome2mOTU[k[1]]:
        v['mOTU'] = genome2mOTU[k[0]]
    else :
        v['mOTU'] = (genome2mOTU[k[0]] + ";" +  genome2mOTU[k[1]]) if genome2mOTU[k[0]]  < genome2mOTU[k[1]] else  (genome2mOTU[k[1]] + ";" +  genome2mOTU[k[0]])


bmft = pandas.DataFrame.from_dict(dat, orient="index")
bmft.index = [";".join(t) for t in bmft.index]
bmft.to_csv("~/temp/pnecs.csv", index_label='pairs')



#print(dat)

#[c for c in tqdm(tt.gene_clusters) if len(c) > 3 and max([len(l[1]) for l in c._genome2genes.items()]) > 1 ]
#trees = [c.tree() for c in tqdm(tt.gene_clusters) if len(c) > 3 and max([len(l[1]) for l in c._genome2genes.items()]) > 1 ]

# test_id = "ERR2098377.12_00548"
# print(tt[test_id.split("_")[0]].get_gene(test_id).translate())
# print(tt[test_id.split("_")[0]].get_gene(test_id))
# print(tt[test_id.split("_")[0]].get_aa(test_id))
#print(gene_clusters[1].get_alignment())
#[print(c.name, len(c)/len(c.get_genomes())) for c in gene_clusters if len(c)/len(c.get_genomes()) > 1.5]
