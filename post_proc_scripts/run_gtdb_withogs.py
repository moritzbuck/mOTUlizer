import shutil
import os
import pandas
from tqdm import tqdm
from os.path import join as pjoin
import sys
import random
import json
from numpy import mean, sqrt
from ete3 import Tree
from Bio import SeqIO


sys.path.append("/home/moritz/projects/0039_mOTUlizer/")
from mOTUlizer.classes.mOTU import mOTU

gtdb_bac_md_file = "/home/moritz/uppmax/dbs/bac120_metadata_r89.tsv"
gtdb_ar_md_file = "/home/moritz/uppmax/dbs/ar122_metadata_r89.tsv"
local_gtdb_path = "/home/moritz/uppmax/proj_folder/uppstore2018126/moritz/gtdb_genomes/"
local_anoxic_path = "/home/moritz/data/data_submit/mOTUs/"
anoxic_md_file = "/home/moritz/uppmax/people/0023_anoxicencyclo/4500_assembly_analysis/magstats.csv"

md = pandas.concat([pandas.read_csv(gtdb_bac_md_file, sep="\t"),pandas.read_csv(gtdb_ar_md_file, sep="\t")])
dd = { v[1]['accession'] : v[1]['gtdb_genome_representative']  for v in md.iterrows()}
mOTU2genome =  {d : []  for d in set(dd.values())}
checkm = { v[1]['accession'] : v[1]['checkm_completeness']  for v in md.iterrows()}
md.index = md.accession

anoxi_mOTUs = {motu : {a[:-4] : pjoin(local_anoxic_path, motu, a) for a in os.listdir(pjoin(local_anoxic_path, motu))} for motu in os.listdir(local_anoxic_path)}
anoxic_md = pandas.read_csv(anoxic_md_file)
anoxic_md.index = anoxic_md['Unnamed: 0']

checkm.update({ v[1]['Unnamed: 0'] : v[1]['completeness']  for v in anoxic_md.iterrows()})
anoxi_mOTUs = { k : {kk : vv for kk, vv in v.items() if anoxic_md.loc[kk,"contamination"] < 5 and anoxic_md.loc[kk,"completeness"] > 40} for k,v in anoxi_mOTUs.items() }
anoxi_mOTUs = {k : v for k,v in tqdm(anoxi_mOTUs.items()) if sum([ checkm[k] for k in v]) > 500}

with open("/home/moritz/data/data_submit/full_taxonomy.tax") as handle:
    taxonomy = {l.split(",")[0] : l[:-1].split(",")[1:] for l in handle}
taxonomy_local = { k : set([tuple(taxonomy[vv]) for vv in v if vv in taxonomy]) for k,v in anoxi_mOTUs.items()}
taxonomy_local.update({k : {vv for vv in v if vv.count('') == min([zz.count('') for zz in v])} for k, v in taxonomy_local.items() if len(v) > 1})
taxonomy_local.update({k : { max(v, key = lambda y : [tuple(taxonomy[zz]) for zz in anoxi_mOTUs[k]].count(y))} for k, v in taxonomy_local.items() if len(v) > 1})
taxonomy_local = {k : list(v)[0] for k,v in taxonomy_local.items()}

for k,v in tqdm(dd.items()):
    mOTU2genome[v] += [k]


with open("/home/moritz/projects/0023_anoxicencyclo/data/bin2og.json") as handle:
    bin2og = json.load(handle)
bin2og = { k : set(v) for k,v in bin2og.items()}

out_dat = {}

for i, t in tqdm(enumerate(anoxi_mOTUs.items())):
    k, v = t
    motu = mOTU( k , v , {vv : bin2og [vv] for vv in v.keys()}, checkm_dict = {vv : checkm [vv] for vv in v.keys()})
    out_dat.update(motu.get_stats())
    if i % 10 == 0:
        print("check pointing")
        with open("og_cores.json", "w") as handle :
            json.dump(out_dat, handle)
with open("og_cores.json", "w") as handle :
    json.dump(out_dat, handle)


genome2comps = lambda k, v: {genome : (len(set(v['core']).intersection(cogs))/len(v['core']), len(set(v['core']).difference(cogs))/len(v['core'])) for genome, cogs in v['cogs']['genome'].items() }


cool_dat =  pandas.DataFrame.from_dict(
    { k :
        {
        "mean_cogs" : mean([100*len(l)/v['completes'][kk] for kk,l in  v['cogs']['genome'].items()]),
        "core_len" : len(v['core']),
        "aux_len" : len(v['aux_genome']),
        "aux_sinsingle" : len(set(v['aux_genome']).difference(v['singleton_cogs'])),
        "mean_single" : mean([len(set(l).intersection(v['singleton_cogs'])) for kk,l in  v['cogs']['genome'].items()]),
        "mean_variable" : mean([100*(len(l)-len(set(l).intersection(v['core']))-len(set(l).intersection(v['singleton_cogs'])))/v['completes'][kk] for kk,l in  v['cogs']['genome'].items()]),
        "nb_genomes" : len(v['cogs']['genome']),
        "set" : "gtdb" if (k.startswith("UBA") or k.startswith("RS_") or k.startswith("GB_")) else "anoxi",
        "rep" : min(sorted(genome2comps(k,v).items(), key = lambda x : x[1][1])[0:5], key = lambda y : y[1][1])[0]
        }
        for k,v in out_dat.items()
    }, orient="index")



gtdb_tax = {k : {a : b for a,b in zip(["domain", "phylum", "class", "order", "family", "genus", "species"], v)} for k,v in taxonomy_local.items()}
taxo = pandas.DataFrame.from_dict(gtdb_tax, orient="index")

cool_dat = pandas.concat([cool_dat, taxo], axis = 1)
cool_dat.to_csv("cool_dat.csv")

anv_data = cool_dat.copy()
anv_data['ID'] = anv_data['rep']
anv_data['variable_fraction'] = anv_data['mean_variable']/anv_data['mean_cogs']
anv_data['aux_genome_div'] = anv_data['aux_sinsingle']/anv_data['mean_variable']/sqrt(anv_data['nb_genomes'])
del anv_data['rep']


tree = Tree("motulizer_reps.tree.nwk")
nodes = list(tree.iter_leaf_names())
anv_data = anv_data.loc[[v in nodes for v in anv_data.ID ]]
tt = anv_data['ID']
del anv_data['ID']
anv_data.insert(0, 'ID', tt)

anv_data.to_csv("for_anvio.tsv", sep="\t", index=None)



[shutil.copy(make_file(k), "motulizer_reps/" + k + ".faa") for k in tqdm(cool_dat.loc[cool_dat.set != "anoxi", "rep"])]
[shutil.copy("proteomics/proteoms/" + k + ".faa", "motulizer_reps/") for k in tqdm(cool_dat.loc[cool_dat.set == "anoxi", "rep"])]

for f in tqdm(os.listdir("motulizer_reps/")):
    seqs = [ s for s in SeqIO.parse(pjoin("motulizer_reps/", f), "fasta")]
    for i,s in enumerate(seqs):
        s.id = f + "_" + str(i)
    SeqIO.write(seqs, pjoin("motulizer_reps/", f), "fasta")
