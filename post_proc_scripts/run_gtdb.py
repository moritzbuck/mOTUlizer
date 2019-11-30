import os
import pandas
from tqdm import tqdm
from os.path import join as pjoin
import sys
import random
import json
from numpy import mean


sys.path.append("/home/moritz/repos/moritz/0039_mOTUlizer/")
from mOTUlizer.classes.mOTU import mOTU

gtdb_bac_md_file = "/home/moritz/dbs/bac120_metadata_r89.tsv"
gtdb_ar_md_file = "/home/moritz/dbs/ar122_metadata_r89.tsv"
local_gtdb_path = "/home/moritz/proj_folder/uppstore2018126/moritz/gtdb_genomes/"
local_anoxic_path = "/home/moritz/people/0023_anoxicencyclo/4500_assembly_analysis/mags/mOTUs/"
anoxic_md_file = "/home/moritz/people/0023_anoxicencyclo/4500_assembly_analysis/magstats.csv"


anoxi_mOTUs = {motu : {a[:-4] : pjoin(local_anoxic_path, motu, a) for a in os.listdir(pjoin(local_anoxic_path, motu))} for motu in os.listdir(local_anoxic_path)}
anoxic_md = pandas.read_csv(anoxic_md_file)
anoxi_mOTUs = {k : v for k,v in tqdm(anoxi_mOTUs.items()) if len(v) > 5 }
for k,v in tqdm(anoxi_mOTUs.items()):
     if len(v) > 50:
         selected = random.sample(list(v.items()),50)
         anoxi_mOTUs[k] = dict(selected)


md = pandas.concat([pandas.read_csv(gtdb_bac_md_file, sep="\t"),pandas.read_csv(gtdb_ar_md_file, sep="\t")])
dd = { v[1]['accession'] : v[1]['gtdb_genome_representative']  for v in md.iterrows()}
mOTU2genome =  {d : []  for d in set(dd.values())}
checkm = { v[1]['accessio'] : v[1]['checkm_completeness']  for v in md.iterrows()}
checkm.update({ v[1]['Unnamed: 0'] : v[1]['completeness']  for v in anoxic_md.iterrows()})
md.index = md.accession

for k,v in tqdm(dd.items()):
    mOTU2genome[v] += [k]

def make_file(genome_id) :
    if genome_id.startswith("UBA") :
        path = pjoin(local_gtdb_path, "UBA/{gtdb_id}/proteom.faa".format(gtdb_id = genome_id))
    else :
        path = pjoin(local_gtdb_path, "{first}/{second}/{third}/{gtdb_id}/proteom.faa".format(first = genome_id.split('_')[-1][0:3],second = genome_id.split('_')[-1][3:6], third = genome_id.split('_')[-1][6:9], gtdb_id = genome_id))
    return path

with open("gtdb_cores.json") as handle:
    out_dat = json.load(handle)


mOTU2aa = {k : {g : make_file(g) for g in v} for k,v in mOTU2genome.items() if len(v) > 5 and k not in out_dat}
mOTU2aa = {k : {g : p for g, p in v.items() if os.path.exists(p)} for k,v in tqdm(list(mOTU2aa.items()))}
mOTU2aa = {k : v for k,v in tqdm(mOTU2aa.items()) if len(v) > 5 }
for k,v in tqdm(mOTU2aa.items()):
     if len(v) > 50:
         selected = random.sample(list(v.items()),50)
         mOTU2aa[k] = dict(selected)


#out_dat = {}
for k, v in tqdm(mOTU2aa.items()):
    if k not in out_dat:
        motu = mOTU( k , v , None, checkm_dict = {vv : checkm [vv] for vv in v.keys()})
        out_dat.update(motu.get_stats())
        with open("gtdb_cores.json", "w") as handle :
            json.dump(out_dat, handle)

for k, v in tqdm(anoxi_mOTUs.items()):
    if k not in out_dat:
        motu = mOTU( k , v , None, checkm_dict = {vv : checkm [vv] for vv in v.keys()})
        out_dat.update(motu.get_stats())
        with open("gtdb_cores.json", "w") as handle :
            json.dump(out_dat, handle)



cool_dat =  pandas.DataFrame.from_dict(
    { k :
        {
        "mean_cogs" : mean([100*len(l)/v['completes'][kk] for kk,l in  v['cogs']['genome'].items()]),
        "core_len" : len(v['core']),
        "aux_len" : len(v['aux_genome']),
        "mean_single" : mean([len(set(l).intersection(v['singleton_cogs'])) for kk,l in  v['cogs']['genome'].items()]),
        "aux_sinsingle" : len(set(v['aux_genome']).difference(v['singleton_cogs'])),
        "nb_genomes" : len(v['cogs']['genome']),
        "set" : "gtdb" if (k.startswith("UBA") or k.startswith("RS_") or k.startswith("GB_")) else "anoxi"
        }
        for k,v in out_dat.items()
    }, orient="index")

taxo = pandas.DataFrame.from_dict({k : {l : vv[3:] for l,vv in zip(["domain", "phylum", "class", "order", "family", "genus", "species"] , v.split(";"))}  for k,v in md.loc[cool_dat.index]['gtdb_taxonomy'].iteritems()}, orient="index")

cool_dat = pandas.concat([cool_dat, taxo], axis = 1)
cool_dat.to_csv("cool_dat.csv")
