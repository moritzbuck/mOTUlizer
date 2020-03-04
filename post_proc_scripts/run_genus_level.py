import shutil
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

md = pandas.concat([pandas.read_csv(gtdb_bac_md_file, sep="\t"),pandas.read_csv(gtdb_ar_md_file, sep="\t")])
dd = { v[1]['accession'] : v[1]['gtdb_taxonomy']  for v in md.iterrows()}
dd = {k : [vv[3:] for vv in v.split(";") if vv.startswith("g__")][0] for k,v in dd.items()}

mOTU2genome =  {d : []  for d in set(dd.values())}
checkm = { v[1]['accession'] : v[1]['checkm_completeness']  for v in md.iterrows()}
md.index = md.accession

anoxi_mOTUs = {motu : {a[:-4] : pjoin(local_anoxic_path, motu, a) for a in os.listdir(pjoin(local_anoxic_path, motu))} for motu in os.listdir(local_anoxic_path)}
anoxic_md = pandas.read_csv(anoxic_md_file)
anoxic_md.index = anoxic_md['Unnamed: 0']

checkm.update({ v[1]['Unnamed: 0'] : v[1]['completeness']  for v in anoxic_md.iterrows()})
anoxi_mOTUs = { k : {kk : vv for kk, vv in v.items() if anoxic_md.loc[kk,"contamination"] < 5 and anoxic_md.loc[kk,"completeness"] > 40} for k,v in anoxi_mOTUs.items() }
anoxi_mOTUs = {k : v for k,v in tqdm(anoxi_mOTUs.items()) if sum([ checkm[k] for k in v]) > 500}

with open("full_taxonomy.tax") as handle:
    taxonomy = {l.split(",")[0] : l[:-1].split(",")[1:] for l in handle}
taxonomy_local = { k : set([tuple(taxonomy[vv]) for vv in v if vv in taxonomy]) for k,v in anoxi_mOTUs.items()}
taxonomy_local.update({k : {vv for vv in v if vv.count('') == min([zz.count('') for zz in v])} for k, v in taxonomy_local.items() if len(v) > 1})
taxonomy_local.update({k : { max(v, key = lambda y : [tuple(taxonomy[zz]) for zz in anoxi_mOTUs[k]].count(y))} for k, v in taxonomy_local.items() if len(v) > 1})
taxonomy_local = {k : list(v)[0] for k,v in taxonomy_local.items()}

for k,v in tqdm(anoxi_mOTUs.items()):
     if len(v) > 50:
         selected = random.sample(list(v.items()),50)
         anoxi_mOTUs[k] = dict(selected)

for k,v in tqdm(dd.items()):
    mOTU2genome[v] += [k]

def make_file(genome_id) :
    if genome_id.startswith("UBA") :
        path = pjoin(local_gtdb_path, "UBA/{gtdb_id}/proteom.faa".format(gtdb_id = genome_id))
    else :
        path = pjoin(local_gtdb_path, "{first}/{second}/{third}/{gtdb_id}/proteom.faa".format(first = genome_id.split('_')[-1][0:3],second = genome_id.split('_')[-1][3:6], third = genome_id.split('_')[-1][6:9], gtdb_id = genome_id))
    return path

out_dat = {}
with open("gtdb_genus_cores.json") as handle:
    out_dat = json.load(handle)


mOTU2aa = {k : {g : make_file(g) for g in v} for k,v in mOTU2genome.items() if len(v) > 5 and k not in out_dat and len(set(md.loc[v].gtdb_taxonomy)) > 4}
mOTU2aa = {k : {g : p for g, p in v.items() if os.path.exists(p)} for k,v in tqdm(list(mOTU2aa.items()))}
mOTU2aa = {k : v for k,v in tqdm(mOTU2aa.items()) if len(v) > 5 }
for k,v in tqdm(mOTU2aa.items()):
     if len(v) > 50:
         selected = random.sample(list(v.items()),50)
         mOTU2aa[k] = dict(selected)

#anoxi_mOTUs = {k : v for k,v in anoxi_mOTUs.items() if k not in out_dat}


for i, t in tqdm(enumerate(mOTU2aa.items())):
    k, v = t
    if k not in out_dat:
        motu = mOTU( k , v , None, checkm_dict = {vv : checkm [vv] for vv in v.keys()})
        out_dat.update(motu.get_stats())
        if i % 10 == 0:
            print("check pointing")
            with open("gtdb_genus_cores.json", "w") as handle :
                json.dump(out_dat, handle)

with open("gtdb_genus_cores.json", "w") as handle :
    json.dump(out_dat, handle)

for i, t in tqdm(enumerate(anoxi_mOTUs.items())):
    k, v = t
    if k not in out_dat:
        motu = mOTU( k , v , None, checkm_dict = {vv : checkm [vv] for vv in v.keys()})
        out_dat.update(motu.get_stats())
        if i % 10 == 0:
            print("check pointing")
            with open("gtdb_genus_cores.json", "w") as handle :
                json.dump(out_dat, handle)

with open("gtdb_genus_cores.json", "w") as handle :
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
        "nb_species" : len(set(md.loc[v['completes'].keys()].gtdb_taxonomy)),
        "set" : "gtdb" if (k.startswith("UBA") or k.startswith("RS_") or k.startswith("GB_")) else "anoxi",
        "rep" : min(sorted(genome2comps(k,v).items(), key = lambda x : x[1][1])[0:5], key = lambda y : y[1][1])[0]
        }
        for k,v in out_dat.items()
    }, orient="index")


gtdb_tax = {k : {l : vv[3:] for l,vv in zip(["domain", "phylum", "class", "order", "family", "genus", "species"] , v.split(";"))}  for k,v in md.loc[cool_dat.rep]['gtdb_taxonomy'].iteritems() if not k.startswith("aniOTU_")}
gtdb_tax.update({k : {a : b for a,b in zip(["domain", "phylum", "class", "order", "family", "genus", "species"], v)} for k,v in taxonomy_local.items()})
taxo = pandas.DataFrame.from_dict(gtdb_tax, orient="index")
del taxo['species']
taxo = taxo.loc[cool_dat.rep]
taxo.index = cool_dat.index

cool_dat = pandas.concat([cool_dat, taxo], axis = 1)
cool_dat.to_csv("cool_genus_dat.csv")

[shutil.copy(make_file(k), "motulizer_reps/" + k + ".faa") for k in tqdm(cool_dat.loc[cool_dat.set != "anoxi", "rep"])]
[shutil.copy("proteomics/proteoms/" + k + ".faa", "motulizer_reps/") for k in tqdm(cool_dat.loc[cool_dat.set == "anoxi", "rep"])]

for f in tqdm(os.listdir("motulizer_reps/")):
    seqs = [ s for s in SeqIO.parse(pjoin("motulizer_reps/", f), "fasta")]
    for i,s in enumerate(seqs):
        s.id = f + "_" + str(i)
    SeqIO.write(seqs, pjoin("motulizer_reps/", f), "fasta")
