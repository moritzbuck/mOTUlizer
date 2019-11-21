import os
import pandas
from tqdm import tqdm
from os.path import join as pjoin
import sys
sys.path.append("/home/moritz/repos/moritz/0039_mOTUlizer/")
from mOTUlizer.classes.mOTU import mOTU

gtdb_bac_md_file = "/home/moritz/dbs/bac120_metadata_r89.tsv"
gtdb_ar_md_file = "/home/moritz/dbs/ar122_metadata_r89.tsv"
local_gtdb_path = "/home/moritz/proj_folder/uppstore2018126/moritz/gtdb_genomes/"

md = pandas.concat([pandas.read_csv(gtdb_bac_md_file, sep="\t"),pandas.read_csv(gtdb_ar_md_file, sep="\t")])
dd = { v[1]['accession'] : v[1]['gtdb_genome_representative']  for v in md.iterrows()}
mOTU2genome =  {d : []  for d in set(dd.values())}
checkm = { v[1]['accession'] : v[1]['checkm_completeness']  for v in md.iterrows()}


for k,v in tqdm(dd.items()):
    mOTU2genome[v] += [k]

def make_file(genome_id) :
    if genome_id.startswith("UBA") :
        path = pjoin(local_gtdb_path, "UBA/{gtdb_id}/proteom.faa".format(gtdb_id = genome_id))
    else :
        path = pjoin(local_gtdb_path, "{first}/{second}/{third}/{gtdb_id}/proteom.faa".format(first = genome_id.split('_')[-1][0:3],second = genome_id.split('_')[-1][3:6], third = genome_id.split('_')[-1][6:9], gtdb_id = genome_id))
    return path

mOTU2aa = {k : {g : make_file(g) for g in v} for k,v in mOTU2genome.items() if len(v) > 5}
mOTU2aa = {k : {g : p for g, p in v.items() if os.path.exists(p)} for k,v in tqdm(list(mOTU2aa.items())[0:10])}
mOTU2aa = {k : v for k,v in tqdm(mOTU2aa.items()) if len(v) > 5 }
for k,v in tqdm(mOTU2aa.items()):
     if len(v) > 50:
         selected = random.sample(list(v.items()),50)
         mOTU2aa[k] = dict(selected)


for k, v in mOTU2aa:
    motu = mOTU( k , v , cog_dict, checkm_dict = {v : checkm [v]for vv in v.keys()})
