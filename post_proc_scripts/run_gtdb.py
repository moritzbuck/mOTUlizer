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
local_anoxic_path = "/home/moritz/uppmax/people/0023_anoxicencyclo/4500_assembly_analysis/mags/mOTUs/"
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

aa2counts = {}
aa2norm_counts = {}
aa2def = {}
path_format_1 = "/home/moritz/uppmax/people/0023_anoxicencyclo/1000_processed_reads/{sample}/assemblies/megahit/binning/metabat/clean_bins/{bin}/{bin}.ffn"
path_format_2 = "/home/moritz/uppmax/people/0023_anoxicencyclo/1500_coasses/{sample}/assemblies/megahit/binning/metabat/clean_bins/{bin}/{bin}.ffn"

for k, v in tqdm(anoxi_mOTUs.items()):
    for vv in v.values():
        bin_name = os.path.basename(vv)[:-4]
        sample = os.path.basename(vv).split("_")[0]
        ffn_name = path_format_1.format(bin = bin_name, sample = sample) if os.path.exists(path_format_1.format(bin = bin_name, sample = sample)) else path_format_2.format(bin = bin_name, sample = sample)
        for s in SeqIO.parse(ffn_name, "fasta"):
            ii = s.id
            aa2def[ii] = " ".join(s.description.split()[1:])
            aa2counts[ii] = {'A' : s.seq[1::3].count('A'), 'T' : s.seq[1::3].count('T'), 'G' : s.seq[1::3].count('G'), 'C' : s.seq[1::3].count('C')}
            aa2norm_counts[ii] = {'A' : s.seq.count('A'), 'T' : s.seq.count('T'), 'G' : s.seq.count('G'), 'C' : s.seq.count('C')}



for k,v in tqdm(dd.items()):
    mOTU2genome[v] += [k]

def get_GC_pair(out_d, count_dict = aa2counts):
    core =  out_d['core']
    aux_sinsingle = set(out_d['aux_genome']).difference(out_d['singleton_cogs'])

    counts = { 'core' : 0, 'aux_sinsingle' : 0 , 'singeltons' :0}
    gc_counts = { 'core' : 0, 'aux_sinsingle' : 0 , 'singeltons' :0}
    for k, v in tqdm(out_d['cogs']['aa'].items()):
        content = count_dict.get(k, { 'A' : 0, 'T' : 0, 'G' : 0, 'C' : 0 })
        if v in core:
            counts['core'] += sum(content.values())
            gc_counts['core'] += content['G'] + content['C']
        elif v in aux_sinsingle:
            counts['aux_sinsingle'] += sum(content.values())
            gc_counts['aux_sinsingle'] += content['G'] + content['C']
        else :
            counts['singeltons'] += sum(content.values())
            gc_counts['singeltons'] += content['G'] + content['C']

    if sum(counts.values()) > 0:
        return {k : gc_counts[k]/counts[k] if counts[k] > 0 else None for k in counts.keys()}
    else :
        return None

def is_CRISPRed(out_d, word = "CRISPR"):
    core =  out_d['core']
    aux = out_d['aux_genome']

    crispr_genes = [k for aa,k in out_d['cogs']['aa'].items() if word in aa2def.get(aa, "blabla")]
    if all([not k in aa2def for k in out_d['cogs']['aa'].keys()]):
        return { 'core_' + word : None, 'aux_' + word : None}
    counts = {  'core_' + word : len(set(crispr_genes).intersection(core) ), 'aux_' + word : len(set(crispr_genes).intersection(aux))}
    return counts


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

anoxi_mOTUs = {k : v for k,v in anoxi_mOTUs.items() if k not in out_dat}


#out_dat = {}
for i, t in tqdm(enumerate(mOTU2aa.items())):
    k, v = t
    if k not in out_dat:
        motu = mOTU( k , v , None, checkm_dict = {vv : checkm [vv] for vv in v.keys()})
        out_dat.update(motu.get_stats())
        if i % 10 == 0:
            print("check pointing")
            with open("gtdb_cores.json", "w") as handle :
                json.dump(out_dat, handle)
with open("gtdb_cores.json", "w") as handle :
    json.dump(out_dat, handle)

for i, t in tqdm(enumerate(anoxi_mOTUs.items())):
    k, v = t
    if k not in out_dat:
        motu = mOTU( k , v , None, checkm_dict = {vv : checkm [vv] for vv in v.keys()})
        out_dat.update(motu.get_stats())
        if i % 10 == 0:
            print("check pointing")
            with open("gtdb_cores.json", "w") as handle :
                json.dump(out_dat, handle)
with open("gtdb_cores.json", "w") as handle :
    json.dump(out_dat, handle)


genome2comps = lambda k, v: {genome : (len(set(v['core']).intersection(cogs))/len(v['core']), len(set(v['core']).difference(cogs))/len(v['core'])) for genome, cogs in v['cogs']['genome'].items() }


gc3_pairs = {k : get_GC_pair(v) for k,v in tqdm(out_dat.items())}
tt3 = pandas.DataFrame.from_dict({k : {'core': None, 'aux_sinsingle' : None, 'singeltons' : None} if not v else v for k, v in gc3_pairs.items()}, orient = "index")
tt3.columns = ['core_gc3', 'aux_sinsingle_gc3', 'singletons_gc3']

gc_pairs = {k : get_GC_pair(v, aa2norm_counts) for k,v in tqdm(out_dat.items())}
tt2 = pandas.DataFrame.from_dict({k : {'core': None, 'aux_sinsingle' : None, 'singeltons' : None} if not v else v for k, v in gc_pairs.items()}, orient = "index")
tt2.columns = ['core_gc', 'aux_sinsingle_gc', 'singletons_gc']

crispr = {k : is_CRISPRed(v) for k,v in tqdm(out_dat.items())}
cricri = pandas.DataFrame.from_dict(crispr, orient = "index")

nucleases = {k : is_CRISPRed(v, "nucleas") for k,v in tqdm(out_dat.items())}
nucnuc = pandas.DataFrame.from_dict(nucleases, orient = "index")


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



gtdb_tax = {k : {l : vv[3:] for l,vv in zip(["domain", "phylum", "class", "order", "family", "genus", "species"] , v.split(";"))}  for k,v in md.loc[cool_dat.index]['gtdb_taxonomy'].iteritems() if not k.startswith("aniOTU_")}
gtdb_tax.update({k : {a : b for a,b in zip(["domain", "phylum", "class", "order", "family", "genus", "species"], v)} for k,v in taxonomy_local.items()})
taxo = pandas.DataFrame.from_dict(gtdb_tax, orient="index")

cool_dat = pandas.concat([cool_dat, taxo, tt2, tt3, cricri, nucnuc], axis = 1)
cool_dat.to_csv("cool_dat.csv")

anv_data = cool_dat.copy()
anv_data['ID'] = anv_data['rep']
anv_data['variable_fraction'] = anv_data['mean_variable']/anv_data['mean_cogs']
anv_data['aux_genome_div'] = anv_data['aux_sinsingle']/anv_data['mean_variable']/sqrt(anv_data['nb_genomes'])
del anv_data['rep']

gc_contents = md.gc_percentage.to_dict()
gc_contents.update((anoxic_md.GC*100).to_dict())
anv_data['gc_content'] = [gc_contents[v] for v in anv_data.ID]

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
