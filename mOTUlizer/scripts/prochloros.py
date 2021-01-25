from os.path import join as pjoin
import sys, os
from Bio import Entrez, SeqIO
from xml.etree import ElementTree as ET
from tqdm import tqdm
import re
from collections import defaultdict
import gzip
import json
from random import choices
from mOTUlizer import __version__
from mOTUlizer.classes import *
from mOTUlizer.utils import *
from mOTUlizer.classes.mOTU import mOTU
from multiprocessing import Pool
from statistics import mean
import h5py
import hdf5plugin
import pandas
import shutil
from subprocess import call

folder = "/home/moritz/projects/0039_mOTUlizer/test_data/prochlos/"
os.chdir(folder)
r95 = pandas.read_csv("/home/moritz/data/gtdb/bac120_metadata_r95.tsv", sep ="\t", index_col=0)
counts = {s : 0 for s in set(r95.loc[r95.scaffold_count < 4].gtdb_taxonomy)}
for s in r95.loc[r95.scaffold_count < 4].gtdb_taxonomy:
     counts[s] += 1
gids = {s : [] for s,v in counts.items() if v > 20}

for k,v in r95.loc[r95.scaffold_count < 4].gtdb_taxonomy.items():
     if v in gids:
         gids[v] += [k]

def make_folder(gid, k):
    gid = gid[3:]
    pat = pjoin("/home/moritz/data/gtdb/root/{k}/{gid}/{gid}.fna.gz".format(k = k.replace(";","/"),gid = gid))
    return pat

gids = {k: [l for l in ll if os.path.exists(make_folder(l, k))] for k, ll in gids.items() if sum([os.path.exists(make_folder(l, k)) for l in ll])  > 20}

def process_species(liss, tax):

    paths = [make_folder(g, tax) for g in liss]
    paths = [g for g in paths if os.path.exists(g)]
    os.makedirs("temp", exist_ok = True)
    for p in paths:
        shutil.copy(p, "temp/")
    paths = [pjoin("temp", g) for g in os.listdir("temp") ]
    for p in paths:
        call("unpigz " + p , shell=True)

    with open("temp_gff.list", "w") as handle:
        handle.writelines([os.path.basename(c)[:-7] + "\t" + c[:-3] + "\n" for c in paths])

    call("""
    ppanggolin annotate --fasta temp_gff.list  -o temp/ --basename temp --use_pseudo --tmpdir temp/ -c 20 -f
    ppanggolin cluster -p temp/temp.h5 -c 20 --tmpdir temp
    ppanggolin graph -p temp/temp.h5 -c 20
    ppanggolin partition -K3 -f  -p temp/temp.h5 --cpu 20

    python ../../mOTUlizer/bin/mOTUconvert.py --in_type ppanggolin temp/temp.h5 > temp/ppanggolin.gid2cog
    """, shell = True)
    call("""
    roary -p 22 -o other_soft/roary/s__Prochlorococcus_A_clusters.txt  -cd 1 -v `cat nucleotides/gff_list  | cut -f2`
    cp accessory* blast_identity_frequency.Rtab  core_accessory* gene_presence_absence.* number_of_* summary_statistics.txt other_soft/roary/s__Prochlorococcus_A
    python ../../mOTUlizer/bin/mOTUconvert.py --in_type roary other_soft/roary/s__Prochlorococcus_A_clusters.txt > static_data/roary_s__Prochlorococcus_A.gid2cog

    """,
    shell = True)
    shutil.rmtree("temp/")


def setup_data():
    os.makedirs(pjoin(folder, "nucleotides", "fnas"), exist_ok = True)
    os.makedirs(pjoin(folder, "proteins", "faas"), exist_ok = True)
    os.makedirs(pjoin(folder, "temp"), exist_ok = True)

    Entrez.email = username
    query = " OR ".join(ids)
    query = '"Prochlorococcus"[Organism]'
    ass_ids = Entrez.read(Entrez.esearch(db="assembly", term = query, retmax=1000000))['IdList']

    def get_assemblies(ass_id):
        ass_info = Entrez.read(Entrez.esummary(db="assembly", id = ass_id))
        genbank = ass_info['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
        refseq = ass_info['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']

        ftp_fold = genbank if refseq == "" else refseq
        ftp_file = "{fold}/{ID}_genomic.gbff.gz".format(fold = ftp_fold, ID = ftp_fold.split("/")[-1])
        call("wget {path} 2> /dev/null".format(path = ftp_file), shell = True)
        shutil.move(ftp_fold.split("/")[-1] + "_genomic.gbff.gz", "nucleotides/ncbi_gbffs")

    for ass_id in tqdm(ass_ids):
        get_assemblies(ass_id)

    genomes = dict()
    for f in tqdm(os.listdir("nucleotides/ncbi_gbffs")):
        with gzip.open("nucleotides/ncbi_gbffs" + f, "rt") as handle:
            genomes["_".join(f.split("_")[0:2])] = [s for s in SeqIO.parse(handle, "genbank")]

    def make_rec(feat, genome, seq):
        entry = ">{locus}:{genome}{ID} {product}\n{seq}\n"
        if 'translation' in feat.qualifiers:
            params = {
                'locus' : feat.qualifiers['locus_tag'][0],
                'genome' : genome,
                'ID' : (":" + feat.qualifiers['protein_id'][0]) if 'protein_id' in feat.qualifiers else "",
                'product' : feat.qualifiers['product'][0],
                'seq' : re.sub("(.{64})", "\\1\n", feat.qualifiers['translation'][0], 0, re.DOTALL)
            }
        else :
            seq = feat.extract(seq).translate(table = feat.qualifiers['transl_table'][0])
            params = {
                'locus' : feat.qualifiers['locus_tag'][0],
                'genome' : genome,
                'ID' : (":" + feat.qualifiers['protein_id'][0]) if 'protein_id' in feat.qualifiers else "",
                'product' : feat.qualifiers['product'][0],
                'seq' : re.sub("(.{64})", "\\1\n", str(seq.seq), 0, re.DOTALL)
            }
        return entry.format(**params)

    for k, v in tqdm(genomes.items()):
        SeqIO.write(v, pjoin(folder, "nucleotides", "fnas", k + ".fna"), "fasta")
        aas = [make_rec(s, k,vv) for vv in v for s in vv.features if s.type == "CDS"]
        if len(aas) > 0 :
            with open(pjoin(folder, "proteins", "faas", k + ".faa"), "w") as handle :
                handle.writelines(aas)

    full_of_Ns = "rm nucleotides/fnas/GCF_000291845.1.fna nucleotides/fnas/GCF_000291925.1.fna"

    def compute_fastamd5(file):
        with open(file) as handle:
            lines = "".join([l.strip() if not l.startswith(">") else ">" for l in handle ]).split(">")
        lines = sorted(lines)
        md5s = [md5(l.upper().encode('utf-8')).hexdigest() for l in lines]
        full_md5 = md5("".join(md5s).encode("utf-8")).hexdigest()

        return {'full_md5' : full_md5, "entry_md5s" : md5s}

    md5s_entrez = {f[:-4]: compute_fastamd5('nucleotides/fnas/' + f) for f in tqdm(os.listdir('nucleotides/fnas/'))}
    md5s_gorg = {f[:-6]: compute_fastamd5('nucleotides/gorgs/' + f) for f in tqdm(os.listdir('nucleotides/gorgs/'))}
    common_md5s = {md5['full_md5'] for md5 in md5s_entrez.values()}.intersection({md5['full_md5'] for md5 in md5s_gorg.values()})

    good_gorgs = [k for k,v in md5s_gorg.items() if v['full_md5'] not in common_md5s]
    for f in good_gorgs:
        shutil.copy('nucleotides/gorgs/' + f + ".fasta", 'nucleotides/fnas/GORG_' + f.replace("_contigs",".1.fna"))

    command = 'prokka --force --cpus 20 --genus Prochlorococcus --genus Prochlorococcus --locustag {id} --prefix {id} --outdir nucleotides/prokkas/  nucleotides/fnas/{id}.fna > /dev/null '

    for f in os.listdir('nucleotides/fnas/'):
        if f.startswith("GORG_"):
            rs_id = f[:-4]
            print(rs_id)
            call(command.format(id = rs_id), shell = True)

    "checkm taxonomy_wf -x fna -t20 genus Prochlorococcus nucleotides/fnas/  checkm   > static_data/checkm.txt"



    os.rmtree("static_data")

faas = { g[:-4] : "nucleotides/prokkas/" + g for g in  os.listdir("nucleotides/prokkas/") if g.endswith(".faa")}
checkm = parse_checkm("static_data/checkm.txt")

with open("static_data/Prochlos_mOTUs.json") as handle:
    mOTUs = json.load(handle)
mOTUS = {k : [vv['name'] for vv in v['genomes']] for k, v in mOTUs.items()}
genome2motu = {g : k for k,v in mOTUS.items() for g in v}

with open("static_data/gtdbtk/gtdbtk.bac120.summary.tsv") as handle:
    taxonomy = {l.split("\t")[0][1:] : l.split("\t")[1] for l in handle}
consensus = {k : [taxonomy.get(vv) for vv in v] for k,v in mOTUS.items()}
consensus = {k : [(s, v.count(s)) for s in set(v)] for k,v in consensus.items()}
consensus = {k: max(v, key = lambda x : x[1])[0] for k,v in consensus.items()}
blacklist = ['GCF_000634395.1']


for v in blacklist:
    del checkm[v]

full_genomes = [g for g in faas if g in checkm and checkm[g]['Completeness'] > 95]
good_genomes = [g for g in faas if g in checkm and checkm[g]['Completeness'] > 70]
decent_genomes = [g for g in faas if g in checkm and checkm[g]['Completeness'] > 40]
useable_genomes = [g for g in faas if g in checkm and checkm[g]['Completeness'] > 40]

p_As = [g for g in mOTUS['Prochlos_mOTU_002'] if g not in blacklist]
full_p_As = [k for k in p_As if k in full_genomes]

with open("nucleotides/gff_list", "w") as handle:
    handle.writelines([c + "\tnucleotides/prokkas/" + c + ".gff\n" for c in p_As])

with open("nucleotides/gff_list_fulls", "w") as handle:
    handle.writelines([c + "\tnucleotides/prokkas/" + c + ".gff\n" for c in full_p_As])


def running_tools():
    run("""
    cat nucleotides/prokkas/*.faa > nucleotides/all_faas.faa
    emapper.py --cpu 20 -m diamond --output nucleotides/all_faas -d bact -i nucleotides/all_faas.faa

#    ppanggolin annotate --anno nucleotides/gff_list  -o static_data/ --basename ppanggolin_clusters --use_pseudo --tmpdir . -c 20 -f
#    ppanggolin cluster -p static_data/ppanggolin_clusters.h5 -c 20 --tmpdir .
#    ppanggolin graph -p static_data/ppanggolin_clusters.h5 -c 20
#    ppanggolin partition -K3 -f  -p static_data/ppanggolin_clusters.h5 --cpu 20
#    ppanggolin rarefaction -c 20 --tmpdir . -K3 -f  --min 15 --depth 1 --max 1054 -o other_soft/ppanggolin/rarefaction -p static_data/ppanggolin_clusters.h5

#    mOTUconvert.py --in_type ppanggolin static_data/ppanggolin_clusters.h5 > static_data/ppanggolin.gid2cog

    ppanggolin annotate --anno nucleotides/gff_list  -o static_data/ --basename ppanggolin_clusters_species --use_pseudo --tmpdir . -c 20 -f
    ppanggolin cluster -p static_data/ppanggolin_clusters_species.h5 -c 20 --tmpdir .
    ppanggolin graph -p static_data/ppanggolin_clusters_species.h5 -c 20
    ppanggolin partition -K3 -f  -p static_data/ppanggolin_clusters_species.h5 --cpu 20
    ppanggolin rarefaction -c 20 --tmpdir . -K3 -f  --min 15 --depth 20 --max 10000 -o other_soft/ppanggolin/rarefaction_species -p static_data/ppanggolin_clusters_species.h5

    python ../../mOTUlizer/bin/mOTUconvert.py --in_type ppanggolin static_data/ppanggolin_clusters_species.h5 > static_data/ppanggolin_species.gid2cog

    ppanggolin annotate --anno nucleotides/gff_list_fulls  -o static_data/ --basename ppanggolin_clusters_goods --use_pseudo --tmpdir . -c 20 -f
    ppanggolin cluster -p static_data/ppanggolin_clusters_goods.h5 -c 20 --tmpdir .
    ppanggolin graph -p static_data/ppanggolin_clusters_goods.h5 -c 20
    ppanggolin partition -K3 -f  -p static_data/ppanggolin_clusters_goods.h5 --cpu 20

    python ../../mOTUlizer/bin/mOTUconvert.py --in_type ppanggolin static_data/ppanggolin_clusters_goods.h5 > static_data/ppanggolin_goods.gid2cog

    mkdir other_soft/roary/s__Prochlorococcus_A
    roary -p 22 -o other_soft/roary/s__Prochlorococcus_A_clusters.txt  -cd 1 -v `cat nucleotides/gff_list  | cut -f2`
    cp accessory* blast_identity_frequency.Rtab  core_accessory* gene_presence_absence.* number_of_* summary_statistics.txt other_soft/roary/s__Prochlorococcus_A
    python ../../mOTUlizer/bin/mOTUconvert.py --in_type roary other_soft/roary/s__Prochlorococcus_A_clusters.txt > static_data/roary_s__Prochlorococcus_A.gid2cog


    mkdir other_soft/roary/goods
    roary -p 20 -o other_soft/roary/goods_clusters.txt  -cd 1 -v `cat nucleotides/gff_list_fulls  | cut -f2`
    mv accessory* blast_identity_frequency.Rtab  core_accessory* gene_presence_absence.* number_of_* summary_statistics.txt other_soft/roary/goods
    mOTUconvert.py --in_type roary other_soft/roary/goods_clusters.txt > static_data/roary_goods.gid2cog

    """, shell=True)

def get_genome_stats(g):
    with open("nucleotides/prokkas/" + g + ".fna") as handle:
        length = 0
        nb_contig = 0
        for l in handle:
            if l.startswith(">") :
                nb_contig +=1
            else :
                length += len(l.strip())
    return {'genome_len' : length, 'nb_contigs' : nb_contig}
genome2stats = {g: get_genome_stats(g) for g in tqdm(genome2motu)}

def run_motupan(genomes, gid2cog, name = "test" , k=15):
    kks = choices(genomes, k = k)
    k_dups = {k : 0 for k in set(kks) if kks.count(k) >1}
    for i,g in enumerate(kks):
        if g in k_dups:
            if k_dups[g] != 0:
                kks[i] = g +"#"+str(k_dups[g])
            k_dups[g] +=1
    faas_loc = {g : faas.get(g,faas[g.split("#")[0]])  for g in kks}
    cog_dict = {k : gid2cog.get(k,gid2cog[k.split("#")[0]])  for  k in kks}
    checkm_loc = {g : checkm.get(g,checkm[g.split("#")[0]])['Completeness']  for g in kks}
    motu = None
    motu = mOTU( name =  name , faas = faas_loc, cog_dict = cog_dict, checkm_dict = checkm_loc, max_it = 20, threads = 20, precluster = True, method = "default")
    stats = motu.get_stats()
#    roc = motu.roc_values()
    eff_genomes = sum([v for v in checkm_loc.values()])/100
    new_eff = sum([g.new_completness for g in motu])/100
    out = { 'nb_genomes' : k, 'core_len' : len(stats['test']['core']), 'aux_len' : len(stats['test']['aux_genome']), "checkm_eff" : eff_genomes, "motu_eff" : new_eff, 'mean_new_complete' : mean([g.new_completness for g in motu]) }
#    out.update(roc)
    return out

def ppanggolin_vs_motulizer():

    with open("static_data/ppanggolin_species.gid2cog") as handle:
        gid2cog_species = json.load(handle)
    gid2cog_species = {k : set(v) for k,v in gid2cog_species.items()}
    j = list(range(15,len(gid2cog_species), 1))*5
    def motupan_species_test_w_ppanggolin_cogs(i):
        return run_motupan(p_As,k=j[i], gid2cog = gid2cog_species)

    motupan_ppanggolin_species = mOTU( name =  "motupan_ppanggolin_species" ,
            faas = {g : v for g,v in faas.items() if g in gid2cog_species} ,
            cog_dict = gid2cog_species,
            checkm_dict = {g : checkm[g]['Completeness'] for g in gid2cog_species}, max_it = 20, method = "default")


    pool = Pool(processes=20)
    motupan_species_rarefact_w_ppanggolin_cogs = pool.map(motupan_species_test_w_ppanggolin_cogs, list(range(len(j))))
    with open("analyses/motupan_species_rarefact_w_ppanggolin_cogs.tsv", "w") as handle:
        head = list(motupan_species_rarefact_w_ppanggolin_cogs[0].keys())
        ppangg2genome2genehandle.writelines(["\t".join(head) + "\n"])
        handle.writelines(["\t".join([str(l[k]) for k in head]) + "\n"  for l in motupan_species_rarefact_w_ppanggolin_cogs])

    handle = h5py.File("static_data/ppanggolin_clusters_species.h5", "r")
    gene2ppangg = {a.decode() : b.decode() for a,b in  tqdm(handle['geneFamilies'])}
    ppang2gene = {s : [] for s in set(gene2ppangg.values())}
    for a,b in gene2ppangg.items():
        ppang2gene[b] += [a]

    ppangg2genome2gene = {k : {genome : [zz for zz in v if zz.startswith(genome)] for genome in set([vv.split("_CDS")[0] for vv in v])}  for k,v in ppang2gene.items()}

    ppanggolin_sets = {a[0].decode() : {'ppanggolin' : 'persistent' if a[1].decode() == "P" else 'shell' if a[1].decode() == "S" else 'cloud'}  for a in tqdm(handle['geneFamiliesInfo'])}
    for k, v in ppanggolin_sets.items():
        ppanggolin_sets[k]['motupan'] = 'core' if k in motupan_ppanggolin_species.core else 'accessory'


    ppanggolin_matrix = {gg : dict({k : 0 for k in gid2cog_species}) for g in gid2cog_species.values() for gg in g}
    for k, gs in gid2cog_species.items():
        for gg in gs:
            ppanggolin_matrix[gg][k] = 1

    pandas.DataFrame.from_dict(ppanggolin_sets, orient = "index").to_csv("analyses/ppanggolin_matrix_species_cogs.csv", index_label = "COG")
    pandas.DataFrame.from_dict(ppanggolin_matrix, orient = "index").to_csv("analyses/ppanggolin_matrix_species.csv", index_label = "cog")
    pandas.DataFrame.from_dict({g.name : {'motupan_completeness' : g.new_completness, 'checkm_completeness' : g.checkm_complet, **genome2stats[g.name]} for g in motupan_ppanggolin_species}, orient="index").to_csv("analyses/ppanggolin_matrix_species_genomes.csv", index_label = "genome")



def pange_dict2roary_classes(gid2cog, mean_complete = 100):
    core_cutoff = 0.99*mean_complete/100
    softcore_cutoff = 0.95*mean_complete/100
    shell_cutoff = 0.15*mean_complete/100
    cloud_cutoff = 0

    nb_genomes = len(gid2cog)
    all_cog_counts = {vv : 0 for v in gid2cog.values() for vv in v}

    for v in gid2cog.values():
        for vv in v:
            all_cog_counts[vv] += 1

    core = {k for k,v in all_cog_counts.items() if v/nb_genomes > core_cutoff}
    for k in core:
        del all_cog_counts[k]

    softcore = {k for k,v in all_cog_counts.items() if v/nb_genomes > softcore_cutoff}
    for k in softcore:
        del all_cog_counts[k]

    shell = {k for k,v in all_cog_counts.items() if v/nb_genomes > shell_cutoff}
    for k in shell:
        del all_cog_counts[k]

    cloud = set(all_cog_counts.keys())
    return {'core' : len(core), 'softcore' : len(softcore), 'shell' : len(shell), 'cloud' : len(cloud) }

def roary_vs_motupan():
    with open("static_data/roary_goods.gid2cog") as handle:
        roary_gi2cog = {k : set(v) for k, v in json.load(handle).items()}
        reps=20
        boots = []
        for i in tqdm(range(3, len(roary_gi2cog))):
            for j in range(reps):
                sub_gid2cog = {k : roary_gi2cog[k] for k in choices(list(roary_gi2cog), k=i)}
                est_complete = mean([checkm[k]['Completeness'] for k in sub_gid2cog])
                tt = pange_dict2roary_classes(sub_gid2cog)
                tt['nb_org'] = i
                tt['rep'] = j
                tt['method'] = "strict"
                tt['est_checkm_complete'] = est_complete
                motupan = run_motupan(list(sub_gid2cog.keys()), k=i, gid2cog = sub_gid2cog)

                tt['motupan_est_checkm'] = motupan['mean_new_complete']
                tt['motupan_core'] = motupan['core_len']
                tt['motupan_cloud'] = motupan['aux_len']

                boots += [tt]

        head = list(boots[0].keys())
        with open("analyses/roary_rarefaction.csv", "w") as handle:
            handle.writelines([",".join(head) + "\n"] + [ ",".join([str(dd[k]) for k in head]) + '\n' for dd in boots])

    with open("static_data/ppanggolin_goods.gid2cog") as handle:
        ppanggolin_gi2cog = {k : set(v) for k, v in json.load(handle).items()}
        boots = []
        for i in tqdm(range(3, len(ppanggolin_gi2cog))):
            for j in range(reps):
                sub_gid2cog = {k : ppanggolin_gi2cog[k] for k in choices(list(ppanggolin_gi2cog), k=i)}
                est_complete = mean([checkm[k]['Completeness'] for k in sub_gid2cog])
                tt = pange_dict2roary_classes(sub_gid2cog)
                tt['nb_org'] = i
                tt['rep'] = j
                tt['method'] = "strict"
                tt['est_checkm_complete'] = est_complete
                motupan = run_motupan(list(sub_gid2cog.keys()),k=i, gid2cog = sub_gid2cog)

                tt['motupan_est_checkm'] = motupan['mean_new_complete']
                tt['motupan_core'] = motupan['core_len']
                tt['motupan_cloud'] = motupan['aux_len']

                boots += [tt]

        head = list(boots[0].keys())
        with open("analyses/ppanggolin_rarefaction.csv", "w") as handle:
            handle.writelines([",".join(head) + "\n"] + [ ",".join([str(dd[k]) for k in head]) + '\n' for dd in boots])



def processing_goods():
    with open("static_data/roary_goods.gid2cog") as handle:
        roary = {k : set(v) for k, v in json.load(handle).items()}
    with open("static_data/ppanggolin_goods.gid2cog") as handle:
        ppanggolin = {k : set(v) for k, v in json.load(handle).items()}

    motupan_roary = mOTU( name =  "motupan_roary" ,
            faas = {g : v for g,v in faas.items() if g in roary} ,
            cog_dict = {g : v for g,v in roary.items() if g in roary},
            checkm_dict = {g : checkm[g]['Completeness'] for g in roary}, max_it = 20, method = "default")

    motupan_ppanggolin = mOTU( name =  "motupan_ppanggolin" ,
            faas = {g : v for g,v in faas.items() if g in ppanggolin} ,
            cog_dict = {g : v for g,v in ppanggolin.items() if g in ppanggolin},
            checkm_dict = {g : checkm[g]['Completeness'] for g in ppanggolin}, max_it = 20, method = "default")

    motupan_roary.roc_values()
    motupan_ppanggolin.roc_values()
    handle = h5py.File("static_data/ppanggolin_clusters_goods.h5", "r")
    gene2ppangg = {a.decode() : b.decode() for a,b in  tqdm(handle['geneFamilies'])}
    ppang2gene = {s : [] for s in set(gene2ppangg.values())}
    for a,b in gene2ppangg.items():
        ppang2gene[b] += [a]

    ppangg2genome2gene = {k : {genome : [zz for zz in v if zz.startswith(genome)] for genome in set([vv.split("_CDS")[0] for vv in v])}  for k,v in ppang2gene.items()}

    ppanggolin_sets = {a[0].decode() : {'ppanggolin' : 'persistent' if a[1].decode() == "P" else 'shell' if a[1].decode() == "S" else 'cloud'}  for a in tqdm(handle['geneFamiliesInfo'])}
    for k, v in ppanggolin_sets.items():
        ppanggolin_sets[k]['motupan'] = 'core' if k in motupan_ppanggolin.core else 'accessory'



    ppanggolin_matrix = {gg : dict({k : 0 for k in ppanggolin}) for g in ppanggolin.values() for gg in g}
    for k, gs in ppanggolin.items():
        for gg in gs:
            ppanggolin_matrix[gg][k] = 1

    pandas.DataFrame.from_dict(ppanggolin_matrix, orient = "index").to_csv("analyses/ppanggolin_matrix_species.csv", index_label = "cog")
    pandas.DataFrame.from_dict({g.name : {'motupan_completeness' : g.new_completness, 'checkm_completeness' : g.checkm_complet, **genome2stats[g.name]} for g in motupan_ppanggolin}, orient="index").to_csv("analyses/ppanggolin_matrix_genomes.csv", index_label = "genome")

    roary_matrix = {gg : dict({k : 0 for k in roary}) for g in roary.values() for gg in g}
