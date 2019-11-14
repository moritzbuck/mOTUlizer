
#    stats = pandas.read_csv("test_data/magstats.csv", index_col = 0)
#
#    with open("test_data/taxonomy.tax") as handle:
#        taxonomy = {l.split(",")[0] : l[:-1].split(",")[1:] for l in handle}
#
#    with open("test_data/mag2cog.tsv") as handle:
#        mag2cog = {l.split()[0] : l[:-1].split()[1:] for l in handle}
#
#    with open("test_data/fastani.tsv") as handle:
#        ani_dict = {( l.split()[0], l.split()[1] ) : float(l.split()[2]) for l in handle}


def main():
    import os
    import pandas
    import shutil
    import sys
    from os.path import join as pjoin
    import json

    sys.path += ["/home/moritz/repos/moritz/0039_mOTUlizer"]
    sys.path += ["/home/moritz/repos/moritz/0042_emapper2json"]

    from mOTUlizer.classes.mOTU import mOTU
    from mOTUlizer.config import *

    from tqdm import tqdm
    from emapper2json import _COG_CATS_, get_blocks

    input_file = sys.argv[1]
    input_data_folder = sys.argv[2]
# /home/moritz/temp/cores

    print("Loading stats")
    stats = pandas.read_csv(input_data_folder  + "/magstats.csv", index_col = 0)

    print("Loading taxonomy")
    with open(input_data_folder , "/full_taxonomy.tax") as handle:
        taxonomy = {l.split(",")[0] : l[:-1].split(",")[1:] for l in handle}

    print("Loading cogs")
    with open(input_data_folder  + "/mags/mag2cogs.tsv") as handle:
        mag2cog = {l.split()[0] : l[:-1].split()[1:] for l in handle}

    print("Loading ANIs")
    with open(input_data_folder  + "/mags/fastani_pairs.csv") as handle:
        handle.readline()
        ani_dict = {( l.split()[0], l.split()[1] ) : float(l.split()[2]) for l in handle}

    print("Loading annotations")
    with open(input_data_folder  + "/proteomics/cog2annot.json") as handle:
        cog2annot = json.load(handle)

    print("Loading kegg_paths")
    with open( DB_FOLDER + "/ec2kegg/kegg_maps.txt") as handle:
        kegg2name = {}
        lines = handle.readlines()
        for l in lines:
            if l.startswith("    "):
                kegg2name["map" + l.split()[0]] = " ".join(l.split()[1:])
                kegg2name["ko" + l.split()[0]] = " ".join(l.split()[1:])

    print("Loading gonames")
    with open(DB_FOLDER + "/go.obo") as handle:
        go2name = {}
        for l in handle:
            if l.startswith("id: "):
                go = l.split()[1]
            if l.startswith("name: "):
                go2name[go] = " ".join(l.split()[1:])

    print("Loading kegg_paths")
    module2name = {}
    module2def = {}
    module2ko = {}
    for m in os.listdir(DB_FOLDER + "/ec2kegg/modules"):
        with open(DB_FOLDER + "/ec2kegg/modules/" + m) as handle:
            lines = handle.readlines()
            for l in lines:
                if l.startswith("DEFINITION"):
                    module2def[m] = " ".join(l.split()[1:])
                if l.startswith("NAME"):
                    module2name[m] = " ".join(l.split()[1:])

    data_pack = {'mag2cog' : mag2cog, 'taxonomy' : taxonomy, 'stats' : stats, 'base_folder' : "/home/moritz/people/0023_anoxicencyclo/4500_assembly_analysis/mags/"}


    otu_list = []
    with open(input_file) as handle:
        for i,l in tqdm(enumerate(handle)):
            if i < 10:
                name = l.split()[0]
                bins = l.split()[1].split(";")
#                print("processing", name )
                if len(bins) > 3 :
                    otu_list += [ mOTU( name = name, members = bins, data_pack = data_pack, precomp_ani = ani_dict, funct_derep = 0.99)]

    otu2mods = {(otu.name, k) : v for otu in tqdm(otu_list) for k,v in otu.get_kos(cog2annot, module2def).items()}

    out_dir = "outputs/" + input_file.split("ani_")[-1][:-4]
    os.makedirs(out_dir, exist_ok=True)

    threashold = 0.3
    otu2def = {}
    for otu in tqdm(otu_list) :
        cogs = otu.core
        kos = [max(cog2annot[c]['KEGG_ko'].items(), key = lambda x: x[1]) for c in cogs if c in cog2annot and cog2annot[c]['KEGG_ko']]
        core_kos = {k.replace("ko:","") for k,v in kos if v > threashold}

        cogs = [c for c,v in otu.cogCounts.items() if c not in otu.core]
        kos = [max(cog2annot[c]['KEGG_ko'].items(), key = lambda x: x[1]) for c in cogs if c in cog2annot and cog2annot[c]['KEGG_ko']]
        pangenome_kos = {k.replace("ko:","")  for k,v in kos if v > threashold}
        for k, v in module2def.items():
            if otu2mods[(otu.name, "in_core")][k] > 0:
                tt = v
                for ko in core_kos:
                    tt = tt.replace(ko, "**" + ko + "**")
                otu2def[(otu.name, k, "in_core")] = tt
            if otu2mods[(otu.name, "in_aux")][k] > 0:
                tt = v
                for ko in pangenome_kos:
                    tt = tt.replace(ko, "**" + ko + "**")
                otu2def[(otu.name, k, "in_aux")] = tt
    otu2def = sorted(otu2def.items(), key = lambda x: (x[0][1], x[0][0], x[0][2]))

    pandas.DataFrame.from_dict(otu2mods).transpose().to_csv(pjoin(out_dir, "mod_matrix.csv"))
    with open(pjoin(out_dir, "mod_defs.txt"), "w") as handle:
        handle.writelines([";".join(k) + ";" + v  + "\n" for k,v in otu2def])

    pandas.DataFrame.from_dict({otu.name + "_" + k : v  for otu in otu_list for k, v in otu.annot_partition(cog2annot=cog2annot, annotation = "GOs")['counts'].items() }).transpose().fillna(0).to_csv("COG_cats.csv")


    os.makedirs(out_dir, exist_ok = True)

    otu_stats = [otu.get_otu_stats() for otu in tqdm(otu_list)]
    pandas.DataFrame.from_records(otu_stats, index="otu").to_csv(pjoin(out_dir, "otu.stats"))

    count_stats = sum([otu.get_count_stats() for otu in tqdm(otu_list)],[])
    pandas.DataFrame.from_records(count_stats).to_csv(pjoin(out_dir,"counts.stats"))

    rarefy_stats = sum([ list(otu.rarefy_pangenome().values()) for otu in tqdm(otu_list)],[])
    rarefy_stats = pandas.DataFrame.from_records(rarefy_stats)
    single_rarefy_stats = sum([ list(otu.rarefy_pangenome(singletons = True).values()) for otu in tqdm(otu_list)],[])
    single_rarefy_stats = pandas.DataFrame.from_records(single_rarefy_stats)
    rarefy_stats['singletons'] = False
    single_rarefy_stats['singletons'] = True
    rarefy_stats = rarefy_stats.append(single_rarefy_stats)
    rarefy_stats.to_csv(pjoin(out_dir,"rarefy.stats"))




if __name__ == "__main__":
    main()
