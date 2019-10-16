
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
    from tqdm import tqdm
    from emapper2json import _COG_CATS_

    input_file = sys.argv[1]

    print("Loading stats")
    stats = pandas.read_csv("/home/moritz/people/0023_anoxicencyclo/4500_assembly_analysis/magstats.csv", index_col = 0)

    print("Loading taxonomy")
    with open("/home/moritz/people/0023_anoxicencyclo/4500_assembly_analysis/full_taxonomy.tax") as handle:
        taxonomy = {l.split(",")[0] : l[:-1].split(",")[1:] for l in handle}

    print("Loading cogs")
    with open("/home/moritz/people/0023_anoxicencyclo/4500_assembly_analysis/mags/mag2cogs.tsv") as handle:
        mag2cog = {l.split()[0] : l[:-1].split()[1:] for l in handle}

    print("Loading ANIs")
    with open("/home/moritz/people/0023_anoxicencyclo/4500_assembly_analysis/mags/fastani_pairs.csv") as handle:
        handle.readline()
        ani_dict = {( l.split()[0], l.split()[1] ) : float(l.split()[2]) for l in handle}

    print("Loading annotations")
    with open("/home/moritz/temp/cores/cog2annot.json") as handle:
        cog2annot = json.load(handle)

    print("Loading kegg_paths")
    with open("/home/moritz/dbs/ec2kegg/kegg_maps.txt") as handle:
        kegg2name = {}
        lines = handle.readlines()
        for l in lines:
            if l.startswith("    "):
                kegg2name["map" + l.split()[0]] = " ".join(l.split()[1:])
                kegg2name["ko" + l.split()[0]] = " ".join(l.split()[1:])

    print("Loading gonames")
    with open("/home/moritz/dbs/go.obo") as handle:
        go2name = {}
        for l in handle:
            if l.startswith("id: "):
                go = l.split()[1]
            if l.startswith("name: "):
                go2name[go] = " ".join(l.split()[1:])


    data_pack = {'mag2cog' : mag2cog, 'taxonomy' : taxonomy, 'stats' : stats, 'base_folder' : "/home/moritz/people/0023_anoxicencyclo/4500_assembly_analysis/mags/"}


    otu_list = []
    with open(input_file) as handle:
        for i,l in tqdm(enumerate(handle)):
            if i < 10:
                name = l.split()[0]
                bins = l.split()[1].split(";")
#                print("processing", name )
                if len(bins) > 5 :
                    otu_list += [ mOTU( name = name, members = bins, data_pack = data_pack, precomp_ani = ani_dict, funct_derep = 0.95)]

    pvs = {(otu.name, k.replace("ko","map")) : v  for otu in otu_list for k,v in otu.annot_partition(cog2annot=cog2annot, annotation = "KEGG_Pathway")['pvals'].items()}
    pvs = {k : v for k, v in zip(pvs.keys(),multipletests(list(pvs.values()), method = 'fdr_bh')[1]) if v < 0.001}
    counts = {v : 0 for k, v in pvs.keys()}
    for k,v in pvs.items():
        counts[k[1]] += 1
    kegg_hits = sorted({kegg2name.get(k, k) :  v for k,v in counts.items() if v > 10}.items(), key = lambda x: x[1])

    pvs = {(otu.name, k.replace("ko","map")) : 1-v  for otu in otu_list for k,v in otu.annot_partition(cog2annot=cog2annot, annotation = "KEGG_Pathway")['pvals'].items()}
    pvs = {k : v for k, v in zip(pvs.keys(),multipletests(list(pvs.values()), method = 'fdr_bh')[1]) if v > 0.001}
    counts = {v : 0 for k, v in pvs.keys()}
    for k,v in pvs.items():
        counts[k[1]] += 1


    pvs = {(otu.name, k) : v  for otu in otu_list for k,v in otu.annot_partition(cog2annot=cog2annot, annotation = "GOs")['pvals'].items()}
    pvs = {k : v for k, v in zip(pvs.keys(),multipletests(list(pvs.values()), method = 'fdr_bh')[1]) if v < 0.001}
    counts = {v : 0 for k, v in pvs.keys()}
    for k,v in pvs.items():
        counts[k[1]] += 1
    sorted({go2name.get(k,k):  v for k,v in counts.items() if v > 10}.items(), key = lambda x: x[1])[-50:]

    pvs = {(otu.name, k) : v  for otu in otu_list for k,v in otu.annot_partition(cog2annot=cog2annot)['pvals'].items()}
    pvs = {k : v for k, v in zip(pvs.keys(),multipletests(list(pvs.values()), method = 'fdr_bh')[1]) if v < 0.001}
    counts = {v : 0 for k, v in pvs.keys()}
    for k,v in pvs.items():
        counts[k[1]] += 1
    sorted({_COG_CATS_.get(k,k):  v for k,v in counts.items() if v > 10}.items(), key = lambda x: x[1])[-50:]



    out_dir = "outputs/" + input_file.split("ani_")[-1][:-4]

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
