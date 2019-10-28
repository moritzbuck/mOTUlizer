


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


    input_file = "test_data/mOTUs.txt"

    stats = pandas.read_csv("test_data/magstats.csv", index_col = 0)

    with open("test_data/taxonomy.tax") as handle:
        taxonomy = {l.split(",")[0] : l[:-1].split(",")[1:] for l in handle}

    with open("test_data/mag2cog.tsv") as handle:
        mag2cog = {l.split()[0] : l[:-1].split()[1:] for l in handle}

    with open("test_data/fastani.tsv") as handle:
        ani_dict = {( l.split()[0], l.split()[1] ) : float(l.split()[2]) for l in handle}

    print("Loading annotations")
    with open("/home/moritz/temp/cores/cog2annot.json") as handle:
        cog2annot = json.load(handle)


    data_pack = {'mag2cog' : mag2cog, 'taxonomy' : taxonomy, 'stats' : stats, 'base_folder' : "/home/moritz/repos/moritz/0039_mOTUlizer/test_data/"}

    otu_list = []
    with open(input_file) as handle:
        for i,l in enumerate(handle):
                name = l.split()[0]
                bins = l.split()[1].split(";")
                print("processing", name )
                if len(bins) > 3 :
                    otu_list += [ mOTU( name = name, members = bins, data_pack = data_pack, precomp_ani = ani_dict)]


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



if __name__ == "__main__":
    main()
