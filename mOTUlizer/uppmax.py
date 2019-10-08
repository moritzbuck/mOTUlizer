
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

    sys.path += ["/home/moritz/repos/moritz/0039_mOTUlizer"]

    from mOTUlizer.classes.mOTU import mOTU
    from tqdm import tqdm


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

    data_pack = {'mag2cog' : mag2cog, 'taxonomy' : taxonomy, 'stats' : stats, 'base_folder' : "/home/moritz/people/0023_anoxicencyclo/4500_assembly_analysis/mags/"}

    otu_list = []
    with open(input_file) as handle:
        for i,l in enumerate(handle):
    #        if i <3 :
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

if __name__ == "__main__":
    main()
