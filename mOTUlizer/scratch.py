
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

    sys.path += ["/home/moritz/repos/moritz/0039_mOTUlizer"]

    from mOTUlizer.classes.mOTU import mOTU
    from tqdm import tqdm

    stats = pandas.read_csv("/home/moritz/temp/cores/magstats.csv", index_col = 0)

    with open("/home/moritz/temp/cores/full_taxonomy.tax") as handle:
        taxonomy = {l.split(",")[0] : l[:-1].split(",")[1:] for l in handle}

    with open("/home/moritz/temp/cores/mag2cogs.tsv") as handle:
        mag2cog = {l.split()[0] : l[:-1].split()[1:] for l in handle}

    with open("/home/moritz/temp/cores/fastani_pairs.csv") as handle:
        ani_dict = {( l.split()[0], l.split()[1] ) : float(l.split()[2]) for l in handle}

    data_pack = {'mag2cog' : mag2cog, 'taxonomy' : taxonomy, 'stats' : stats, 'base_folder' : "/home/moritz/repos/moritz/0039_mOTUlizer/test_data/"}

    otu_list = []
    with open("test_data/mOTUs.txt") as handle:
        for l in handle:
            name = l.split()[0]
            bins = l.split()[1].split(";")
            print("processing", name )
            if len(bins) > 3 :
                otu_list += [ mOTU( name = name, members = bins, data_pack = data_pack, precomp_ani = ani_dict)]


    otu_stats = [otu.get_otu_stats() for otu in tqdm(otu_list)]
    pandas.DataFrame.from_records(otu_stats, index="otu").to_csv("test.stats")

    count_stats = sum([otu.get_count_stats() for otu in tqdm(otu_list)],[])
    pandas.DataFrame.from_records(count_stats).to_csv("counts.stats")

    rarefy_stats = sum([ list(otu.rarefy_pangenome().values()) for otu in tqdm(otu_list)],[])
    pandas.DataFrame.from_records(rarefy_stats).to_csv("rarefy.stats")

    single_rarefy_stats = sum([ list(otu.rarefy_pangenome(singletons = True).values()) for otu in tqdm(otu_list)],[])
    pandas.DataFrame.from_records(single_rarefy_stats).to_csv("single_rarefy.stats")

if __name__ == "__main__":
    main()
