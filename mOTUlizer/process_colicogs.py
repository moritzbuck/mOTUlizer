
def main():
    import os
    import pandas
    import shutil
    import sys
    from os.path import join as pjoin
    from random import choice, sample


    sys.path += ["/home/moritz/repos/moritz/0039_mOTUlizer"]

    from mOTUlizer.classes.mOTU import mOTU
    from tqdm import tqdm

    print("Loading stats")
    stats = pandas.read_csv("/home/moritz/temp/cores/ecolis_magstats.csv", index_col = 0)

    print("Loading taxonomy")
    taxonomy = {l : "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli".split(";") for l in stats.index}

    print("Loading cogs")
    with open("/home/moritz/temp/cores/ecolis_mags2cogs.csv") as handle:
        handle.readline()
        mag2cog = {l.split()[0] : l[:-1].split()[1].split(";") for l in handle if len(l[:-1].split()[1:]) > 0 }

    print("Loading ANIs")
    ani_dict = {}

    data_pack = {'mag2cog' : mag2cog, 'taxonomy' : taxonomy, 'stats' : stats, 'base_folder' : "/home/moritz/repos/moritz/0039_mOTUlizer/test_data/"}

    otu = mOTU( name = "EColi", members = [m for m in mag2cog.keys() if stats.loc[m,"completeness"] > 70], data_pack = data_pack, precomp_ani = ani_dict, funct_derep = 0.95)


    out_dir = "outputs/" + "colicogs"

    os.makedirs(out_dir, exist_ok = True)

    otu_stats = otu.get_otu_stats()
    pandas.DataFrame.from_records(otu_stats, index="otu").to_csv(pjoin(out_dir, "otu.stats"))

    count_stats = otu.get_count_stats()
    pandas.DataFrame.from_records(count_stats).to_csv(pjoin(out_dir,"counts.stats"))

    rarefy_stats = list(otu.rarefy_pangenome().values())
    rarefy_stats = pandas.DataFrame.from_records(rarefy_stats)
    single_rarefy_stats = list(otu.rarefy_pangenome(singletons = True).values())
    single_rarefy_stats = pandas.DataFrame.from_records(single_rarefy_stats)
    rarefy_stats['singletons'] = False
    single_rarefy_stats['singletons'] = True
    rarefy_stats = rarefy_stats.append(single_rarefy_stats)
    rarefy_stats.to_csv(pjoin(out_dir,"rarefy.stats"))

if __name__ == "__main__":
    main()
