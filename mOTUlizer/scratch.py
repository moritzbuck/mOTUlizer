import os
import pandas
import shutil
import sys
sys.path += ["/home/moritz/repos/moritz/0039_mOTUlizer"]

from mOTUlizer.classes.mOTU import mOTU
from tqdm import tqdm

def main():
    stats = pandas.read_csv("test_data/magstats.csv", index_col = 0)

    with open("test_data/taxonomy.tax") as handle:
        taxonomy = {l.split(",")[0] : l[:-1].split(",")[1:] for l in handle}

    with open("test_data/mag2cog.tsv") as handle:
        mag2cog = {l.split()[0] : l[:-1].split()[1:] for l in handle}

    data_pack = {'mag2cog' : mag2cog, 'taxonomy' : taxonomy, 'stats' : stats, 'base_folder' : "/home/moritz/repos/moritz/0039_mOTUlizer/test_data/"}

    otu_list = []
    with open("test_data/mOTUs.txt") as handle:
        for l in handle:
            name = l.split()[0]
            bins = l.split()[1].split(";")
            otu_list += [ mOTU( name = name, members = bins, data_pack = data_pack) ]
    print({t.name : t.mean_overlap() for t in tqdm(otu_list)})

if __name__ == "__main__":
    main()
