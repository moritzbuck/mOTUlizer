import tempfile
import sys, os
import gzip
from mOTUlizer.config import *

def compute_COGs(faas, name = "silixCOGs", precluster = False):
    temp_file = tempfile.NamedTemporaryFile().name
    temp_out = tempfile.NamedTemporaryFile().name
    temp_clust = tempfile.NamedTemporaryFile().name

    if precluster:
        sys.exit("Preclustering for COGs, not implemented yet!")

    prot_ids = set()
    prot2faa ={}
    for k,v in faas.items():
        if not v.endswith(".gz"):
            with open(v) as handle:
                ids = [l[:-1].split()[0][1:] for l in handle if l[0] == ">"]
        else :
            with gzip.open(v) as handle:
                ids = [l.decode()[:-1].split()[0][1:] for l in handle if l.decode()[0] == ">"]
        for i in ids:
            if i in prot2faa:
                prot2faa[i] += [k]
            else :
                prot2faa[i] = [k]

    print("all v all diamond for silix", file = sys.stderr)

    if all([faa.endswith(".gz") for faa in faas.values()]):
        cat = "zcat "
    elif all([not faa.endswith(".gz") for faa in faas.values()]) :
        cat = "cat "
    else :
        sys.exit("Please make my life easy, either all faas gziped or none gzipped ...")

    os.system(cat + " ".join([*faas.values()]) + " > " + temp_file)
    os.system("diamond makedb --db {faas} --in {faas} > /dev/null 2> /dev/null".format(faas = temp_file))
    os.system("diamond blastp --more-sensitive -p {threads} -f 6 -q {faas} --db {faas} -o {out} 2> /dev/null > /dev/null".format(faas = temp_file, out = temp_out, threads = THREADS))

    print("running silix", file = sys.stderr)
    os.system("silix {faas} {out} > {clust_temp} #2> /dev/null".format(faas = temp_file, out = temp_out, clust_temp = temp_clust))

    print("parsing silix", file = sys.stderr)
    with open(temp_clust) as handle:
        recs = {l[:-1].split()[1] : l[:-1].split()[0]  for l in handle}

    #pretty formating names
    fill = max([len(v) for v in recs.values()])
    recs = {k : name + "_" + v.zfill(fill) for k, v in recs.items()}
    genome2cog = {k : set() for k in faas.keys()}
    for k,v in recs.items():
        for vv in prot2faa[k]:
            genome2cog[vv].update([v])
    return { 'genome2cogs' : genome2cog, 'aa2cog' : recs}
