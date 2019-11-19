import tempfile
import sys, os
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
        with open(v) as handle:
            ids = [l[:-1].split()[0] for l in handle if l[0] == ">"]
            prot2faa.update({i[1:] : k for i in ids})
            set_ids = set(ids)
            assert len(ids) == len(set_ids), "the faa {name} has duplicated ids".format(name = k)
            assert len(prot_ids.intersection(set_ids)) == 0, "the faa {name} has duplicated ids with someone else, the id(s) are {bla}[...]".format(name = k, bla = list(prot_ids.intersection(set_ids))[0:10])
            prot_ids.update(set_ids)

    print("all v all diamond for silix", file = sys.stderr)
    os.system("cat " + " ".join([*faas.values()]) + " > " + temp_file)
    os.system("diamond makedb --db {faas} --in {faas} > /dev/null".format(faas = temp_file))
    os.system("diamond blastp --more-sensitive -p {threads} -f 6 -q {faas} --db {faas} -o {out} > /dev/null".format(faas = temp_file, out = temp_out, threads = THREADS))

    print("running silix", file = sys.stderr)
    os.system("silix {faas} {out} > {clust_temp} 2> /dev/null".format(faas = temp_file, out = temp_out, clust_temp = temp_clust))

    print("parsing silix", file = sys.stderr)
    with open(temp_clust) as handle:
        recs = {l[:-1].split()[1] : l[:-1].split()[0]  for l in handle}

    #pretty formating names
    fill = max([len(v) for v in recs.values()])
    recs = {k : name + "_" + v.zfill(fill) for k, v in recs.items()}
    genome2cog = {k : set() for k in faas.keys()}
    for k,v in recs.items():
        genome2cog[prot2faa[k]].update([v])
    return { 'genome2cogs' : genome2cog, 'aa2cog' : recs}
