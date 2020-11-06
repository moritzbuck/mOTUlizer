import tempfile
import sys, os
import gzip
from mOTUlizer.config import *
from Bio import SeqIO

def compute_COGs(faas, name = "silixCOGs", precluster = False, threads = 4):
    temp_file = tempfile.NamedTemporaryFile().name
    temp_out = tempfile.NamedTemporaryFile().name
    temp_clust = tempfile.NamedTemporaryFile().name


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

    print("concatenating all faas", file = sys.stderr)

    if all([faa.endswith(".gz") for faa in faas.values()]):
        cat = "zcat "
    elif all([not faa.endswith(".gz") for faa in faas.values()]) :
        cat = "cat "
    else :
        sys.exit("Please make my life easy, either all faas gziped or none gzipped ...")

    os.system(cat + " ".join([*faas.values()]) + " > " + temp_file)

    if precluster:
        cdhit_file = tempfile.NamedTemporaryFile().name

        exec = "cd-hit -i {input} -o {output} -c 0.95 -M 0 -T {threads} -d 0 -s 0.95 >> {log} 2>&1".format(input = temp_file, output = cdhit_file, threads = threads, log = "/dev/null")
        print("Running cd-hit preclustering", file = sys.stderr)
        os.system(exec)

        print("parsing cd-hit", file = sys.stderr)

        with open(cdhit_file + ".clstr") as handle:
            clusters = "\n".join(handle.readlines()).split("Cluster ")

        os.remove(cdhit_file + ".clstr")
        os.remove(cdhit_file)

        clusters = [c.split("\n\n") for c in clusters[1:] if "*" in c]
        clusters = [[cc.split(">")[1].split("... ") for cc in c if ">" in cc and cc != ">"] for c in clusters ]
        clusters = {[cc[0] for cc in c if cc[1] == "*" or cc[1] == "*\n"][0] : [cc[0] for cc in c] for c in clusters}

        print("For", len(prot2faa), "CDSes we got ", len(clusters), "preclusters")
        seqs = [s for s in SeqIO.parse(temp_file, "fasta") if s.id in clusters]
        SeqIO.write(seqs, temp_file, "fasta")

    print("all v all diamond for silix", file = sys.stderr)

    os.system("diamond makedb --db {faas} --in {faas} > /dev/null 2> /dev/null".format(faas = temp_file))
    os.system("diamond blastp --more-sensitive -p {threads} -f 6 -q {faas} --db {faas} -o {out} 2> /dev/null > /dev/null".format(faas = temp_file, out = temp_out, threads = threads))

    print("running silix", file = sys.stderr)
    os.system("silix {faas} {out} > {clust_temp} #2> /dev/null".format(faas = temp_file, out = temp_out, clust_temp = temp_clust))

    print("parsing silix", file = sys.stderr)
    with open(temp_clust) as handle:
        if precluster:
            recs = {g : l[:-1].split()[0]  for l in handle for g in clusters[l[:-1].split()[1]]}
        else :
            recs = {l[:-1].split()[1] : l[:-1].split()[0]  for l in handle}

    #pretty formating names
    fill = max([len(v) for v in recs.values()])
    recs = {k : name + "_" + v.zfill(fill) for k, v in recs.items()}
    genome2cog = {k : set() for k in faas.keys()}
    for k,v in recs.items():
        for vv in prot2faa[k]:
            genome2cog[vv].update([v])

    os.remove(temp_file)
    os.remove(temp_out)
    os.remove(temp_clust)

    return { 'genome2cogs' : genome2cog, 'aa2cog' : recs}
