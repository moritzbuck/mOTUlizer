from tqdm import tqdm
import os
import time
from numpy import mean
import pandas
import multiprocessing
from Bio import SeqIO
from multiprocessing.dummy import Pool as ThreadPool
from itertools import repeat

from mOTUlizer.db.SeqDb import SeqDb
from mOTUlizer.utils import parse_fasta
from timeit import default_timer as timer
from mOTUlizer.classes.mOTU import mOTU
from mOTUlizer.classes.MetaBin import MetaBin
from mOTUlizer.classes.tools.Muscle import Muscle
SeqDb.init_db("/home/moritz/projects/0064_bis/data/all_genomes.sqlite")
db = SeqDb.get_global_db()

print("Number of processors: ", multiprocessing.cpu_count())


folder = "example_files/aquadb_mOTU_00113/"
ffolder = "/home/moritz/temp/pnecs/fnas/"
gfolder = "/home/moritz/temp/pnecs/gffs/"


fnas = { f[:-4] : ffolder + f for f in os.listdir(ffolder) if f.endswith(".fna")}
gffs = { f[:-4] : gfolder + f for f in os.listdir(gfolder) if f.endswith(".gff")}

name, fna = list(fnas.items())[0]
gff = gffs[name]
#for name,fna in tqdm(fnas.items()):
#    contigs = [(s.id, str(s.seq)) for s in SeqIO.parse(fna, "fasta")]
#    db.add_genome(name, contigs = contigs)

cursor = SeqDb.seq_db._connection.cursor()

taxtoadd = [(g[0], json.dumps( { 'gtdbtk_r207' :  gtdb_md.gtdb_taxonomy[g[0]]} ).replace('"', '""') )   for g in data]
for g, t in taxtoadd:
    cursor.execute(f""" UPDATE genomes SET taxonomy = '{t}' WHERE genome_name='{g}' ;""")
conn.commit()

with open("/home/moritz/temp/geneome.list") as handle:
    genomes = [MetaBin(name = l.split("/")[-2], gbk_file = l.strip()) for l in tqdm(list(handle))]

if False:
    print("loading genomes")
    genomes = []
    for g in tqdm(fnas):
        genomes += [MetaBin(g, nucleotide_file = fnas[g], gff_file = gffs[g])]

    print("exporting genomes")
    for g in tqdm(genomes):
        g.export_genome_fasta("test.fna")

    start = timer()
    cust = parse_fasta("test.fna")
    end = timer()

    print("Custom parser time :", end - start)


    start = timer()
    seqio_ = [(s.id, str(s.seq)) for s in SeqIO.parse("test.fna", "fasta")]
    end = timer()

    print("SeqIO parser time :", end - start)
