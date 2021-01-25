import json
from random import sample
from mOTUlizer import __version__
from mOTUlizer.classes import *
from mOTUlizer.utils import *
from mOTUlizer.classes.mOTU import mOTU
from multiprocessing import Pool

init_compl = 95
max_it = 100

with open("g__Escherichia.silix.clusters.gid2cog") as handle:
    gid2cog = json.load(handle)

def run_motupan(i, name = "test" , k=100):
    print("Starting iteration", i)
    kks = sample(list(gid2cog.keys()), k = k)
    sub_gid = {k : set(gid2cog[k]) for k in kks}
    motu = mOTU( name =  name , faas = {} , cog_dict = sub_gid, checkm_dict = {g : init_compl for g in sub_gid}, max_it = max_it)
    stats = motu.get_stats()
    print("Ending iteration", i)
    return { 'nb_genomes' : k, 'core_len' : len(stats['test']['core']), 'aux_len' : len(stats['test']['aux_genome'])}

pool = Pool(processes=16)

p = pool.map(run_motupan, list(range(100)))
all_res = p

def run_k1000(i):
    return run_motupan(i, k=1000)

pool = Pool(processes=16)
p = pool.map(run_k1000, list(range(100)))
all_res += p

def run_k10000(i):
    return run_motupan(i, k=10000)
pool = Pool(processes=16)
p = pool.map(run_k10000, list(range(16)))
all_res += p
