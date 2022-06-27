#!/usr/bin/env python

import os
import shutil
import sys
from os.path import join as pjoin
import argparse
import json
from random import uniform
from mOTUlizer import __version__
from mOTUlizer.classes import *
from mOTUlizer.utils import *
from mOTUlizer.errors import *
from mOTUlizer.classes.mOTU import mOTU
import multiprocessing
from mOTUlizer.classes.MetaBin import MetaBin

#from mOTUlizer.config import *

description_text = """
From a buch of amino-acid sequences or COG-sets, computes concensus AA/COG sets.

Returns all to stdout by default.
"""
__checkm_default__ = 95

def motupan(args):

    #parse and check your amino-acid files
    if args.txt and args.faas:
        with open(args.faas[0]) as handle:
            faas = {os.path.splitext(os.path.basename(f.strip().rstrip(".gz")))[0] : f.strip() for f in handle.readlines()}
    elif args.faas:
        faas = {os.path.splitext(os.path.basename(f.strip().rstrip(".gz")))[0] : f for f in args.faas}
    else :
        faas = {}

    for f in faas.values():
        if not os.path.exists(f):
            raise FileError(f"one or some of your faas don't exists, more specifically {f}")

    genomes = list(faas)

    out_json = args.output
    checkm = {}
    if args.checkm :
        assert os.path.exists(args.checkm), "The file for checkm does not exists"

        print("Parsing the checkm-file", file = sys.stderr)
        checkm = {k : v['Completeness'] for k,v in parse_checkm(args.checkm).items()}
        checkm = {g : checkm.get(g) for g in genomes}
        if not all([v for v in checkm.values()]):
            print("you did not have completeness for all bins/genomes, abscent values have been put to default value, e.g. " + str(__checkm_default__), file = sys.stderr)
            checkm = {k : v if v  else __checkm_default__ for k,v in checkm.items()}
    elif args.length_seed :
        for f in genomes:
            checkm[f] = None
    elif args.seed :
        for f in genomes:
            checkm[f] = args.seed
    elif args.random_seed :
        for f in genomes:
            checkm[f] = uniform(50,80)
    else :
        for f in genomes:
            checkm[f] = __checkm_default__

    nb_boots = args.boots

    threads = args.threads
    precluster = args.precluster

    name = args.name if args.name else random_name()
    name = (name + "_") if not name.endswith("_") else name
    max_it = args.max_iter
    if faas is None and gene_clusterss is None:
        sys.exit("at least one of --faas and --gene_clusters_file is required")

    print(len(genomes), " genomes in your clade", file = sys.stderr)

    all_bins = [MetaBin(name = g, amino_acid_file = None if faas is None else faas[g], genome_completeness = checkm[g]) for g in genomes]

    motu = mOTU( genomes = all_bins, name = name, threads = threads, precluster = precluster, make_gene_clustering = False if args.load_gene_clusters else True)

    if args.load_gene_clusters:
        print("Loading imported gene-clusters")
        motu.load_gene_clusters(args.load_gene_clusters)

    motu.compute_core(max_it = max_it)

    if args.output:
        out_handle = open(out_json, "w")
    else :
        out_handle = sys.stdout
    stats = motu.get_stats()

    if args.save_gene_clusters:
        motu.export_gene_clusters(args.save_gene_clusters)

    if args.checkm:
        abs = lambda x : x if x > 0 else -x
        mean = lambda l : sum(l)/len(l)
        mean_complet_diff = mean([k.checkm_complet - k.new_completness for k in motu])
        max_complete_diff = max([k.checkm_complet - k.new_completness for k in motu])
        if mean_complet_diff > 5:
            print(f"**WARNING** : the mean difference between you prior and posterior completeness estimates is pretty high ({mean_complet_diff:.2f} %), This doesn't have to be a problem, but it could be due to a highly unbalanced set of genomes (e.g. many of one strain and few of an other), some genomes in the input set that shouldn't be there, or maybe your gene-clusters are too 'strict'...", file = sys.stderr)
        if max_complete_diff > 60:
            print(f"**WARNING** : the largest difference between you prior and posterior completeness estimates is pretty high ({max_complete_diff:.2f} %), This doesn't have to be a problem, but it could be due to a highly unbalanced set of genomes (e.g. many of one strain and few of an other), some genomes in the input set that shouldn't be there, or maybe your gene-clusters are too 'strict'...", file = sys.stderr)

    nb_boots = args.boots
    motu.roc_values(boots=nb_boots)

    if args.long and not args.genome2gene_clusters_only:
        json.dump(stats, out_handle, indent=4, sort_keys=True)
    else :
        print(motu.pretty_pan_table(), file = out_handle)
    if args.output:
        out_handle.close()

    return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog = "mOTUpan.py", description=description_text, epilog = "Let's do this")
    parser.add_argument('--output', '-o', nargs = None, help = "send output to this file")
    parser.add_argument('--force', '-f', action='store_true', help = "force execution answering default answers")
    parser.add_argument('--checkm', '-k',nargs = None, help = "checkm file if you want to seed completnesses with it, accepts concatenations of multiple checkm-files, check manual for more formating options")
    parser.add_argument('--seed', '-s', type = float , nargs = '?', help = "seed completeness, advice a number around 90 ({} default), this is the default completeness prior".format(__checkm_default__))
    parser.add_argument('--length_seed', '--ls', action='store_true', help = "seed completeness by length fraction [0-100]")
    parser.add_argument('--random_seed', '--rs', action='store_true', help = "random seed completeness between 5 and 95 percent")
    parser.add_argument('--save_gene_clusters', nargs = None, help = "saves the gene_clusters in an appropriate json-formated file")
    parser.add_argument('--precluster', action='store_true', default = False, help = "precluster proteomes with cd-hit, mainly for legacy reasons, mmseqs2 is faaaaaast")
    parser.add_argument('--faas','-F', nargs = '*', help = "list of amino-acids faas of MAGs or whatnot, or a text file with paths to the faas (with the --txt switch)")
    parser.add_argument('--txt', action='store_true', help = "the '--faas' switch indicates a file with paths")
    parser.add_argument('--load_gene_clusters', '--gene_clusters', '-c', nargs = '?', help = "file with COG-sets (see doc)")
    parser.add_argument('--name', '-n', nargs = None, help = "if you want to name this bag of bins")
    parser.add_argument('--long', action='store_true', help = "longform output, a JSON-file with a lot more information (might be cryptic...)")
    parser.add_argument('--boots', '-b', nargs = None, type = int , default = 0 , help = "number of bootstraps for fpr and recall estimate (default 0), careful, slows down program linearly")
    parser.add_argument('--max_iter', '-m', nargs = None, type = int , default = 20 , help = "max number of iterations (set to one if you have only few traits, e.g. re-estimation of completeness is nonsensical)")
    parser.add_argument('--threads', '-t', nargs = None, type = int , default = multiprocessing.cpu_count(), help = "number of threads (default all, e.g. {}), only gene clustering is multithreaded right now, the rest is to come".format(multiprocessing.cpu_count()))
    parser.add_argument('--version','-v', action="store_true", help = "get version")

    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        print("{script} Version {version}".format(script = __file__, version = __version__), file = sys.stderr)
        sys.exit(1)

    if args.version:
        print("{script} Version {version}".format(script = __file__, version = __version__), file = sys.stderr)
        sys.exit()


#    print(args, file=sys.stderr)

    motupan(args)


#for tt in `sed 's/\t/:/' scratch/test_data/mOTUs.txt` ;
#do
#    echo $tt | cut -f1 -d":"
#    fs=`echo $tt | cut -f2 -d":"| sed 's#;#.faa scratch/test_data/proteoms/#g'`;
#    mOTUlizer/bin/__main__.py -n `echo $tt | cut -f1 -d":"` --faas scratch/test_data/proteoms/${fs}.faa >> test;
#done
