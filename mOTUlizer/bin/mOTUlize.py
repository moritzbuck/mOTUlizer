#!/usr/bin/env python

import os
import shutil
import sys
from os.path import join as pjoin
import argparse
import json
from random import uniform

#print("This is temporary, fix the hard-path once all is clean", file=sys.stderr)
sys.path.append("/home/moritz/projects/0039_mOTUlizer/")

from mOTUlizer.classes import *
from mOTUlizer.utils import *
from mOTUlizer.classes.MetaBin import MetaBin
from mOTUlizer.classes.mOTU import mOTU

#from mOTUlizer.config import *

description_text = """
From a set of genomes, makes metagenomic Operational Taxonomic Units (mOTUs). By default it makes a graph of 95%
(reciprocal) ANI (with fastANI) connected MAGs (with completeness > 40%, contamination < 5%). The mOTUs will be the
connected components of this graph, to which smaller "SUBs" with ANI > 95%  are
recruited.

If similarities provided, it should be a TAB-separated file with columns as query, subject and similarity (in percent, e.g. [0-100])
if you also provide fasta-files (for stats purpouses) query and names should correspond to the fasta-files you provide.
If the columns are file names, the folders are removed (mainly so it can read fastANI output directly).

"""


def main(args):
    #parse and check your amino-acid files

    if args.txt and args.fnas:
        with open(args.fnas[0]) as handle:
            fnas = {os.path.splitext(os.path.basename(f.strip().rstrip(".gz")))[0] : f.strip() for f in handle.readlines()}
    elif args.fnas:
        fnas = {os.path.splitext(os.path.basename(f.strip().rstrip(".gz")))[0] : f for f in args.fnas}
    else :
        fnas = {}

    ani_cutoff = args.similarity_cutoff
    similarities = args.similarities
    checkm_file = args.checkm
    out_json = args.output
    prefix = args.prefix
    mag_complete = args.MAG_completeness
    mag_contamin = args.MAG_contamination
    sub_complete = args.SUB_completeness
    sub_contamin = args.SUB_contamination
    threads = args.cpus
    keep_simi = args.keep_simi_file

    assert 0 < ani_cutoff < 100, "similarity cutoff needs to be between 0 and 100 (percent similarity)"
    assert all([os.path.exists(f) for f in fnas.values()]), "one or some of your fnas don't exists"
    assert not (len(fnas) == 0 and similarities is None), "you need to give at least one "

    if similarities:
        assert os.path.exists(similarities), "The file for similarities does not exists"
    if checkm_file:
        assert os.path.exists(checkm_file), "The file for checkm does not exists"

    print("Parsing the checkm-file")

    checkm_info = parse_checkm(checkm_file)

    if similarities:
        print("Loading similarities")
        dist_dict = {}
        with open(similarities) as handle:
            for l in handle:
                if "query" not in l:
                    ll = l.split("\t")
                    g1 = ".".join(os.path.basename(ll[0]).split(".")[:-1]) if ll[0].endswith(".fna") else ll[0]
                    g2 = ".".join(os.path.basename(ll[1]).split(".")[:-1]) if ll[1].endswith(".fna") else ll[1]
                    dist = float(ll[2])
                    dist_dict[(g1,g2)] = dist
    else :
        dist_dict = None

    if len(fnas) > 0 :
        genomes = set(fnas.keys())
    else :
        genomes = {a for k in dist_dict for a in k}

    if dist_dict:
        print("{nb_gen} genomes found with {counts} ANI mOTU edges".format(nb_gen = len(genomes), counts = len(dist_dict)), file=sys.stderr)
    else :
        print("{nb_gen} genomes found".format(nb_gen = len(genomes)), file=sys.stderr)



    assert all([g in checkm_info  for g in genomes]), "you do not have completness/contamination info for all you bins, values missing for :" + ", ".join([g for g in genomes if g not in checkm_info][0:10]) + "... (only 10 first shown   )"

    if fnas == {}:
        fnas = {g : None for g in genomes}

    print("making bin-objects")

    all_bins = [MetaBin(name = g, cogs = None, faas = None, fnas = fnas[g], complet = checkm_info[g]['Completeness'], contamin = checkm_info[g]['Contamination'], max_complete = 100) for g in genomes]

    if dist_dict is None:
        print("Similarities not provided, will compute them with fastANI", file = sys.stderr)
        dist_dict = MetaBin.get_anis(all_bins, threads = threads, outfile = keep_simi)

    print("making mOTUs")

    mOTUs = mOTU.cluster_MetaBins(all_bins, dist_dict, ani_cutoff, prefix, mag_complete, mag_contamin, sub_complete, sub_contamin)

    print("making stats")

    out_dict = {}
    for m in mOTUs:
        out_dict.update(m.get_stats())

    if args.output:
        out_handle = open(out_json, "w")
    else :
        out_handle = sys.stdout
    json.dump(out_dict, out_handle, indent=4, sort_keys=True)
    if args.output:
        out_handle.close()

    return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog = "mOTUlize", description=description_text, epilog = "Let's do this")
    parser.add_argument('--output', '-o', nargs = '?', help = "send output to this file")
    parser.add_argument('--force', '-f', action='store_true', help = "force execution answering default answers")
    parser.add_argument('--checkm', '-k',nargs = '?', help = "checkm file (or whatever you want to use as completness)", required=True)
    parser.add_argument('--similarities', '-I', nargs = '?', help = "file containing similarities between MAGs, if not provided, will use fastANI to compute one")
    parser.add_argument('--fnas','-F', nargs = '*', help = "list of nucleotide fasta-files of MAGs or whatnot")
    parser.add_argument('--prefix', '-n', nargs = '?', default = "mOTU_", help = "prefix for the mOTU names, default : mOTU_ ")
    parser.add_argument('--MAG-completeness', '--MC', '-M', nargs = '?', type=float, default = 40, help = "completeness cutoff for seed MAGs, default : 40")
    parser.add_argument('--MAG-contamination', '--Mc', '-m', nargs = '?', type=float, default = 5, help = "contamination cutoff for seed MAGs, default : 5")
    parser.add_argument('--SUB-completeness', '--SC', '-S', nargs = '?', type=float, default = 0, help = "completeness cutoff for recruited SUBs, default : 0")
    parser.add_argument('--SUB-contamination', '--Sc', '-s', nargs = '?', type=float, default = 100, help = "contamination cutoff for recruited SUBs, default : 100")
    parser.add_argument('--similarity-cutoff', '-i', nargs = '?', type=float, default = 95, help = "distance cutoff for making the graph, default : 95")
    parser.add_argument('--cpus', '-c', nargs = '?', type=int, default = 1, help = "number of threads, default : 1")
    parser.add_argument('--keep-simi-file', '-K', nargs = '?', default = None, help = "keep generated similarity file if '--similarities' is not procided")
    parser.add_argument('--txt', '-t', action='store_true', help = "the '--fnas' switch indicates a file with paths")

    args = parser.parse_args()

#    print(args, file=sys.stderr)

    main(args)
