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

from mOTUlizer import __version__
from mOTUlizer.classes import *
from mOTUlizer.utils import *
from mOTUlizer.classes.mOTU import mOTU

#from mOTUlizer.config import *

description_text = """
From a buch of amino-acid sequences or COG-sets, computes concensus AA/COG sets.

Returns all to stdout by default.
"""
__checkm_default__ = 95

def main(args):
    if args.version:
        print("{script} Version {version}".format(script = __file__, version = __version__))
        sys.exit()

    if args.cog_file:
        try :
            if args.cog_file.endswith(".json") or args.cog_file.endswith(".gid2cog"):
                with open(args.cog_file) as handle:
                    cog_dict = json.load(handle)
            else :
                with open(args.cog_file) as handle:
                    cog_dict = {l.split("\t")[0] : l[:-1].split("\t")[1:] for l in handle}
            cog_dict = {k : set(v) for k,v in cog_dict.items()}
        except :
            print("Either the cog_file does not exists or it is not formated well")

        if all([len(v) == 0 for v in cog_dict]):
            print("None of your bins have any cogs in them, that's weird, you probably have wrong delimiter in you file, use tab.\nIf you do not have COGs you can also just run it without the --cog_file option and mOTUlizer will automatically compute some! (slower)")
    else :
        cog_dict = {}

    #parse and check your amino-acid files
    if args.txt and args.faas:
        with open(args.faas[0]) as handle:
            faas = {os.path.splitext(os.path.basename(f.strip().rstrip(".gz")))[0] : f.strip() for f in handle.readlines()}
    elif args.faas:
        faas = {os.path.splitext(os.path.basename(f.strip().rstrip(".gz")))[0] : f for f in args.faas}
    else :
        faas = {}

    assert all([os.path.exists(f) for f in faas.values()]), "one or some of your faas don't exists"


    if len(faas) > 0:
        genomes = set(faas.keys()).intersection(set(cog_dict.keys()))
    else :
        genomes = set(cog_dict.keys())
    print("len cog dict" , len(cog_dict))

    if cog_dict and len(faas) > 0:
        if len(genomes) != len(faas) or len(faas) != len(cog_dict):
            print("your faas and cog_drct are not the same length,\nit might not matter just wanted to let you know.", file = sys.stderr)

    if len(cog_dict) > 0 :
        cog_dict = {g : cog_dict.get(g) for g in genomes if g in cog_dict}

    out_json = args.output
    checkm = {}
    if args.checkm :
        assert os.path.exists(args.checkm), "The file for checkm does not exists"

        print("Parsing the checkm-file")
        checkm = {k : v['Completeness'] for k,v in parse_checkm(args.checkm).items()}
        checkm = {g : checkm.get(g) for g in genomes}
        if not all([v for v in checkm.values()]):
            print("you did not have completeness for all bins/genomes, abscent values have been put to default value, e.g. " + str(__checkm_default__))
            checkm = {k : v if v  else __checkm_default__ for k,v in checkm.items()}
    if args.seed :
        for f in genomes:
            checkm[f] = args.seed
    if args.random_seed :
        for f in genomes:
            checkm[f] = uniform(50,80)
    if args.length_seed :
        checkm = "length_seed"


    name = args.name if args.name else random_name()
    max_it = args.max_iter
    if faas is None and cogs is None:
        sys.exit("at least one of --faas and --cog_file is required")

    print(len(genomes), " genomes in your clade")

    motu = mOTU( name = name , faas = faas , cog_dict = cog_dict, checkm_dict = checkm, max_it = max_it)

    if args.output:
        out_handle = open(out_json, "w")
    else :
        out_file = sys.stdout
    if not args.genome2cog_only:
        json.dump(motu.get_stats(), out_handle, indent=4, sort_keys=True)
    else :
        json.dump({k : list(v) for k,v in motu.cog_dict.items()}, out_handle)
    if args.output:
        out_handle.close()

    return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog = "mOTUpan.py", description=description_text, epilog = "Let's do this")
    parser.add_argument('--output', '-o', nargs = '?', help = "send output to this file")
    parser.add_argument('--force', '-f', action='store_true', help = "force execution answering default answers")
    parser.add_argument('--checkm', '-k',nargs = '?', help = "checkm file if you want to see completnesses with it")
    parser.add_argument('--seed', '-s', type = float , nargs = '?', help = "seed completeness, advice a number around 90 ({} default)".format(__checkm_default__))
    parser.add_argument('--length_seed', '--ls', action='store_true', help = "seed completeness by length fraction [0-100]")
    parser.add_argument('--random_seed', '--rs', action='store_true', help = "random seed completeness between 5 and 95 percent")
    parser.add_argument('--genome2cog_only', action='store_true', help = "returns genome2cog only")
    parser.add_argument('--faas','-F', nargs = '*', help = "list of amino-acids faas of MAGs or whatnot, or a text file with paths to the faas (with the --txt switch)")
    parser.add_argument('--txt', '-t', action='store_true', help = "the '--faas' switch indicates a file with paths")
    parser.add_argument('--cog_file', '--cogs', '-c', nargs = '?', help = "file with COG-sets (see doc)")
    parser.add_argument('--name', '-n', nargs = '?', help = "if you want to name this bag of bins")
    parser.add_argument('--max_iter', '-m', nargs = '?', type = int , default = 20 , help = "if you want to name this bag of bins")
    parser.add_argument('--version','-v', action="store_true", help = "get version")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

#    print(args, file=sys.stderr)

    main(args)


#for tt in `sed 's/\t/:/' scratch/test_data/mOTUs.txt` ;
#do
#    echo $tt | cut -f1 -d":"
#    fs=`echo $tt | cut -f2 -d":"| sed 's#;#.faa scratch/test_data/proteoms/#g'`;
#    mOTUlizer/bin/__main__.py -n `echo $tt | cut -f1 -d":"` --faas scratch/test_data/proteoms/${fs}.faa >> test;
#done
