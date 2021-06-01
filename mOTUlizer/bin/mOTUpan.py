#!/usr/bin/env python

import os
import shutil
import sys
from os.path import join as pjoin
import argparse
import json
from random import uniform

#print("This is temporary, fix the hard-path once all is clean", file=sys.stderr)
#sys.path.append("/home/moritz/projects/0039_mOTUlizer/")

from mOTUlizer import __version__
from mOTUlizer.classes import *
from mOTUlizer.utils import *
from mOTUlizer.classes.mOTU import mOTU
import multiprocessing

#from mOTUlizer.config import *

description_text = """
From a buch of amino-acid sequences or COG-sets, computes concensus AA/COG sets.

Returns all to stdout by default.
"""
__checkm_default__ = 95

def motupan(args):
    if args.cog_file:
        try :
            with open(args.cog_file) as handle:
                cog_dict = json.load(handle)
        except :
            try :
                with open(args.cog_file) as handle:
                    cog_dict = {l.split("\t")[0] : l[:-1].split("\t")[1:] for l in handle}
            except :
                print("Either the cog_file does not exists or it is not formated well", file = sys.stderr)
        cog_dict = {k : set(v) for k,v in cog_dict.items()}
        if all([len(v) == 0 for v in cog_dict]):
            print("None of your bins have any cogs in them, that's weird, you probably have wrong delimiter in you file, use tab.\nIf you do not have COGs you can also just run it without the --cog_file option and mOTUlizer will automatically compute some!", file = sys.stderr)
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


    if len(faas) > 0 and args.cog_file:
        genomes = set(faas.keys()).intersection(set(cog_dict.keys()))
    elif len(faas) > 0:
        genomes = set(faas.keys())
    else :
        genomes = set(cog_dict.keys())

    if cog_dict and len(faas) > 0:
        if len(genomes) != len(faas) or len(faas) != len(cog_dict):
            print("your faas and cog_drct are not the same length,\nit might not matter just wanted to let you know.", file = sys.stderr)

    if len(cog_dict) > 0 :
        cog_dict = {g : cog_dict.get(g) for g in genomes if g in cog_dict}

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
        checkm = "length_seed"
    elif args.seed :
        for f in genomes:
            checkm[f] = args.seed
    elif args.random_seed :
        for f in genomes:
            checkm[f] = uniform(50,80)
    else :
        for f in genomes:
            checkm[f] = __checkm_default__

    threads = args.threads
    precluster = args.precluster
    method = None if args.genome2cog_only else "default"

    name = args.name if args.name else random_name()
    name = (name + "_") if not name.endswith("_") else name
    max_it = args.max_iter
    if faas is None and cogs is None:
        sys.exit("at least one of --faas and --cog_file is required")

    print(len(genomes), " genomes in your clade", file = sys.stderr)

    motu = mOTU( name = name , faas = faas , cog_dict = cog_dict, checkm_dict = checkm, max_it = max_it, threads = threads, precluster = precluster, method = method)

    if args.output:
        out_handle = open(out_json, "w")
    else :
        out_handle = sys.stdout
    stats = motu.get_stats()

    nb_boots = args.boots
    if args.long and not args.genome2cog_only:

        stats[name].update(motu.roc_values(boots=nb_boots))
        json.dump(stats, out_handle, indent=4, sort_keys=True)
    elif args.genome2cog_only :
        json.dump({k : list(v) for k,v in motu.cog_dict.items()}, out_handle, indent=4, sort_keys=True)
    else :
        out_dict = {}

        stats = list(stats.values())[0]
        stats.update(motu.roc_values(boots = nb_boots))
        cogs = set([cc for g,c in stats['cogs'].items() for cc in c]) if 'aa' not in stats['cogs'] else set(stats['cogs']['aa'].values())
        for k in cogs:
            out_dict[k] = {}
            out_dict[k]['type'] = 'core' if k in stats['core'] else 'accessory'
            out_dict[k]['genome_occurences'] = 0
            out_dict[k]['log_likelihood_to_be_core'] = stats['likelies'][k]
            out_dict[k]['genomes'] = []
            out_dict[k]['genes'] = [] if 'aa' in stats['cogs'] else ["NA"]
            out_dict[k]['trait_name'] = k
        if 'aa' in stats['cogs']:
            for k,v in stats['cogs']['aa'].items():
                out_dict[v]['genes'] += [k]
        for k,v in stats['cogs'].items() if 'aa' not in stats['cogs'] else stats['cogs']['genome'].items():
            for vv in v:
                out_dict[vv]['genomes'] += [k]
                out_dict[vv]['genome_occurences'] += 1

        for k,v in out_dict.items():
            v['mean_copy_per_genome'] = "NA" if not v['genes'] else len(v['genes'])/len(v['genomes'])
            v['genes'] = ";".join(v['genes'])
            v['genomes'] = ";".join(v['genomes'])

        header = ['trait_name','type', 'genome_occurences', 'log_likelihood_to_be_core', 'mean_copy_per_genome','genomes', 'genes']
        genome_line = "genomes=" + ";".join(["{}:prior_complete={}:posterior_complete={}".format(k['name'], k['checkm_complet'], k['new_completness']) for k in stats['genomes']])
        mean = lambda l : sum([ll for ll in l])/len(l)

        if stats['mean_recall'] != "NA":
            bootsy = """
#bootstrapped_mean_false_positive_rate={fpr:.2f};bootstrapped_sd_false_positive_rate={sd_fpr:.2f}
#bootstrapped_mean_recall={recall:.2f};bootstrapped_sd_recall={sd_recall:.2f}
#bootstrapped_mean_lowest_false_positive={lowest:.2f};bootstrapped_sd_lowest_false_positive={sd_lowest:.2f}
#bootstrapped_nb_reps={boots}
#""".format( boots=nb_boots,
            fpr=stats['mean_fpr'],
            recall = stats['mean_recall'], lowest = stats['mean_lowest_false'],sd_fpr=stats['sd_fpr'],
            sd_recall = stats['sd_recall'], sd_lowest = stats['sd_lowest_false'])
        else :
            bootsy=""


        outformat ="""#mOTUlizer:mOTUpan:{version}
#run_name={name}
#
#genome_count={nb}
#core_length={core_len}
#mean_prior_completeness={prior_complete:.2f}
#mean_posterior_completeness={post_complete:.2f}
#sum_abs_loglikelihood_ratios={SALLHR:.2f}
#mean_est_genome_size={size:.2f};traits_per_genome
#{genomes}
#{boostrap}
{header}
{data}
"""
        print(outformat.format(version = __version__ , nb = stats['nb_genomes'],
            name = motu.name.strip("_"),
            core_len = len(motu.core),
            genomes=genome_line,
            prior_complete=mean([b.checkm_complet for b in motu]),
            post_complete=mean([b.new_completness for b in motu]),
            SALLHR=sum([l if l > 0 else -l for l in motu.likelies.values()]),
            size=mean([100*len(b.cogs)/b.new_completness for b in motu]),
            boostrap=bootsy,
            header = "\t".join(header), data = "\n".join(["\t".join([str(v[hh]) for hh in header]) for v in out_dict.values()])),
            file=out_handle)

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
    parser.add_argument('--genome2cog_only', action='store_true', help = "returns genome2cog only")
    parser.add_argument('--precluster', action='store_true', default = False, help = "precluster proteomes with cd-hit, mainly for legacy reasons, mmseqs2 is faaaaaast")
    parser.add_argument('--faas','-F', nargs = '*', help = "list of amino-acids faas of MAGs or whatnot, or a text file with paths to the faas (with the --txt switch)")
    parser.add_argument('--txt', action='store_true', help = "the '--faas' switch indicates a file with paths")
    parser.add_argument('--cog_file', '--cogs', '-c', nargs = '?', help = "file with COG-sets (see doc)")
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
