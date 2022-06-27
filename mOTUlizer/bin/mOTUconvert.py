#!/usr/bin/env python

import os
import shutil
import sys
from os.path import join as pjoin
import argparse
import json
from mOTUlizer import __version__

from mOTUlizer.utils import *
from mOTUlizer.classes.Parser import EmapperParse, PPanGGolinParse, RoaryParse, MmseqsParse, AnvioParse, OldEmapperParse

description_text = """
Converts output of diverse software generatig COGs, or genetically encoded traits into a genome2gene_clusters-JSON file useable by mOTUpan
"""

list_text = """
emapper/old_emapper :
\tdescribes each genomes as a set of eggNOGs based on the output of eggNOGmapper,
\tonly the deepest (e.g. closest to the root) eggNOG for each hit is taken,
\tinput-file is just the  '.emapper.annotations'-file.
\tThe old_emapper is for older versions of eggnoggmapper can be spoted by the
\tfirst column being #query_name not #query

ppanggolin :
\textracts the gene families from the hdf5-file (e.g. '.h5' file in the output folder)

roary :
\tjust converts the cluster output file, clustered_proteins by default, the file
\tof the '-o' switch, to a json file readable by mOTUpan

mmseqs2 :
\tjust converts the cluster output file, <o:cluster_Prefix>_cluster.tsv by default, generated
\twhen running mmseq2 easy-cluster, where <o:cluster_Prefix> is the output prefix you give to
\tthe command.

anvi'o :
\textracts the gene-clusters from a pangenome-database

"""

def motuconvert(args):

    kwargs = dict()
    if args.gene2genome:
        print("Parsing the --gene2genome file", file = sys.stderr )
        with open(args.gene2genome) as handle:
            kwargs['gene_id2genome'] = {l.split("\t")[0] : l.split("\t")[1].strip().split(";") for l in handle}

    method = args.in_type
    if method == "emapper":
        converter = EmapperParse(**kwargs)
    elif method == "old_emapper":
        converter = OldEmapperParse(**kwargs)
    elif method == "ppanggolin":
        converter = PPanGGolinParse(**kwargs)
    elif method == "roary":
        converter = RoaryParse(**kwargs)
    elif method == "mmseqs2":
        converter = MmseqsParse(**kwargs)
    elif method == "anvio":
        converter = AnvioParse(**kwargs)
    else :
        print("This is not implemented yet run with '--list' to see available options", file = sys.stderr )
        sys.exit(0)


    print("Doing the convertion", file = sys.stderr )

    genome2gene_clusterss = converter.convert(infile = args.input, count = args.count)

    if args.output:
        out_handle = open(args.output, "w")
    else :
        out_handle = sys.stdout

    json.dump(genome2gene_clusterss, out_handle, indent=4, sort_keys=True)

    if args.output:
        out_handle.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog = "mOTUconverts.py", description=description_text, epilog = "Let's do this")
    parser.add_argument('--output', '-o', nargs = '?', help = "send output to this file defaults to stdout")
    parser.add_argument('input', nargs = '?', metavar = "INPUT", help = "input file(s), check '--list' for specifics")
    parser.add_argument('--gene2genome', nargs='?', help = "if gene names not '${genome_name}_[0-9]*', a tab separated file with id of gene in the fist column and a semi-column separated second column containing genomes_id of genomes containing it")
    parser.add_argument('--in_type', '-I', nargs = '?', default = "emapper" , help = "software generating the input")
    parser.add_argument('--version','-v', action="store_true", help = "get version and exit")
    parser.add_argument('--list','-l', action="store_true", help = "list tools available and exit")
    parser.add_argument('--count','-c', action="store_true", default = False, help = "count the occurences of COG/trait")

    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help()
        print("{script} Version {version}".format(script = __file__, version = __version__), file = sys.stderr)
        sys.exit(1)

    if args.version:
        print("{script} Version {version}".format(script = __file__, version = __version__), file = sys.stderr)
        sys.exit()

    if args.list:
        print(list_text)
        print("{script} Version {version}".format(script = __file__, version = __version__), file = sys.stderr)
        sys.exit()

    motuconvert(args)
