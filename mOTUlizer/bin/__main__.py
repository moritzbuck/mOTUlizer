#!/usr/bin/env python

import os
import shutil
import sys
from os.path import join as pjoin
import argparse

print("This is temporary, fix the hard-path once all is clean", file=sys.stderr)
sys.path.append("/home/moritz/repos/moritz/0039_mOTUlizer/")

from mOTUlizer.classes import *


#from mOTUlizer.classes.mOTU import mOTU

#from mOTUlizer.config import *

description_text = """
From a buch of amino-acid sequences or COG-sets, computes concensus AA/COG sets.

Returns all to stdout by default.
"""

def main(args):
    faas = args['faas']
    cogs = args['cogs']
    out_json = args['output']
    checkm = args['checkm']

    return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog = "mOTUlizer", description=description_text, epilog = "Let's do this")
    parser.add_argument('--faas','-F', nargs = '*', help = "list of amino-acids faas of MAGs or whatnot")
    parser.add_argument('--cog_file', '--cogs', '-c', nargs = '?', help = "file with COG-sets (see doc)")
    parser.add_argument('--output', '-o', nargs = '?', help = "send output to this file")
    parser.add_argument('--force', '-f', help = "force execution answering default answers")
    parser.add_argument('--checkm', '-k',nargs = '?', help = "checkm file if you want to see completnesses with it")

    args = parser.parse_args()

    print(args, file=sys.stderr)

    main(args)
