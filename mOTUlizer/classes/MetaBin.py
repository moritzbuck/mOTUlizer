from os.path import join as pjoin
import tempfile
import sys
from subprocess import Popen, PIPE, call
import os
import shutil

class MetaBin:
    def __repr__(self) :
        return "< bin {name} with {n} cogs>".format(n = len(self.cogs) if cogs else "NA", name = self.name)

    def __init__(self, name, cogs,fnas, faas, complet = None, contamin = 0, max_complete = 99.9):
        self.name = name
        self.cogs = cogs if type(cogs) != str else set([cogs])
        self.faas = faas
        self.fnas = fnas
        self.checkm_complet = complet
        self.checkm_contamin = contamin
        if not self.checkm_complet is None:
            if self.checkm_complet > max_complete:
                self.checkm_complet =  max_complete
        self.new_completness = None

    def get_data(self):
        return { 'name' : self.name,
                 'faa-file' : self.faas,
                 'fna-file' : self.fnas,
                 'checkm_complet' : self.checkm_complet,
                 'checkm_contamin' : self.checkm_contamin,
                 'new_completness' : self.new_completness,
        }


    def overlap(self, target):
        return self.cogs.intersection(target.cogs)

    def estimate_nb_cogs(self):
        assert self.new_completness != None, "new_completness not computed, please do"
        return 100*len(self.cogs)/self.new_completness

    @classmethod
    def get_anis(cls, bins, outfile = None, method = "fastANI", block_size = 500, threads=1):
        if method == "fastANI":
            if not shutil.which('fastANI'):
                print("You need fastANI if you do not provide a file with pairwise similarities, either install it or provide pairwise similarities (see doc...)", file = sys.stderr)
                sys.exit(-1)
            fastani_file = tempfile.NamedTemporaryFile().name if outfile is None else outfile

            mags = [b.fnas for b in bins]

            mag_blocks = [mags[i:(i+block_size)] for i in list(range(0,len(mags), block_size))]

            if len(mag_blocks) > 1:
                print("You have more then {bsize} bins, so we will run fastANI in blocks, if it crashes due to memory, make smaller blocks".format(bsize = block_size), file=sys.stderr)

            with open(fastani_file, "w") as handle:
                handle.writelines(["query\tsubject\tani\tsize_q\tsize_s\n"])

            for i,bloc1 in enumerate(mag_blocks):
                b1_tfile = tempfile.NamedTemporaryFile().name

                with open(b1_tfile, "w") as handle:
                    handle.writelines([l +"\n" for l in bloc1])

                for j,bloc2 in enumerate(mag_blocks):
                        print("doing bloc {i} and {j}".format(i = i, j=j), file = sys.stderr)
                        b2_tfile = tempfile.NamedTemporaryFile().name
                        with open(b2_tfile, "w") as handle:
                            handle.writelines([l  +"\n" for l in bloc2])

                        out_tfile = tempfile.NamedTemporaryFile().name
                        call("fastANI --ql {b1} --rl {b2} -o {out} -t {threads} 2> /dev/null".format(b1 = b1_tfile, b2 = b2_tfile, out = out_tfile, threads = threads), shell = True)
                        with open(out_tfile) as handle:
                            new_dat = ["\t".join([ll for ll in l.split()]) +"\n" for l in handle.readlines()]
                        with open(fastani_file, "a") as handle:
                            handle.writelines(new_dat)

                        os.remove(out_tfile)
                        os.remove(b2_tfile)

            os.remove(b1_tfile)
            with open(fastani_file) as handle:
                handle.readline()
                out_dists = {(l.split()[0], l.strip().split()[1]) : float(l.split()[2]) for l in handle}
#                tfile = lambda k : ".".join(k.split(".")[:-1]) if (k.endswith(".fna") or k.endswith(".fa") or k.endswith(".fasta") or k.endswith(".fna") or k.endswith(".ffn")) else k
#                out_dists = {(tfile(k[0]),tfile(k[1])) : v for k,v in out_dists.items()}
            if outfile is None:
                os.remove(fastani_file)
        else :
            print("No other method for ani computation implemented yet")
            sys.exit()

        return out_dists
