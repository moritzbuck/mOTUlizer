from mOTUlizer.classes.tools import Wrapper
import tempfile, os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from ete3 import Tree
from mOTUlizer.utils import parse_fasta
from os.path import join as pjoin
import shutil

class SuperPang(Wrapper):

    __executables__ = ['SuperPang.py']

    def get_command(self):
        files = " ".join([g.nucleotide_file for g in self.motu])
        self.motu.write_checkm_file(self.temp_checkm, for_temp_files = True)
        return f"SuperPang.py --force-overwrite {self.options} -i {self.identity} -q {self.temp_checkm} -f {files} -t {self.threads} -o {self.out_file}" + (" > /dev/null 2> /dev/null" if self.quiet else "")

    def __init__(self, motu, out_file = None, identity = 0.75, full_assembly = False, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.motu = motu
        self.full_assembly = full_assembly
        self.identity = identity
        self.temp_checkm = tempfile.NamedTemporaryFile(delete = False).name

        if not out_file :
            self.temp_out = True
            self.out_file = pjoin(tempfile.mkdtemp(), "superpang_out")
        else :
            self.temp_out = False
            self.out_file = out_file

    def parse_output(self):
        if self.full_assembly:
            seqs = parse_fasta(pjoin(self.out_file, "assembly.fasta"))
        else :
            seqs = parse_fasta(pjoin(self.out_file, "NBPs.fasta"))
        return seqs

    def __del__(self):
        if self.temp_out:
            shutil.rmtree(self.out_file)
        os.remove(self.temp_checkm)
