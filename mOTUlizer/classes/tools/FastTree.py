from mOTUlizer.classes.tools import Wrapper
import tempfile, os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from ete3 import Tree

class FastTree(Wrapper):

    __executables__ = ['fasttree']

    def get_command(self):
        return f"fasttree {self.options} -out {self.out_file} {self.in_file}" + (" > /dev/null 2> /dev/null" if self.quiet else "")

    def __init__(self, seqsdict_or_infile, out_file = None, options = None , quiet = True):
        super().__init__()
        if options:
            self.options = options
        else :
            self.options = ""
        self.quiet = quiet

        if type(seqsdict_or_infile) == dict:
            self.temp_in = True
            self.in_file = tempfile.NamedTemporaryFile(delete = False).name
            SeqIO.write([SeqRecord(id = id_, seq = seq, description = "") for id_,seq in seqsdict_or_infile.items()], self.in_file, "fasta")
        else :
            self.temp_in = False
            self.in_file = seqsdict_or_infile
        if not out_file :
            self.temp_out = True
            self.out_file = tempfile.NamedTemporaryFile(delete = False).name
        else :
            self.temp_out = False
            self.out_file = out_file

    def parse_output(self):
        return Tree(self.out_file)

    def __del__(self):
        if self.temp_out:
            os.remove(self.out_file)
        if self.temp_in:
            os.remove(self.in_file)
