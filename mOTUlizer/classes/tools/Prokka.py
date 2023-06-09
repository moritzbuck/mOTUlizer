from mOTUlizer.classes.tools import Wrapper
from mOTUlizer.errors import CantGFFError, AlreadyAnnotatedError
import tempfile, os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from mOTUlizer.utils import message, parse_gff
from mOTUlizer import get_quiet
import shutil
import mOTUlizer

class Prokka(Wrapper):

    __executables__ = ['prokka']
    __source__ = "prokka v1.13"

    def get_command(self):
        quiet = " > /dev/null 2> /dev/null" if mOTUlizer._quiet_ else ""
        if self.annotate:
            annot = ""
        else :
            annot = "--noanno"
        return f"prokka --outdir {self.out_file} --prefix {self.prefix} --locustag {self.prefix} --cpus {self.threads} {annot} {self.fna}  {quiet}"

    def __init__(self, metabin, threads = 20, annotate = True):
        super().__init__()

        if metabin.has_features():
            raise AlreadyAnnotatedError(f"You're trying to reannotate {metabin.name} frist remove previous gene-calls")
        self.fna = metabin.nucleotide_file
        self.prefix = metabin.name
        self.out_file = metabin.temp_folder + "/prokka"
        self.threads = 24
        self.annotate = annotate

    def parse_output(self):
        if not os.path.exists(self.out_file + "/" + self.prefix + ".gff"):
            raise CantGFFError(f"Your prokka run of {self.prefix} somehow didn't work")
        else :
            features = parse_gff(self.out_file + "/" + self.prefix + ".gff")
            feat_source =  "prokka v1.13"
        return features

    def __del__(self):
        if os.path.exists(self.out_file):
            shutil.rmtree(self.out_file)
