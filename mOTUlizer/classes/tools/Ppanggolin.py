from mOTUlizer.classes.tools import Wrapper
from mOTUlizer.errors import CantGFFError, AlreadyAnnotatedError
import tempfile, os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from mOTUlizer.utils import message, parse_gff
from mOTUlizer import get_quiet, get_threads
import shutil
import tempfile
from os.path import join as pjoin
import h5py
import hdf5plugin
from mOTUlizer.db.SeqDb import SeqDb
import sys
from subprocess import call

class Ppanggolin(Wrapper):

    __executables__ = ['ppanggolin']
    __source__ = "prokka v1.13"

    def get_command(self):
        quiet = " > /dev/null 2> /dev/null" if get_quiet() else ""
        command = f"""
        ppanggolin annotate --anno {self.out_folder}/gff_list  -o {self.out_folder}/ppangg --basename cluster --tmpdir {self.out_folder} -c {self.threads} {quiet}
        ppanggolin cluster -p {self.out_folder}/ppangg/cluster.h5 -c {self.threads} --tmpdir {self.out_folder} {quiet}
        """
        if self.function == "cluster":
            return command
        elif self.function == "partition":
            return command + self.get_partitioning_command()

    def get_partitioning_command(self):
        quiet = " > /dev/null 2> /dev/null" if get_quiet() else ""
        return f"""
        ppanggolin graph -p {self.out_folder}/ppangg/cluster.h5 -c {self.threads} --tmpdir {self.out_folder} {quiet}
        ppanggolin partition -K {self.k} -p {self.out_folder}/ppangg/cluster.h5 -c {self.threads} --tmpdir {self.out_folder} {quiet}
        """

    def run_just_partitioning(self):
        try :
            call(self.get_partitioning_command(), shell = True)
        except :
            raise CantRunError(f"Some error happened while executing an {type(self)} object.\nThe command it was trying to run was {self.get_partitioning_command()}")


    def __init__(self, motu, folder = None, function = "cluster"):
        super().__init__()
        self.db = SeqDb.get_global_db()
        self.motu = motu
        self.temp =  not folder
        self.k = 3
        if not folder:
            self.out_folder = tempfile.mkdtemp()
        else :
            self.out_folder = folder
        self.function = function

        os.makedirs(pjoin(self.out_folder, "gffs"), exist_ok=True)

        for genome in motu:
            genome.export_genome_gff(path = pjoin(self.out_folder, "gffs", genome.name + ".gff"), CDS_only = True)

        with open(pjoin(self.out_folder,"gff_list"), "w") as handle:
            handle.writelines([f"{genome.name}\t{self.out_folder}/gffs/{genome.name}.gff\n" for genome in motu])

    def remove_superpang(self):
        superpang = list(self.db.get_genome_feats([m.name for m in self.motu if "_SuperPang" in m.name][0]).keys())
        handle = h5py.File(pjoin(self.out_folder, "ppangg", "cluster.h5"), "a")
        tt = handle['geneFamilies'][:]
        new_data = tt[ [int(v[0].decode()) not in superpang for v in tt ] ]
        del handle['geneFamilies']
        handle.create_dataset('geneFamilies', data = new_data)

        survivings_clsts = {t[1].decode() for t in new_data}
        survivings_genes = {t[0].decode() for t in new_data}

        tt2 = handle['geneFamiliesInfo'][:]
        new_info = tt2[ [v[0].decode() in survivings_clsts for v in tt2 ] ]
        del handle['geneFamiliesInfo']
        handle.create_dataset('geneFamiliesInfo', data = new_info)

        tt3 = handle['annotations/genes'][:]
        new_genes = tt3[ [v[1][0].decode() in survivings_genes for v in tt3 ] ]
        del handle['annotations/genes']
        handle.create_dataset('annotations/genes', data = new_genes)
        handle.close()


    def parse_output(self):
        if not self.quiet:
            print("Parse and decode gene-clusters from the h5py-file", file = sys.stderr)
        with h5py.File(pjoin(self.out_folder, "ppangg", "cluster.h5"), "r") as handle:
            gene_id2family = {a.decode() : b.decode() for a,b in handle['geneFamilies']}
            clst2partition = { v[0].decode() : v[1].decode() for v in handle['geneFamiliesInfo']}
        clst2feature_id = {k : [] for k in set(gene_id2family.values())}
        for k,v in gene_id2family.items():
            clst2feature_id[v] += [k]
        clst_name2clst_id = {k : i for i, k in enumerate(clst2feature_id.keys())}
        clst2rep = { clst_name2clst_id[k] : k for k,v in clst2feature_id.items()}
        clst2feature_id = { clst_name2clst_id[k] : v for k,v in clst2feature_id.items()}
        clst2partition = { clst_name2clst_id[k] : v for k,v in clst2partition.items()}

        return {'clst2feature_id' : clst2feature_id, 'clst2rep' : clst2rep, 'clst2partition' : clst2partition}

    def __del__(self):
        if self.temp and os.path.exists(self.out_folder):
            shutil.rmtree(self.out_folder)
