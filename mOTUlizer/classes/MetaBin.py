from os.path import join as pjoin
import tempfile
import sys
from subprocess import Popen, PIPE, call
import os
import shutil
from mOTUlizer.config import FASTA_EXTS, MAX_COMPLETE
from mOTUlizer.errors import *
from Bio import SeqIO

class MetaBin:
    def __repr__(self) :
        return "< bin {name} with {n} gene_clusters >".format(n = len(self.gene_clustering) if self.gene_clustering else "NA", name = self.name)

    def __init__(self, name,nucleotide_file = None, amino_acid_file = None, gff_file = None, complet = 100, contamin = 0):
        self.name = name
        self.gene_clustering = None
        self.amino_acid_file = amino_acid_file
        self.amino_acids = None
        self.nucleotide_file = nucleotide_file
        self.contigs = None
        self.gff_file = gff_file
        self.gff = None
        self.original_complet = complet
        self.original_contamin = contamin
        if self.original_complet > MAX_COMPLETE:
            self.original_complet =  MAX_COMPLETE
        self.new_completness = None

        self.faas = None

    def get_amino_acids(self):
        if not self.amino_acids:
            if not self.amino_acid_file:
                if not self.gff_file:
                    raise CantAminoAcidsError("You need either a gff or a amino-acid fasta if you want to use amino-acids")
                if not self.nucleotide_file:
                    raise CantAminoAcidsError("You need an nucleotide fasta if you want to use a gff to get amino-acids")
                self.amino_acids = { feat.get_id() : feat.get_amino_acids() for feat in self.get_gff() if feat.feature == "CDS"}
            else :
                self.amino_acids = {s.id : s.seq for s in SeqIO.parse(self.amino_acid_file, "fasta")}
        return self.amino_acids

    def get_gff(self):
        if not self.gff_file:
            CantGFFError("can't parse a gff if there ain't one")
        self.gff = GFF(self, self.gff_file)

    def get_contigs(self):
        if not self.nucleotide_file:
            CantNucleotideError("can't parse a nucleotide if there ain't one")
        self.gff = GFF(self, self.gff_file)

    def get_data(self):
        return { 'name' : self.name,
                 'faa-file' : self.amino_acid_file,
                 'fna-file' : self.nucleotide_file,
                 'original_complet' : self.original_complet,
                 'original_contamin' : self.original_contamin,
                 'new_completness' : self.new_completness,
        }


    def overlap(self, target):
        return self.gene_clustering.intersection(target.gene_clustering)

    def estimate_nb_gene_clusters(self):
        assert self.new_completness != None, "new_completness not computed, please do"
        return 100*len(self.gene_clustering)/self.new_completness
