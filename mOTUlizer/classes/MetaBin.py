from os.path import join as pjoin
import tempfile
import sys
import os
import shutil
from mOTUlizer.config import FASTA_EXTS, MAX_COMPLETE
from mOTUlizer.errors import *
from Bio import SeqIO
from mOTUlizer.classes.GeneClusters import GeneClusters
from mOTUlizer.classes.GFF import GFF

class MetaBin:
    def __repr__(self) :
        return "< bin {name} with {n} gene_clusters >".format(n = len(self.gene_clustering) if self.gene_clustering else "NA", name = self.name)

    def __init__(self, name, nucleotide_file = None, amino_acid_file = None, gff_file = None, genome_completeness = 100, genome_redundancy = 0):
        self.can_haz_gene_clusters = True
        self.name = name
        self.gene_clusters = None
        self._genes = None
        self.amino_acid_file = amino_acid_file
        self.amino_acids = None
        self.nucleotide_file = nucleotide_file
        self._contigs = None
        self.gff_file = gff_file
        self.gff = None
        self.original_complet = genome_completeness
        self.original_redundancy = genome_redundancy
        if self.original_complet and self.original_complet > MAX_COMPLETE:
            self.original_complet =  MAX_COMPLETE
        self.new_completness = None
        self.gene_clustering = None

    def get_clust_set(self):
        return self.gene_clusters.get_clust_set()

    def set_gene_clusters(self, clusters):
        self.gene_clusters = GeneClusters(self, clusters)

    def get_aa(self, id, default = None):
        aa_dict = self.get_amino_acids()
        return aa_dict.get(id, default)

    def get_gene(self, id, default = None):
        gene_dict = self.get_genes()
        return gene_dict.get(id, default)

    def get_genes(self):
        if not self._genes:
            if not self.gff_file:
                raise CantGenesError("You need a gff if you want to use genes")
            if not self.nucleotide_file:
                raise CantGenesError("You need an nucleotide fasta if you want to use a gff")
            self._genes = { feat.get_id() : feat.get_gene() for feat in self.get_gff() if feat.feature == "CDS"}
        return self._genes


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
        if not self.gff:
            self.gff = GFF(self, self.gff_file)
        return self.gff

    def get_contigs(self):
        if not self.nucleotide_file:
            raise CantNucleotideError("can't parse a nucleotide if there ain't one")
        if not self._contigs:
            self._contigs = {s.id : s.seq for s in SeqIO.parse(self.nucleotide_file, "fasta")}
        return self._contigs

    def get_contig(self, contig_id):
        contigs = self.get_contigs()
        if contig_id not in contigs:
            raise CantNucleotideError("This contig ID ain't un your nucleotide data")
        else :
            return contigs[contig_id]

    def get_data(self):
        return { 'name' : self.name,
                 'faa-file' : self.amino_acid_file,
                 'fna-file' : self.nucleotide_file,
                 'original_complet' : self.original_complet,
                 'original_redundancy' : self.original_redundancy,
                 'new_completness' : self.new_completness,
        }


    def overlap(self, target):
        return self.gene_clustering.intersection(target.gene_clustering)

    def estimate_nb_gene_clusters(self):
        assert self.new_completness != None, "new_completness not computed, please do"
        return 100*len(self.gene_clustering)/self.new_completness
