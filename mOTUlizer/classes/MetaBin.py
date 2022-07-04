from os.path import join as pjoin
import tempfile
import sys
import os
import shutil
import gzip
from mOTUlizer.config import FASTA_EXTS, MAX_COMPLETE
from mOTUlizer.errors import *
from Bio import SeqIO
from mOTUlizer.classes.GeneClusters import GeneClusters
from mOTUlizer.classes.GFF import GFF

class MetaBin:
    def __repr__(self) :
        return "< bin {name} with {n} gene_clusters >".format(n = len(self.gene_clusters) if self.gene_clusters else "NA", name = self.name)

    def __init__(self, name, nucleotide_file = None, amino_acid_file = None, gff_file = None, genome_completeness = None, genome_redundancy = 0, gbk_file = None):
        self.can_haz_gene_clusters = True
        self.name = name
        self.gene_clusters = None
        self._genes = None
        self.amino_acid_file = amino_acid_file
        self.amino_acids = None
        self.nucleotide_file = nucleotide_file
        self.gbk_file =  gbk_file
        if self.gbk_file :
            _ = self.get_nucleotide_file()
        self._contigs = None
        self.gff_file = gff_file
        self.gff = None
        self._original_complet = genome_completeness
        self.original_redundancy = genome_redundancy
        if self._original_complet and self._original_complet > MAX_COMPLETE:
            self._original_complet =  MAX_COMPLETE
        self.new_completness = None
        self.gene_clusters = None

    def __del__(self):
        if self.gbk_file and self.nucleotide_file:
            os.remove(self.nucleotide_file)


    def set_completeness(self, value):
        self._original_complet = value if value < MAX_COMPLETE else MAX_COMPLETE

    def get_completeness(self):
        return self._original_complet

    def get_redundancy(self):
        return self.original_redundancy

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

    def get_nucleotide_file(self):
        if self.nucleotide_file:
            return nucleotide_file
        if self.gbk_file :
            temp_fna = tempfile.NamedTemporaryFile().name
            if self.gbk_file.endswith(".gz"):
                handle = gzip.open(self.gbk_file, 'rt', encoding='utf-8')
            else :
                handle = open(self.gbk_file, 'r')
            seqs = SeqIO.parse(handle, "genbank")
            with open(temp_fna, "w") as handle2:
                SeqIO.write(seqs,handle2,"fasta")
            handle.close()
            self.nucleotide_file = temp_fna
            return  temp_fna
        raise CantNucleotideError("You need an nucleotide fasta if you want to use a gff")


    def get_amino_acids(self):
        if not self.amino_acids:
            if not self.amino_acid_file:
                if not self.gff_file:
                    raise CantAminoAcidsError("You need either a gbk, gff or a amino-acid fasta if you want to use amino-acids")
                if not self.nucleotide_file:
                    raise CantAminoAcidsError("You need an nucleotide fasta if you want to use a gff to get amino-acids")
                self.amino_acids = { feat.get_id() : feat.get_amino_acids() for feat in self.get_gff() if feat.feature == "CDS"}
            else :
                self.amino_acids = {s.id : s.seq for s in SeqIO.parse(self.amino_acid_file, "fasta")}
            self.amino_acids = {k : v for k, v in self.amino_acids.items() if "*" not in v or v.endswith("*")}
            self.amino_acids = {k : s if s.endswith("*") else (s + "*") for k,s in self.amino_acids.items()}

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
                 '_original_complet' : self._original_complet,
                 'original_redundancy' : self.original_redundancy,
                 'new_completness' : self.new_completness,
        }
    def bin_size(self):
        return sum([len(c) for c in self.get_contigs().values()])

    def overlap(self, target):
        return self.gene_clusters.intersection(target.gene_clusters)

    def estimate_nb_gene_clusters(self):
        assert self.new_completness != None, "new_completness not computed, please do"
        return 100*len(self.gene_clusters)/self.new_completness

    def gene_cluster_ani(self, other_bin):
        tocompare = self.gene_clusters.intersection(other_bin.gene_clusters)
        tocompare = [c for c in tocompare if len(c.get_genes(self)) > 0 and len(c.get_genes(other_bin)) > 0 and len(c) > 2]
        out_pairs = []
        for c in tocompare:
            simis = c.within_cluster_mutations()
            genes1 = c.get_genes(self)
            genes2 = c.get_genes(other_bin)
            pairs = {frozenset({g1,g2}) for g1 in genes1 for g2 in genes2 if g1 != g2}
            pairs = [tuple(p) for p in pairs]
            simis = {p : simis.get(p, simis.get((p[1],p[0]))) for p in pairs if p in simis or (p[1],p[0]) in simis}
            if len(simis) >1:
                for i in range(min(len(genes1), len(genes2))):
                    if len(simis) >1 :
                        best_pair = min(simis.items(), key = lambda ite: (ite[1]['synonymous_muts']+ite[1]['non_synonymous_muts'])/ite[1]['counted_bases'])
                        out_pairs.append(best_pair[1])
                        for p in list(simis):
                            if best_pair[0][0] in p or best_pair[0][1] in p:
                                del simis[p]
            if len(simis) == 1 :
                out_pairs.append(list(simis.items())[0][1])
        out_data = {
        'nb_gene_pairs'                  : len(out_pairs),
        'total_bases'                    : sum([k['counted_bases'] for k in out_pairs]),
        'total_non_synonymous_mutations' : sum([k['non_synonymous_muts'] for k in out_pairs]),
        'total_synonymous_mutations'     : sum([k['synonymous_muts'] for k in out_pairs])
         }
        if out_data['total_bases'] != 0:
            out_data['ani'] = 1-(out_data['total_synonymous_mutations']+out_data['total_non_synonymous_mutations'])/out_data['total_bases']
            out_data['nonsyn_ani'] = (out_data['total_non_synonymous_mutations']/out_data['total_bases'])/(1-out_data['ani'])
            out_data['syn_cont'] = (out_data['total_synonymous_mutations']/out_data['total_bases'])/(1-out_data['ani'])
        else :
            out_data['ani']= None
            out_data['nonsyn_cont']= None
            out_data['syn_cont'] = None
        return out_data
