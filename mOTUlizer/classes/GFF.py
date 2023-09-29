from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from mOTUlizer.errors import CantAminoAcidsError
from math import floor
from mOTUlizer.db.SeqDb import SeqDb

class GFF():
    def __len__(self):
        return len(self.entries)

    def __getitem__(self, i):
        if type(i) == int and i < len(self):
            return self.entries[i]
        elif i in self.cdss:
            return self.cdss[i]
        else :
            raise KeyError(str(i) + "is not a valid entry")

    def __iter__(self):
       return GFFIterator(self)

    def __init__(self, genome , gff_file : str = "", db = None):
        if not db:
            self.db = SeqDb.get_global_db()
            if not SeqDb.seq_db:
                raise DataBaseNotInitialisedError("The database has not been initialised")
        else : 
            self.db = db

        self.entries = []
        self.genome = genome

        with open(gff_file) as handle:
            for l in handle:
                if l.startswith("##FASTA") or l == "\n":
                    break
                elif not l.startswith("#"):
                    self.entries.append(GFFentry(l.strip(), self))
        self.cdss = {e.attributes["ID"] : e for e in self.entries if e.feature == "CDS"}

    def make_fasta(self, file, type = "nucl"):
        if type == "aas":
            records = [SeqRecord(id = s.attributes['ID'], description = "", seq = s.get_amino_acids()) for s in self if s.feature == "CDS"]
        elif type == "nucl":
            records = [SeqRecord(id = s.attributes['ID'], description = "", seq = s.get_gene()) for s in self]
        else :
            raise ValueError("type should be aas or nucl for amino-acids or nucleotides respectively")
        with open(file, "w") as handle:
            SeqIO.write(records, handle, "fasta")

class GFFIterator:
    def __init__(self, gff):
        self._gff = gff
        self.index = -1

    def __next__(self):
        if self.index < len(self._gff)-1:
            self.index += 1
            return self._gff[self.index]
        raise StopIteration

class GFFentry():
    def __init__(self, line, gff):
        self.gff = gff
        line = line.split("\t")
        self.seqname = line[0]
        self.source = line[1]
        self.feature = line[2]
        self.start = int(line[3])
        self.end = int(line[4])
        try :
            self.score = float(line[5])
        except :
            self.score = None
        self.strand = line[6]
        try :
            self.frame = int(line[7])
        except :
            self.frame = None

        self.attributes = {l.split('=')[0] : l.split('=')[1]for l in line[8].split(";") if l != ""}

    def get_gene(self):
        if self.strand == "-":
            seq = self.gff.genome.get_contig(self.seqname)[(self.start-1):(self.end)]
            seq = seq.reverse_complement()
            if seq[-3:].translate() != "*":
                seq += self.gff.genome.get_contig(self.seqname)[(self.start-4):(self.start-1)].reverse_complement()
        else :
            seq = self.gff.genome.get_contig(self.seqname)[(self.start-1):(self.end)]
            if seq[-3:].translate() != "*":
                seq += self.gff.genome.get_contig(self.seqname)[(self.end):(self.end+4)]
#        assert len(seq) % 3 == 0
        return seq

    def get_id(self):
        if "ID" not in self.attributes:
            raise ValueError("ID is not in the attributes of the feature you are looking for")
        return self.attributes["ID"]

    def get_amino_acids(self):
        if self.feature != "CDS":
            raise CantAminoAcidsError("the feature you are getting can't be translated")
        gene = self.get_gene()
        gene = gene[:(3*floor(len(gene)/3))]
        return gene.translate()
