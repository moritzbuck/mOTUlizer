import random
import string
import re
import gzip
from enum import Enum, auto
from Bio.Seq import Seq
from Bio import SeqIO
from mOTUlizer.errors import *
import json
import sys
from mOTUlizer import get_quiet

class FASTAtypes(Enum):
    FEATURES = auto()
    CONTIGS = auto()
    RAW = auto()


def message(text, outp = sys.stdout):
    if get_quiet():
        print(text, file = outp, flush = True)


def random_name(stringLength=8):
    """Generate a random string of fixed length """
    letters = string.ascii_uppercase
    return ''.join(random.choice(letters) for i in range(stringLength))

def compute_theoretical_diff_ratio():
    from Bio.Seq import Seq
    alphabet = {'A', 'T', 'C', 'G'}
    all_codons = [Seq(a + b + c) for a in alphabet for b in alphabet for c in alphabet]
    possible_muts = lambda codon: [codon[:i] + base + codon[i+1:] for i in range(3) for base in alphabet if base != codon[i]]
    syn = 0
    non_syn = 0
    for codon in all_codons:
        aa = str(codon.translate())
        for p in possible_muts(codon):
            if str(p.translate()) == aa:
                syn += 1
            else :
                non_syn += 1

def parse_checkm(file):
    # Should be able to read old checkm files, new checkm files  (e.g. checkm lineage_wf infolder datfolder > checkm_file)
    # or any TAB-serparated files with a column called 'Bin Id', 'Completeness', and 'Contamination'
    with open(file) as handle:
        all_lines = [l.strip() for l in  handle.readlines() if " INFO:" not in l]

    all_lines = [re.sub(r"  +","\t", a).split("\t") for a in all_lines]

    header_lines = [i for i,l in enumerate(all_lines) if 'Bin Id' in l and 'Completeness' in l and 'Contamination' in l]
    assert len(header_lines) >= 1, "your completness files is badly formed, it should be TAB-separated (multispaced...) and needs a header line with 'Bin Id', 'Completeness', and 'Contamination' in it"

    header_lines = header_lines[0]
    header_line = all_lines[header_lines]

    lines = [l for i,l in enumerate(all_lines) if i != header_lines and len(l) == len(header_line) and 'Completeness' not in l]

    assert len(lines) > 0, "there are no lines in the file, you sure the serparators are 'OK'"

    lines = [{a : b if a not in ['Completeness', 'Contamination'] else float(b) for a,b in zip(header_line,l) }for l in lines]
    return {l['Bin Id'] : l for l in lines}

def parse_fasta(file, gzipped = False, type = FASTAtypes.RAW):
    if not isinstance(type, FASTAtypes) or type not in FASTAtypes.__members__.values():
        raise CantParseError("the type of FASTA file you want to parse isn't known to me")
    records = []
    id = None
    seq = ""
    if gzipped :
        handle = gzip.open(file, mode = "rt", encoding = "utf-8")
    else :
        handle = open(file)
    id = handle.readline()[1:-1]
    for l in handle:
        if l[0] == ">" :
            records += [(id, seq)]
            id = l[1:].rstrip()
            seq = ""
        else :
            seq += l.rstrip()
    records += [(id, seq)]
    handle.close()

    if type == FASTAtypes.CONTIGS:
        return [{ 'contig_name' : r[0].split(" ")[0] , 'sequence' : r[1], 'annotations' : {} if " " not in r[0] else {"fasta_description" : "".join(r[0].split(" ")[1])} } for r in records]
    elif type == FASTAtypes.FEATURES:
         return [{'name' : r[0].split(" ")[0], 'amino_acids' : r[1], 'annotations' : {} if " " not in r[0] else {"fasta_descritption" : "".join(r[0].split(" ")[1])} } for r in records]
    return records

def parse_genbank(file, gzipped = False):
    if gzipped :
        handle = gzip.open(file, mode = "rt", encoding = "utf-8")
    else :
        handle = open(file)

    records = list(SeqIO.parse(handle, "genbank"))
    handle.close()

    def clean_record_annotation(annot):
        if "references" in annot:
            if len(annot['references']) > 0:
                annot['references'] = [str(ref) for ref in annot['references']]
        return annot

    enumerators = dict()
    def make_uniq_name(contig, ftype):
        if ftype not in enumerators:
            enumerators[ftype] = 0
        enumerators[ftype] += 1
        return f"{contig}_{ftype}_{enumerators[ftype]}"

    contigs = [{ 'contig_name' : r.id ,
                 'sequence' : str(r.seq),
                 'annotations' : json.dumps({"fasta_description" : r.description, "gbk_annotation" : clean_record_annotation(r.annotations)}).replace('"', '""')
                } for r in records]

    c2seq = {c['contig_name'] : Seq(c['sequence']) for c in contigs}

    feats = []
    for s in records:
        for f in s.features:
            if 'translation' in f.qualifiers:
                aas = f.qualifiers['translation'][0]
                del f.qualifiers['translation']
            else :
                aas = ""

            for k,v in list(f.qualifiers.items()):
                f.qualifiers[k] = v[0]

            if f.type not in {'gene', 'source'}:
                feat = {
                "contig_name" : s.id,
                "source" : f"file={file}",
                "feature" : f.type,
                "start" : None if f.location_operator == "join" else int(f.location.start)+1,
                "end" : None if f.location_operator == "join" else  int(f.location.end),
                "score" : "",
                "strand" : "+" if f.location.strand == 1 else "-",
                "frame" : "",
                "name" : make_uniq_name(s.id, f.type) if 'locus_tag' not in f.qualifiers else f.qualifiers['locus_tag'],
                "annotations" : json.dumps(f.qualifiers).replace('"', '""')
                }
                if feat['start']:
                    nucs = c2seq[feat['contig_name']][(feat['start']-1):feat['end']]
                    if feat['strand'] == "-":
                        nucs = nucs.reverse_complement()
                else :
                    nucs = ""
                feat['nucleotides'] = str(nucs)
                if  feat['feature'] == "CDS" and feat['start']:
                    if not aas:
                        aas = nucs.translate()
                    feat['amino_acids'] = str(aas)
                else :
                    feat['amino_acids'] = ""
                feats += [feat]

    handle.close()

    return (contigs, feats)


def parse_gff(file, gzipped = False):
    if gzipped :
        handle = gzip.open(file, mode = "rt", encoding = "utf-8")
    else :
        handle = open(file)

    feats = []
    lines = []
    for l in handle:
        if l.startswith("##FASTA") or l == "\n":
            break
        elif not l.startswith("#"):
            line = l.strip()
            line = line.split("\t")
            attributes = {l.split('=')[0] : l.split('=')[1]for l in line[8].split(";") if l != ""}
            name = attributes.get("ID", "")
            if name in attributes:
                del attributes['name']
            feat = {
            "contig_name" : line[0],
            "source" : f"{line[1]}:file={file}",
            "feature" : line[2],
            "start" : int(line[3]),
            "end" : int(line[4]),
            "score" : line[5],
            "strand" : line[6],
            "frame" : line[7],
            "name" : name,
            "annotations" : json.dumps(attributes).replace('"', '""')
            }
            feats += [feat]
    handle.close()
    contigs = {v['contig_name'] for v in feats}
    from mOTUlizer.db.SeqDb import SeqDb
    if SeqDb.seq_db:
        contigs = {k : Seq(v) for k,v in SeqDb.seq_db.get_contigs(contigs).items()}
    for feat in feats:
        nucs = contigs[feat['contig_name']][(feat['start']-1):feat['end']]
        if feat['strand'] == "-":
            nucs = nucs.reverse_complement()
        aas = nucs.translate()
        feat['nucleotides'] = str(nucs)
        if  feat['feature'] == "CDS":
            feat['amino_acids'] = str(aas)
        else :
            feat['amino_acids'] = ""
    return feats
