import random
import string
import re

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
