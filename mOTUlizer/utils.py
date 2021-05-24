import random
import string
import re

def random_name(stringLength=8):
    """Generate a random string of fixed length """
    letters = string.ascii_uppercase
    return ''.join(random.choice(letters) for i in range(stringLength))


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
