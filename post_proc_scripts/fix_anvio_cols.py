import json
import os

os.system("anvi-export-state -p profile.db -o to_fix.json -s default")
os.system("anvi-export-collection -p profile.db -C default -O to_fix_coll")

default_col = '#dcdcdc'
with open("to_fix.json") as handle:
    anvio_profile = json.load(handle)

with open("to_fix_coll.txt") as handle:
    anvio_coll = {l.split()[0] : l.split()[1] for l in handle}

with open("class2col.csv") as handle:
    class2col = {l.split(",")[0].split(";")[1] : l[:-1].split(",")[1] for l in handle}

with open("phylum2col.csv") as handle:
    phylum2col = {l.split(",")[0] : l[:-1].split(",")[1] for l in handle}

for k, v in anvio_profile['categorical_data_colors']['class'].items():
    if k in class2col:
        anvio_profile['categorical_data_colors']['class'][k] = class2col[k]
    else :
        anvio_profile['categorical_data_colors']['class'][k] = default_col

for k, v in anvio_profile['categorical_data_colors']['phylum'].items():
    if k in phylum2col:
        anvio_profile['categorical_data_colors']['phylum'][k] = phylum2col[k]
    else :
        anvio_profile['categorical_data_colors']['phylum'][k] = default_col

with open("collection_info.tsv", "w") as handle:
    handle.writelines([v+ "\tgtdb\t" + phylum2col[v] + "\n" for v in set(anvio_coll.values())])

with open("fixed.json", "w") as handle:
    json.dump(anvio_profile, handle)

os.system("anvi-import-state --state fixed.json  --name default -p profile.db")
os.system("anvi-import-collection --bins-info collection_info.tsv  -C default -p profile.db to_fix_coll.txt")
