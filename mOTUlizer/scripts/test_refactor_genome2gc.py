
import os, sys

from tqdm import tqdm
import shutil
import traceback
import shutil
import traceback
import pandas
from statistics import mean, median
import json 
sys.path = ['/home/moritz/projects/0039_mOTUlizer/mOTUlizer'] + sys.path


from mOTUlizer.classes.MetaBin import MetaBin
from mOTUlizer.classes.mOTU import mOTU
from mOTUlizer.db.SeqDb import SeqDb
from mOTUlizer.classes.tools.Prokka import Prokka
import mOTUlizer
from mOTUlizer.classes.GeneClusters import compute_GeneClusters
mOTUlizer._quiet_ = False

if os.path.exists(f"test/genome2gc_test.db"):
	os.remove(f"test/genome2gc_test.db")
SeqDb.init_db(f"test/genome2gc_test.db")
db = SeqDb.get_global_db()



gene_clusters = [
	{'name' :"GC1"},
	{'name' :"GC2"},
	{'name' :"GC3"},
	{'name' :"GC6"}
]

gene_clusters = [
	{'name' :"GC1"},
	{'name' :"GC2"},
	{'name' :"GC3"},
	{'name' :"GC4"},
	{'name' :"GC5"}

]



genome = MetaBin(name = "test1", gene_clusters = gene_clusters)
genome2 = MetaBin(name = "test2", gene_clusters = gene_clusters)
tt = mOTU(name = "test_mOTU", genomes = [genome, genome2])
tt.compute_core()

AGs = { "GB_GCA_000006155.2" : {
    "name": "GB_GCA_000006155.2",
    "gbk_file": "/home/moritz/data/gtdb//files/GCA/000/006/155/GCA_000006155.2/GCA_000006155.2_genomic.gbff.gz",
    "genome_completeness": 93.12,
    "genome_redundancy": 0.0,
    "taxonomy": {
      "r207": "d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae_G;g__Bacillus_A;s__Bacillus_A anthracis"
    	}
    },  "GB_GCA_000160995.1": {
    "name": "GB_GCA_000160995.1",
    "gbk_file": "/home/moritz/data/gtdb//files/GCA/000/160/995/GCA_000160995.1/GCA_000160995.1_genomic.gbff.gz",
    "genome_completeness": 99.28,
    "genome_redundancy": 1.36,
    "taxonomy": {
      "r207": "d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae_G;g__Bacillus_A;s__Bacillus_A luti"
    }
  },  "GB_GCA_000161075.1": {
    "name": "GB_GCA_000161075.1",
    "gbk_file": "/home/moritz/data/gtdb//files/GCA/000/161/075/GCA_000161075.1/GCA_000161075.1_genomic.gbff.gz",
    "genome_completeness": 98.61,
    "genome_redundancy": 0.19,
    "taxonomy": {
      "r207": "d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae_G;g__Bacillus_A;s__Bacillus_A paranthracis"
    }
  },  "GB_GCA_000161475.1": {
    "name": "GB_GCA_000161475.1",
    "gbk_file": "/home/moritz/data/gtdb//files/GCA/000/161/475/GCA_000161475.1/GCA_000161475.1_genomic.gbff.gz",
    "genome_completeness": 98.88,
    "genome_redundancy": 0.41,
    "taxonomy": {
      "r207": "d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae_G;g__Bacillus_A;s__Bacillus_A tropicus"
    }
  },  "GB_GCA_000161535.1": {
    "name": "GB_GCA_000161535.1",
    "gbk_file": "/home/moritz/data/gtdb//files/GCA/000/161/535/GCA_000161535.1/GCA_000161535.1_genomic.gbff.gz",
    "genome_completeness": 98.68,
    "genome_redundancy": 0.15,
    "taxonomy": {
      "r207": "d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae_G;g__Bacillus_A;s__Bacillus_A thuringiensis_S"
    }
  }
  }





genomes = [MetaBin(**gd, temp_dir = "/home/moritz/temp/") for gd in AGs.values()]
motu2 = mOTU('test_nucs', genomes)
bla = compute_GeneClusters(motu2, "tes_clust")
motu2.compute_core()
motu2.get_stats()



import os, sys

from tqdm import tqdm
import shutil
import traceback
import shutil
import traceback
import pandas
from statistics import mean, median
import json 
sys.path = ['/home/moritz/projects/0039_mOTUlizer/mOTUlizer'] + sys.path


from mOTUlizer.classes.MetaBin import MetaBin
from mOTUlizer.classes.mOTU import mOTU
from mOTUlizer.db.SeqDb import SeqDb
from mOTUlizer.classes.tools.Prokka import Prokka
import mOTUlizer
from mOTUlizer.classes.GeneClusters import compute_GeneClusters
from random import shuffle, choice, choices

mOTUlizer._quiet_ = False

SeqDb.init_db(f"test/genome2gc_test.db")
db = SeqDb.get_global_db()

motu2 = mOTU('test_nucs')

motu2.roc_values()