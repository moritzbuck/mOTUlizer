from mOTUlizer.db.SeqDb import SeqDb
from random import choices, choice, uniform, normalvariate, seed
from numpy import mean, log10
from sourmash import MinHash
from subprocess import call
from tqdm import tqdm
import os
import shutil
import json

BASES = {'A','T','C','G'}
CHOICES = { b : BASES.difference(b) for b in BASES}
min_hash_params =  {
    'n' : 0,
    'ksize' : 21,
    'scaled' : 10000
}
class Genome(object):

     def __get__(self, i):
         return self.genome[i]

     def __iter__(self): return iter(self.genome)
     def __getitem__(self, key): return self.genome[key]
     def __len__(self) : return self.length

     def __init__(self, sequence, mutation_rate, var_fract):
         self.length = len(sequence)
         self.mutation_rate = mutation_rate
         self.var_fract = var_fract
         self.seq = sequence
         self._compute_minhash()

     def _new_base(self,position, base):
         if base not in BASES:
             return base
         if position < (self.length*self.var_fract):
             return choice(list(self.seq))
         if uniform(0,1) < self.mutation_rate:
             return choice(list(CHOICES[base]))
         return base

     def mutate(arg):
        pass

     def mutate(self):
         self.seq = "".join([self._new_base(i,s) for i,s in enumerate(self.seq)])

         self._compute_minhash()

     def copy(self):
         return Genome(self.seq, self.mutation_rate, self.var_fract)

     @property
     def minhash(self):
         return self._minhash

     def _compute_minhash(self):
         self._minhash = MinHash(**min_hash_params)
         self.minhash.add_sequence(self.seq, force = True)

     def real_ani(self, other):
         simi = lambda a,b : sum([xy[0] == xy[1] for i,xy in enumerate(zip(a,b)) if i> (self.length*self.var_fract)])/(len(a)*(1-self.var_fract))

         return simi(self.seq,other.seq)

     def minhash_ani(self,other, method = "containment"):
         if method == "containment":
             return 1-self.minhash.containment_ani(other.minhash).dist
         if method == "max_containment":
             return 1-self.minhash.max_containment_ani(other.minhash).dist
         if method == "avg_containment":
             return self.minhash.avg_containment_ani(other.minhash)
         if method == "jaccard":
             return 1-self.minhash.jaccard_ani(other.minhash).dist
     def fastani_ani(self,other):
         with open("genome1.fa", "w") as handle:
             handle.write(f">1\n{self.seq}\n")
         with open("genome2.fa", "w") as handle:
             handle.write(f">1\n{other.seq}\n")
         call("fastANI -t 24 -q genome1.fa -r genome2.fa -o output.txt 2> /dev/null", shell= True)
         with open("output.txt") as handle:
             return [float(l.split()[2]) for l in handle][0]/100

     def pyani(self, other):
         os.makedirs("g", exist_ok = True)
         with open("g/genome1.fa", "w") as handle:
             handle.write(f">1\n{self.seq}\n")
         with open("g/genome2.fa", "w") as handle:
             handle.write(f">1\n{other.seq}\n")
         call("average_nucleotide_identity.py -i g -o pyaniout -m ANIb --nocompress", shell= True)
         with open("pyaniout/ANIb_percentage_identity.tab") as handle:
            tt = [float(l.split()[3-i]) for i,l in enumerate(handle) if i >0]
         shutil.rmtree("pyaniout")
         return sum(tt)/2

SeqDb.init_db("/home/moritz/projects/0064_bis/data/all_genomes.sqlite")
db = SeqDb.get_global_db()

test1 =[]
seed(42)
for j in range(10):
     genome = db.get_random_genome()
 #    genome = MetaBin(j)
     contigs = "".join(genome.get_raw_seqs(clean = False))
     g = Genome(contigs, 0.0, 0.0)
     print(f"===== {genome.name} =====")
     for mut in [1 - l/1000 for l in range(850, 1000, 1)]:
         for var in [0]: #,0.05,0.10,0.20]:
             print(f"mutation rate : {mut}, variable_fraction : {var}")
             out = dict()
             g2 = g.copy()
             g2.mutation_rate = mut
             g2.var_fract = var
             g2.mutate()
             out['genome'] = genome.name
             out['ani'] = g2.real_ani(g)
             out['fastANI'] = g2.fastani_ani(g)
             out['cont'] = g2.minhash_ani(g)
             out['anib'] = g2.pyani(g)
             out['var_fract'] = g2.var_fract
             out['mut_rate'] = g2.mutation_rate
             out['phylum'] = json.loads(db.get_genome_info(out['genome'])['taxonomy'].replace('""','"'))['r207'].split(";c__")[0]
             out['family'] = json.loads(db.get_genome_info(out['genome'])['taxonomy'].replace('""','"'))['r207'].split(";s__")[0]
             test1 += [out]
             pandas.DataFrame.from_records(test1).to_csv("/home/moritz/temp/themievo2.csv", index=None)
