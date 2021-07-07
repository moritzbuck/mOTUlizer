import abc
from tqdm import tqdm
import sys

class Parser(metaclass=abc.ABCMeta):

    def __init__(self, **kwargs):
        if 'gene_id2genome' in kwargs and gene_id2genome:
            if type(gene_id2genome) == dict:
                self.gene_id2genome = gene_id2genome
            else :
                print("TODO")
                sys.exit(0)
        else :
            self.gene_id2genome = None

    @abc.abstractmethod
    def convert(self, infile, outfile = None, **kwargs):
        pass

class EmapperParse(Parser):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        try :
            from ete3 import NCBITaxa
        except :
            print("""
            To parse eggNOGs properly I need the ete3 python package,
            you could for example install it with pip as in :
                pip install ete3

            or, if conda is more your cupa tea, you can do :
                conda install -c conda-forge -c etetoolkit ete3=3.1.2

            or whichever other way you like. I don't judge.
            """, file = sys.stderr)
            sys.exit(0)

    def convert(self, infile, count = False):
        from ete3 import NCBITaxa
        tax_db = NCBITaxa()

        with open(infile) as handle:
            header = None
            counter = 0
            while not header and counter < 15:
                head = handle.readline().rstrip().split("\t")
                if head[0] == "#query_name" :
                    header = head
                counter += 1
            if not header :
                print("\nYou sure this file is good? Like, is it the '.emapper.annotations' you got from running eggnoggmapper?\n")
                sys.exit(0)

            idx = [i for i,v in enumerate(header) if v == "eggNOG OGs"][0]
            print("Loading eggNOGs from file.", file = sys.stderr)
            gene_id2eggs = {l.split("\t")[0] : l.split("\t")[idx] for l in tqdm(handle.readlines()) if not l.startswith("#")}


        print("Parsing taxonomies, and simplifying to deepest eggNOG.", file = sys.stderr)
        taxos = {vv.split("@")[1] for v in gene_id2eggs.values() for vv in v.split(",")}
        tax2level = {k : len(v) for k,v in tax_db.get_lineage_translator(list(taxos)).items()}
        lowest_level = lambda x : tax2level.get(int(x[1]),1000)
        gene_id2deepest_egg = {k:  min([vv.split('@') for vv in v.split(",")], key = lowest_level)[0] for k,v in tqdm(gene_id2eggs.items())}

        print("Stratify it to genome", file = sys.stderr)
        if not self.gene_id2genome:
            self.gene_id2genome = {k : "_".join(k.split("_")[:-1]) for k in gene_id2deepest_egg.keys()}

        genome2nog = {k : [] for k in set(self.gene_id2genome.values())}
        for k,v in gene_id2deepest_egg.items():
            genome2nog[self.gene_id2genome[k]] += [v]
        if count:
            return {k : {vv : v.count(vv) for vv in set(v)} for k,v in genome2nog.items()}
        else :
            return {k : list(set(v)) for k,v in genome2nog.items()}



class PPanGGolinParse(Parser):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        try :
            import h5py
            import hdf5plugin
        except :
            print("""
            To parse PPanGGolin's h5 file properly I need the h5py and hdf5plugin
            python package, you could for example install it with pip as in :
                pip install h5py hdf5plugin

            or, if conda is more your cupa tea, you can do :
                conda install -c anaconda -c conda-forge ete3

            or whichever other way you like. I don't judge.
            """, file = sys.stderr)
            sys.exit(0)


    def convert(self, infile, count = False):
        import h5py
        import hdf5plugin

        print("Parse and decode gene-clusters from the h5py-file", file = sys.stderr)
        with h5py.File(infile, "r") as handle:
            gene_id2family = {a.decode() : b.decode() for a,b in tqdm(handle['geneFamilies'])}

        print("Stratify it to genome", file = sys.stderr)
        if not self.gene_id2genome:
            self.gene_id2genome = {k : "_".join(k.split("_")[:-2]) for k in gene_id2family.keys()}

        genome2family = {k : [] for k in set(self.gene_id2genome.values())}
        for k,v in gene_id2family.items():
            genome2family[self.gene_id2genome[k]] += [v]
        if count:
            return {k : {vv : v.count(vv) for vv in set(v)} for k,v in genome2family.items()}
        else :
            return {k : list(set(v)) for k,v in genome2family.items()}


class AnvioParse(Parser):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        try :
            import anvio.dbops as dbops
        except :
            print("""
            To parse anvi'o's pangenome database you will need to have anvio installed
            for now at least, until I move my ass and just dig in the SQLite stuff...
            But if you have a pangenome database, you probably have anvio installed anyhow.
            Otherwise, best check the anvi'o install page :
                    https://merenlab.org/2016/06/26/installation-v2/
            """, file = sys.stderr)
            sys.exit(0)


    def convert(self, infile, count = False):
        class Mock:
            def __init_(self):
                self.__dict__ = {}
        import anvio.dbops as dbops
        args = Mock()
        args.__dict__['pan_db'] = infile
        pan = dbops.PanSuperclass(args)

        gene_cluster_ids = pan.gene_cluster_names

        pan.init_gene_clusters(gene_cluster_ids)
        gene_cluster = set(pan.gene_clusters.keys())

        genomes = {g for k,v in pan.gene_clusters.items() for g in v}
        genome2genecluster = {g : [] for g in genomes}
        for gc, hits in pan.gene_clusters.items():
            for genome, genes in hits.items():
                if len(genes) > 0:
                    if gc not in genome2genecluster[genome]:
                        genome2genecluster[genome].append(gc)
        if count:
            return {k : {vv : v.count(vv) for vv in set(v)} for k,v in genome2genecluster.items()}
        else :
            return {k : list(set(v)) for k,v in genome2genecluster.items()}


class RoaryParse(Parser):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def convert(self, infile, count = False):
        print("Parse the cluster-file", file = sys.stderr)
        with open(infile, "r") as handle:
            gene_id2family = {ll : l.split(':')[0] for l in tqdm(handle) for ll in l.split(":")[1].strip().split()}

        print("Stratify it to genome", file = sys.stderr)
        if not self.gene_id2genome:
            self.gene_id2genome = {k : "_".join(k.split("_")[:-1]) for k in gene_id2family.keys()}

        genome2family = {k : [] for k in set(self.gene_id2genome.values())}
        for k,v in gene_id2family.items():
            genome2family[self.gene_id2genome[k]] += [v]
        if count:
            return {k : {vv : v.count(vv) for vv in set(v)} for k,v in genome2family.items()}
        else :
            return {k : list(set(v)) for k,v in genome2family.items()}

class MmseqsParse(Parser):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def convert(self, infile, count = False):
        print("Parse the cluster-file", file = sys.stderr)
        with open(infile, "r") as handle:
            gene_id2family = {l.strip().split('\t')[1] : l.strip().split('\t')[0] for l in tqdm(handle)}

        print("Stratify it to genome", file = sys.stderr)
        if not self.gene_id2genome:
            self.gene_id2genome = {k : "_".join(k.split("_")[:-1]) for k in gene_id2family.keys()}

        genome2family = {k : [] for k in set(self.gene_id2genome.values())}
        for k,v in gene_id2family.items():
            genome2family[self.gene_id2genome[k]] += [v]
        if count:
            return {k : {vv : v.count(vv) for vv in set(v)} for k,v in genome2family.items()}
        else :
            return {k : list(set(v)) for k,v in genome2family.items()}
