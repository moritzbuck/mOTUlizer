import sqlite3
import os
from mOTUlizer import _quiet_
from mOTUlizer.utils import parse_fasta
from enum import Enum, auto
import json
from sourmash.signature import save_signatures, load_signatures
from sourmash.signature import SourmashSignature
from sourmash import MinHash
import re
from tqdm import tqdm

seq_types = {'amino_acids', 'nucleotides'}

class SeqDb:
    seq_db = None

    @classmethod
    def init_db(cls, path = None):
        cls.seq_db = SeqDb(path)

    @classmethod
    def get_global_db(cls):
        if not cls.seq_db:
            raise DataBaseNotInitialisedError("The database has not been initialised")
        return cls.seq_db

    def __repr__(self) :
        return f"< SQLite3-based Sequence DB located at {self.path} >"

    def __init__(self, path):
        if not path:
            self.path = tempfile.NamedTemporaryFile().name
            self.temp = True
        else :
            self.path = path
            self.temp = False

        if not os.path.exists(path):
            if not _quiet_:
                print("creating db")
            self._connection = sqlite3.connect(self.path)
            self.create_db()
        else :
            if not _quiet_:
                print("loading db")
            self._connection = sqlite3.connect(self.path)
            self.check_db()
        self._open_cursor = None

    def __del__(self):
        self._connection.close()
        
        if self.temp:
            os.remove(self.path)

    def create_db(self):
        cursor = self._connection.cursor()

        cursor.execute('''
                    CREATE TABLE "genomes" (
                    	"genome_name"	TEXT,
                    	"source"	TEXT, "annotations" TEXT, "taxonomy" TEXT, "completeness" REAL, "redundancy" REAL, "signature" TEXT,
                    	PRIMARY KEY("genome_name")
                    );
                  ''')

        cursor.execute('''
                    CREATE TABLE "anis" (
                    	"query_name"	TEXT,
                    	"subject_name"	TEXT, "ani"  REAL, "query_chunks" INTEGER, "reference_chunks" INTEGER,
                    	PRIMARY KEY("query_name", "subject_name"),
                    	FOREIGN KEY("query_name") REFERENCES genomes("genome_name"),
                    	FOREIGN KEY("subject_name") REFERENCES genomes("genome_name")
                    );
                  ''')


        cursor.execute('''
                    CREATE TABLE "contigs"  (
                    	"contig_name"	TEXT,
                    	"sequence"	TEXT,
                    	"genome_name"	TEXT, "annotations" TEXT,
                    	PRIMARY KEY("contig_name"),
                    	FOREIGN KEY("genome_name") REFERENCES genomes("genome_name")
                    );
                  ''')

        cursor.execute('''
                    CREATE TABLE "gene_clusters"  (
                    	"gc_name"	TEXT,
                    	"annotations"	TEXT, "alignment" TEXT, "tree" TEXT, "representative" TEXT,
                    	PRIMARY KEY("gc_name"), FOREIGN KEY("representative") REFERENCES features("feature_id")
                    );
                  ''')

        cursor.execute('''
                    CREATE TABLE "cores"  (
                    	"gc_name"	TEXT, "motu_name" TEXT, "loglikelihood" TEXT, FOREIGN KEY("gc_name") REFERENCES gene_clusters("gc_name"), FOREIGN KEY("motu_name") REFERENCES mOTUs("motu_name"),
                    	PRIMARY KEY("gc_name", "motu_name")
                    );
                  ''')

        cursor.execute('''
                    CREATE TABLE "gc2feature"  (
                    	"gc_name"	TEXT,
                    	"feature_id"	TEXT UNIQ,
                    	FOREIGN KEY("gc_name") REFERENCES gene_clusters("gc_name")
                    	FOREIGN KEY("feature_id") REFERENCES features("feature_id")
                    );
                  ''')


        cursor.execute('''
                    CREATE TABLE "genome2gcs"  (
                        "gc_name"   TEXT,
                        "genome_name"    TEXT,
                        FOREIGN KEY("gc_name") REFERENCES gene_clusters("gc_name")
                        FOREIGN KEY("genome_name") REFERENCES genomes("genome_name")
                    );
                  ''')

        cursor.execute('''
                    CREATE INDEX "gc2genomes"
                    ON genome2gcs("gc_name");
                  ''')

        cursor.execute('''
                    CREATE INDEX "genomes2gcs"
                    ON genome2gcs("genome_name");
                  ''')


        cursor.execute('''
                    CREATE INDEX "gc_name_index"
                    ON gc2feature("gc_name");
                  ''')


        cursor.execute('''
                    CREATE INDEX "feature2gc"
                    ON gc2feature("feature_id");
                  ''')



        cursor.execute('''
                    CREATE TABLE "mOTUs"  (
                    	"motu_name"	TEXT,
                    	"taxonomy"	TEXT,
                    	PRIMARY KEY("motu_name")
                    );
                  ''')

        cursor.execute('''
                    CREATE TABLE "genome2motu"  (
                    	"motu_name"	TEXT,
                    	"genome_name"	TEXT,
                    	FOREIGN KEY("genome_name") REFERENCES genomes("genome_name")
                    	FOREIGN KEY("motu_name") REFERENCES mOTUs("motu_name")
                    );
                  ''')

        cursor.execute('''
                    CREATE INDEX "motu2genome"
                    ON genome2motu("genome_name");
                  ''')

        cursor.execute('''
                    CREATE INDEX "genome_to_contigs"
                    ON contigs("genome_name");
                  ''')

        cursor.execute('''
                    CREATE TABLE "features" (
                        "contig_name"   TEXT,
                    	"feature_id"	INTEGER,
                    	"name"	TEXT,
                    	"source"	TEXT,
                    	"feature"	TEXT,
                    	"start"	INTEGER,
                    	"end"	INTEGER,
                    	"score"	TEXT,
                    	"strand"	TEXT,
                    	"frame"	TEXT,
                    	"amino_acids"	TEXT,
                    	"nucleotides"	TEXT, 
                        "annotations" TEXT,
                    	PRIMARY KEY("feature_id" AUTOINCREMENT), FOREIGN KEY("contig_name") REFERENCES contigs("contig_name")

                    );
                  ''')

        cursor.execute('''
                    CREATE INDEX "contig_to_features"
                    ON features("contig_name");
                  ''')

        cursor.execute('''
                    CREATE INDEX "name_to_features"
                    ON features("name");
                  ''')


        self._connection.commit()



    @property
    def open_cursor(self):
        if not self._open_cursor:
            self._operation_counter = 1
            self._open_cursor = self._connection.cursor()
        else :
            self._operation_counter += 1
        return self._open_cursor

    def commit(self):
        self._connection.commit()
        self._open_cursor = None


    def bunched_query(func):
        def wrapper(*args, **kwargs):
            func(*args, **kwargs)
            if not 'bunched' in kwargs or not kwargs['bunched']:
                args[0].commit()
        return wrapper

    @bunched_query
    def add_genome(self, genome_name : str, source : str = '', contigs : list = [], annotations : dict = None, taxonomy = None, completeness = None, redundancy = None, signature = None, bunched = False):
        # contigs is list of dicts dicts with keys  'contig_name' 'sequence' and optinally 'annotations'
        # features is list of dicts with keys 'sequence' and optionaly 'annotations'


        if self.has_genome(genome_name):
            raise DataBaseBadValuesError(f"Genome {genome_name} is already in the database")

        if annotations:
            annotations = json.dumps(annotations).replace('"', '""')

        if taxonomy:
            taxonomy = json.dumps(taxonomy).replace('"', '""')

        if signature:
            if type(signature) == dict:
                mh = MinHash(**signature)
                for c in contigs:
                    mh.add_sequence(c['sequence'], force = True)
                signature = save_signatures([SourmashSignature(mh)])

        self.open_cursor.execute(f"""
        INSERT INTO genomes (genome_name,source, annotations, taxonomy, completeness, redundancy, signature)
        VALUES (?, ?, ?, ?, ?, ?, ?);
        """, (genome_name, source, annotations, taxonomy, completeness, redundancy, signature))
        if contigs:
            self.add_contigs(genome_name, contigs, bunched = bunched)

    @bunched_query
    def add_contigs(self, genome_name : str, contigs : list , bunched = False):
        contigs = [(v['contig_name'], v['sequence'], genome_name, v.get('annotations', {})) for v in contigs]
        values = " , ".join([ f'("{a}", "{b}", "{c}", "{d}")' for a,b,c,d in contigs ])
        self.open_cursor.execute(f"""
        INSERT INTO contigs (contig_name,sequence,genome_name, annotations)
        VALUES {values};
        """)

    @bunched_query
    def add_features(self, genome_name : str, features : dict, source : str = "" , bunched = False):
        if features :
            fields = ["contig_name" , "source" , "feature" ,"start","end","score" ,"strand" , "frame", "name", "annotations", "amino_acids", "nucleotides"]
            values = " , ".join([ "(" + ",".join(['"{}"'.format(v.get(f,"")) for f in fields]) + ")" for v in features ])
            self.open_cursor.execute(f"""
            INSERT INTO features ({",".join(fields)})
            VALUES {values};
            """)

    @bunched_query
    def add_anis(self, ani_dict : dict, bunched = False):
        fields = ["query_name", "subject_name" , "ani" , "query_chunks" , "reference_chunks"]

        values = " , ".join([ f'( "{k[0]}", "{k[1]}" , ' + ",".join(['"{}"'.format(v.get(f,"")) for f in fields[2:]]) + ")" for k,v in ani_dict.items() ])
        self.open_cursor.execute(f"""
        INSERT INTO anis ({",".join(fields)})
        VALUES {values};
        """)

    @bunched_query
    def add_mOTU(self, motu_name : str, genomes : list, bunched = False):
        self.open_cursor.execute(f"""
        INSERT INTO mOTUs (motu_name)
        VALUES (?);
        """, (motu_name,))
        self.open_cursor.executemany(f"""
        INSERT INTO genome2motu (motu_name, genome_name)
        VALUES (?,?);
        """, [(motu_name, g) for g in genomes])


    @bunched_query
    def create_gcs(self, gcs : list, bunched = False):
        """
        { 
            'name' : str,
            'feature' : [ids],
            'genomes' : [ids],
            'representative' : id,
            'annotations' : dict
        }
        """
        data = [ (gc['name'], gc.get('representative'),json.dumps(gc['annotations']).replace('"', '""') if gc.get('annotations') else "{}") for gc in gcs]
        self.open_cursor.executemany(f"""
        INSERT INTO gene_clusters (gc_name, representative, annotations)
        VALUES (?, ?, ?);
        """, data)



    @bunched_query
    def add_gene_cluster(self, name :str , features : list = None , genomes : list = None, representative : int = None, annotations = None, bunched = False):
        if annotations:
            annotations = json.dumps(annotations).replace('"', '""')
        else :
            annotations = "{}"
        self.open_cursor.execute(f"""
        INSERT INTO gene_clusters (gc_name, representative, annotations)
        VALUES (?, ?, ?);
        """, (name,representative, annotations))
        
        if genomes:
            self.add_genome_tuple2gc(name, [(name, g) for g in genomes], bunched)

        if features:
            self.add_feature_tuple2gc(name, [(name, f) for f in features], bunched)


    @bunched_query
    def add_feature_tuple2gc(self, gc_feature_tuples : list, bunched = False):
        self.open_cursor.executemany(f"""
        INSERT INTO gc2feature (gc_name, feature_id)
        VALUES (?, ?);
        """, [f for f in gc_feature_tuples])

    @bunched_query
    def add_genome_tuple2gc(self, gc_genome_tuples : list, bunched = False):
        self.open_cursor.executemany(f"""
        INSERT INTO genome2gcs (gc_name, genome_name)
        VALUES (?, ?);
        """, [f for f in gc_genome_tuples])


    @bunched_query
    def update_anis(self, ani_dict : dict, bunched = False):
        genomes = {kk for k in ani_dict for kk in k}
        values = " , ".join([ f'"{kk}"' for kk in genomes ])
        self.open_cursor.execute(f"""
        DELETE FROM anis
        WHERE "query_name" IN ({values});
        """)
        self.open_cursor.execute(f"""
        DELETE FROM anis
        WHERE "subject_name" IN ({values});
        """)
        self.add_anis(ani_dict)

    @bunched_query
    def update_taxonomy(self, tax_dict : dict, bunched = False):
        current_tax = self.get_taxonomy(tax_dict.keys())
        for k,v in tax_dict.items():
            current_tax[k].update(v)
        for c in current_tax:
            self.set_value("genome", c, "taxonomy", json.dumps(current_tax[c]).replace('"', '""'))

    @bunched_query
    def set_value(self, database, key_name, key, column, value):
        self.open_cursor.execute(f'UPDATE {database} SET {column}={value} WHERE {key_name}="{key}";')


    def delete_genome(self, genome_name : str):
        pass

    def has_genome(self, genome_name):
        self.open_cursor.execute(f"""
        SELECT genome_name FROM genomes WHERE genome_name=? LIMIT 1;
        """, (genome_name,))
        return len(self.open_cursor.fetchall()) > 0

    def has_gc(self, gc_name):
        self.open_cursor.execute(f"""
        SELECT gc_name FROM gene_clusters WHERE gc_name=? LIMIT 1;
        """, (gc_name,))
        return len(self.open_cursor.fetchall()) > 0

    def has_features(self, genome_name):
        self.open_cursor.execute(f"""
        SELECT feature_id FROM features INNER JOIN contigs on contigs.contig_name = features.contig_name WHERE genome_name = ? AND feature = "CDS" LIMIT 1 ;
        """, (genome_name,))
        return len(self.open_cursor.fetchall()) > 0

    def has_any_gcs(self):
        self.open_cursor.execute(f"""
        SELECT gc_name FROM gene_clusters LIMIT 1 ;
        """)
        return len(self.open_cursor.fetchall()) > 0




    def has_mOTU(self, motu_name):
        self.open_cursor.execute(f"""
        SELECT motu_name FROM mOTUs WHERE motu_name=? LIMIT 1;
        """, (motu_name,))
        return len(self.open_cursor.fetchall()) > 0



    def get_anis(self, genomes):
        genome_names = " ( " + ",".join([f' "{g1.name}" ' for g1 in genomes]) + " ) "

        self.open_cursor.execute(f"""
        SELECT query_name, subject_name, ani,query_chunks, reference_chunks  FROM anis
        WHERE query_name IN {genome_names} AND subject_name in {genome_names} ;
        """)
        anis = self.open_cursor.fetchall()
        return { (v[0],v[1]) : {'ani' : v[2], 'query_chunks' : v[3], 'reference_chunks' : v[4]} for v in anis}

    def get_mOTU(self, motu):
        self.open_cursor.execute(f"""
        SELECT genome_name FROM genome2motu
        WHERE motu_name=? ;
        """, (motu,))
        genomes = self.open_cursor.fetchall()
        return [g[0] for g in genomes]

    def get_genomesFromFeat(self, feat):
        self.open_cursor.execute(f"""
        SELECT genome_name FROM contigs
        INNER JOIN features ON contigs.contig_name = features.contig_name
        WHERE feature_id = ? LIMIT 1;
        """, (feat,))
        genomes = self.open_cursor.fetchall()
        return [g[0] for g in genomes][0]

    def list_mOTUs(self):
        self.open_cursor.execute(f"""
        SELECT motu_name FROM mOTUs ;
        """)
        motus = [m[0] for m in self.open_cursor.fetchall()]
        return {m : self.get_mOTU(m) for m in motus}

    def list_GCs(self):
        self.open_cursor.execute(f"""
        SELECT gc_name FROM gene_clusters ;
        """)
        gcs = [m[0] for m in self.open_cursor.fetchall()]
        return gcs

    def get_gc(self, name):
        feats = self.get_features(name)
        fields = ["annotations", "representative"]
        self.open_cursor.execute(f"""
        SELECT {' , '.join(fields)} FROM gene_clusters WHERE gc_name=?;
        """, (name,))
        data = self.open_cursor.fetchall()
        data = {k : v for k,v in zip(fields, data[0])}
        data['annotations'] = json.loads(data['annotations'].replace('""', '"'))
        data['gc_name'] = name 
        data['features'] =  feats 
        return data 


    def get_features(self, gene_cluster):
        self.open_cursor.execute(f"""
        SELECT feature_id FROM gc2feature WHERE gc_name=?;
        """, (gene_cluster,))
        features = [m[0] for m in self.open_cursor.fetchall()]
        return features


    def get_gene_id2feature_id_mapping(self, gene_ids):
        gene_id2feature_id = {}
        for g in gene_ids:
            self.open_cursor.execute(f"""
            SELECT feature_id FROM features WHERE name=?;
            """, (g,))
            gene_id2feature_id[g] = self.open_cursor.fetchall()[0][0]
        return gene_id2feature_id


    def get_features_data(self, feature_id):
        fields = ("contig_name", "feature_id", "name", "source","feature","start", "end", "score", "strand",  "frame", "amino_acids", "nucleotides", "annotations")
        ff = ",".join([ f'"{f}"' for f in fields ])
        self.open_cursor.execute(f"""
        SELECT {ff} FROM features WHERE feature_id=?;
        """, (feature_id,))
        features = [ { k: v for k,v in zip(fields,m)}  for m in self.open_cursor.fetchall()][0]

        return features


    def get_features_seq(self, features, type = "amino_acids"):
        if type not in seq_types:
            raise WrongSeqTypeError("Trying to get a sequence-'type' from the feature table that don't exist, if should be in '{seq_types}' ")

        values = ",".join([f"'{c}'" for c in features])
        if len(features) == 1:
            self.open_cursor.execute(f"""
        SELECT feature_id,{type} FROM features WHERE feature_id ={values};
        """)
        else :
            self.open_cursor.execute(f"""
        SELECT feature_id,{type} FROM features WHERE feature_id = IN ({values});
        """)
        seqs = {m[0] : m[1] for m in self.open_cursor.fetchall()}
        return seqs
    
    def get_contigs(self, contigs):
        values = ",".join([f"'{c}'" for c in contigs])
        if len(contigs) == 1:
            self.open_cursor.execute(f"""
            SELECT contig_name,sequence FROM contigs WHERE contig_name={values};
            """)
        else :
            self.open_cursor.execute(f"""
            SELECT contig_name,sequence FROM contigs WHERE contig_name IN ({values});
            """)
        data = self.open_cursor.fetchall()
        return {c[0] : "".join(c[1]) for c in data}

    def get_genome_feats(self, genome_name, type = "amino_acids"):
        if type not in seq_types:
            raise WrongSeqTypeError("Trying to get a sequence-'type' from the feature table that don't exist, if should be in '{seq_types}' ")
        self.open_cursor.execute(f"""
        SELECT feature_id,{type} FROM features INNER JOIN contigs  on contigs.contig_name = features.contig_name WHERE genome_name = ? ;
        """, (genome_name,))
        seqs = {m[0] : m[1] for m in self.open_cursor.fetchall() if m[1]}
        return seqs

    def get_genome_GCs(self, genome_name):
        self.open_cursor.execute(f"""
    SELECT gc_name FROM genome2gcs
     WHERE genome2gcs.genome_name = ?;
        """, (genome_name,))
        return [m[0] for m in self.open_cursor.fetchall()]

    def get_genome_from_GC(self, gc_name):
        self.open_cursor.execute(f"""
    SELECT genome_name FROM genome2gcs
     WHERE genome2gcs.gc_name = ?;
        """, (gc_name,))
        return [m[0] for m in self.open_cursor.fetchall()]

    def get_genes_from_GC(self, gc_name):
        self.open_cursor.execute(f"""
    SELECT feature_id FROM gc2feature
     WHERE gc2feature.gc_name = ?;
        """, (gc_name,))
        return [m[0] for m in self.open_cursor.fetchall()]



    def genome_has_GCs(self, genome_name):
        self.open_cursor.execute(f"""
    SELECT gc_name FROM genome2gcs WHERE
    genome2gcs.genome_name = ? LIMIT 1);
        """, (genome_name,))
        return len(self.open_cursor.fetchall()) > 0


    def get_mOTU_GCs(self, motu_name):
        self.open_cursor.execute(f"""
    SELECT genome2gcs.genome_name, genome2gcs.gc_name 
    FROM genome2gcs INNER JOIN genome2motu ON genome2motu.genome_name = genome2gcs.genome_name 
    WHERE motu_name = ?;
        """, (motu_name,))
        data = self.open_cursor.fetchall()
        odat = {m[0] : [] for m in data}
        for m in data:
            odat[m[0]] += [m[1]]
        return odat



    def get_genome(self, genome_name):
        self.open_cursor.execute(f"""
        SELECT contig_name,sequence FROM contigs WHERE genome_name='{genome_name}';
        """)
        data = self.open_cursor.fetchall()
        return data

    def get_random_genome(self):
        from mOTUlizer.classes.MetaBin import MetaBin

        self.open_cursor.execute(f"""
        SELECT genome_name FROM genomes ORDER BY RANDOM() LIMIT 1;
        """)
        data = self.open_cursor.fetchall()

        return MetaBin(name = data[0][0])


    def get_signature(self, genome_name):
        self.open_cursor.execute(f"""
        SELECT signature FROM genomes WHERE genome_name='{genome_name}';
        """)
        data = self.open_cursor.fetchall()[0][0]
        return list(load_signatures(data))[0]



    def get_genome_info(self, genome_name):
        fields = ["genome_name", "source", "annotations", "taxonomy", "completeness", "redundancy"]
        self.open_cursor.execute(f"""
        SELECT {' , '.join(fields)} FROM genomes WHERE genome_name=?;
        """, (genome_name,))
        data = self.open_cursor.fetchall()
        return {k : v for k,v in zip(fields, data[0])}


    def write_genome_fasta(self, name : str, file = None, width = 60, type = "nucleotides"):
        if type == "nucleotides":
            data = self.get_genome(name)
        else :
            data = list(self.get_genome_feats(name).items())
        if width:
            chop = lambda d : "".join([d[1][i:i+width] + "\n" for i in range(0, len(d[1]), width)])
        else :
            chop = lambda d : d[1] + "\n"
        text = [f">{d[0]}\n{chop(d)}" for d in data]
        if file:
            with open(file, "w") as handle:
                handle.writelines(text)
        else :
            return "\n".join(text)

    def write_genome_gff(self, name : str, file = None, with_ID_prefix = "", CDS_only = True, ignore_compounded = True) :
        fields = ["features.contig_name", "source", "feature", "start", "end", "score", "strand", "frame", "features.annotations", "feature_id", "name"]
        keep_fields = ["contig_name", "source", "feature", "start", "end", "score", "strand", "frame", "annotations"]

        self.open_cursor.execute(f"""
        SELECT {' , '.join(fields)}  FROM features INNER JOIN contigs  on contigs.contig_name = features.contig_name WHERE genome_name = ? ;
        """, (name,))
        seqs = [{k.replace('features.','') : v for k,v in zip(fields, m)} for m in self.open_cursor.fetchall()]
        
        if CDS_only:
            seqs = [v for v in seqs if v['feature'] == "CDS"]
        if ignore_compounded: 
            seqs = [v for v in seqs if v['start'] and v['start'] != "None"]            
        for s in seqs: 
            s['annotations'] = json.loads(s['annotations']) 
            s['annotations']['original_ID'] = s['annotations'].get('ID', s['annotations'].get('locus_tag', "None")) 
            s['annotations']['internal_id'] = s['feature_id']
            s['annotations']['ID'] = f"{with_ID_prefix}{s['feature_id']}"
            s['annotations']['locus_tag'] = s['feature_id']
            s['source'] = "-"
            if s['score'] == "":
                s['score'] = "."
            if s['frame'] == "":
                s['frame'] = "0"    
            del s['feature_id']
            del s['name']
            if "protein_id" in s['annotations']:
                del s['annotations']['protein_id']
            keys_ = list(s['annotations'].keys())
            keys_.remove("ID")
            keys_ = ["ID"] 
            s['start'] = s['start'] 
            s['annotations'] = ";".join([ f"{k}={s['annotations'][k]}" for k in keys_ ])
        contigs_details = "\n".join(sorted([f"##sequence-region {k} 1 {len(v)}"  for k,v in self.get_genome(name)]))
        features = "\n".join([ "\t".join([f"{s[f]}" for f in keep_fields])  for s in seqs])
        contigs = self.write_genome_fasta(name = name)
        text = f"""##gff-version 3
{contigs_details}
{features}
##FASTA
{contigs}
        """
        if file:
            with open(file, "w") as handle:
                handle.writelines(text)
        else :
            return text
    

    def get_taxonomy(self, genome_ids = None):
        if genome_ids:
            where_line = "WHERE genome_name IN ({names})".format(names = " , ".join([f'"{g}"' for g in genome_ids]))
        else :
            where_line = ""

        self.open_cursor.execute(f"""
        SELECT genome_name, taxonomy FROM genomes {where_line};
        """)
        data = {k : v for k,v in  self.open_cursor.fetchall()}
        data = {k : json.loads(v.replace('""','"')) for k,v in data.items() if v}
        return data

    @bunched_query
    def add_core(self, motu, gcs, likelies, bunched = False):
        self.open_cursor.executemany(f"""
        INSERT INTO cores (motu_name, gc_name, loglikelihood)
        VALUES (?, ?, ?);
        """, [(motu,f, likelies[f]) for f in gcs])

    def get_core(self, motu_name):
        self.open_cursor.execute(f"""
        SELECT gc_name FROM cores WHERE motu_name = ?;
        """, (motu_name,))
        data = [ gc[0] for gc in self.open_cursor.fetchall()]
        return data if len(data) > 0 else None

    def get_likelies(self, motu_name):
        self.open_cursor.execute(f"""
        SELECT gc_name, loglikelihood FROM cores WHERE motu_name = ?;
        """, (motu_name,))
        data = {gc[0] : float(gc[1]) for gc in self.open_cursor.fetchall()}
        return data if len(data) > 0 else None


    def reset_gcs(self):
        self.open_cursor.execute("DELETE FROM cores;")
        self.open_cursor.execute("DELETE FROM gc2feature;")
        self.open_cursor.execute("DELETE FROM gene_clusters;")
        self.open_cursor.execute("DELETE FROM genome2gcs;")
        self.commit()

    def check_db(self):
        pass
