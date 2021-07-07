#!/usr/bin/env python
# -*- coding: utf-8

import sys

import anvio
import anvio.dbops as dbops
import anvio.terminal as terminal
import anvio.filesnpaths as filesnpaths
import anvio.genomestorage as genomestorage

from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError
from anvio.tables.miscdata import TableForItemAdditionalData


__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2018, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Moritz Buck"
__email__ = "moritz.buck@slu.se"
__requires__ = ["pan-db", "genomes-storage-db"]
__provides__ = []
__resources__ = [("mOTUpan github", "https://github.com/moritzbuck/mOTUlizer/")]
__description__ = "runs mOTUpan for your genome-set"


run = terminal.Run()
progress = terminal.Progress()


def main(args):
    if not args.output_file and not args.store_in_db:
        if args.just_do_it:
            run.warning("Lol. You are making anvi'o compute stuff that will not be stored anywhere. Weird you :(", lc="green")
        else:
            raise ConfigError("But why no select a way to report this stuff? :/ You can ask anvi'o to store your results as a "
                              "a TAB-delimited file. Or add them to directly to the database. If you insist that you just want "
                              "to run this analysis without really storing anything anywhere you can always use `--just-do-it` "
                              "(because why not, it is your computer).")

    if args.output_file:
        filesnpaths.is_output_file_writable(args.output_file)

    pan = dbops.PanSuperclass(args, r=terminal.Run(verbose=False))

    # get a pan instance
    pan = dbops.PanSuperclass(args)

    # complain if there is no genomes storage
    if not pan.genomes_storage_is_available:
        raise ConfigError("The anvi'o pan class does not see a genomes storage. No genomes storage no cake.")

    # get a genome storage instance
    genome_storage = genomestorage.GenomeStorage(pan.genomes_storage_path, run=terminal.Run(verbose=False))

    # recover genome completeness values
    genome_completeness_dict = {g: v.get('percent_completion', None) for g, v in genome_storage.genomes_info.items()}
    genome_lengths_dict = {g: v.get('total_length', None) for g, v in genome_storage.genomes_info.items()}

    # initialize gene clusters
    pan.init_gene_clusters()

    # build the gene clusters dict for mOTU
    gene_clusters_dict = {g : set() for g in pan.genome_names}
    for gc, hits in pan.gene_clusters.items():
        for genome, genes in hits.items():
            if len(genes) > 0:
    completess = { 'empty' : None }
    if pan.genomes_storage_is_available:
        genome_storage = genomestorage.GenomeStorage(pan.genomes_storage_path, run=terminal.Run(verbose=False))
        completess = { g : v.get('percent_completion', None) for g,v in genome_storage.genomes_info.items()}

    if any([v is None for v in completess.values()]):
        run.warning("Lol. You are running mOTUpan wihout prior completeness estimate. This is probably either because you didn't run anvi-run-scg-something or you run this wihout genome-storage. It's ok, we will use the longest genome as '100%' complete and the completeness of the others with be estimated by their relative lenghts to that one. It might be a bit worse, but you'll get new completnesses as estiamted by mOTUpan anyhow!", lc="green")
        completess = "length_seed"
                gene_clusters_dict[genome].add(gc)


    motu = mOTU( name = "mOTUpan_core_prediction" , faas = {} , cog_dict = genome2genecluster, checkm_dict = completess, max_it = 100, threads = args.num_threads, precluster = False, method = 'default')
    if args.num_bootstraps :
        motu.roc_values(boots=args.num_bootstraps)
    if args.store_in_db:
        data2add = { g : {'partition' : 'core' if g in motu.core else 'accessory', 'log_likelihood_ratio' : motu.likelies[g]} for g in all_gene_clusters}
        panditional_data = TableForItemAdditionalData(args)
        panditional_data.add(data2add, data_keys_list = ['partition','log_likelihood_ratio'])

    if args.output_file:
        with open(args.output_file, "w") as handle:
            print(motu.pretty_pan_table(), file = handle)
        run.info("Output file", args.output_file, mc="green")

    if args.num_bootstraps > 0 :
        print(motu.roc_values(0))
#        run.info("test test", value = , mc="green")


if __name__ == '__main__':
    parser = ArgumentParser(description=__description__)

    groupA = parser.add_argument_group('INPUT FILES', "Input files from the pangenome analysis.")
    groupA.add_argument(*anvio.A('pan-db'), **anvio.K('pan-db'))
    groupA.add_argument(*anvio.A('genomes-storage'), **anvio.K('genomes-storage', {'required': False}))

    groupB = parser.add_argument_group('REPORTING', "How do you want results to be reported? Anvi'o can produce a TAB-delimited output file for\
                                                     you (for which you would have to provide an output file name). Or the results can be stored\
                                                     in the pan database directly, for which you would have to explicitly ask for it. You can get\
                                                     both as well in case you are a fan of redundancy and poor data analysis practices. Anvi'o\
                                                     does not judge.")
    groupB.add_argument(*anvio.A('output-file'), **anvio.K('output-file'))
    groupB.add_argument(*anvio.A('store-in-db'), **anvio.K('store-in-db'))

    groupC = parser.add_argument_group('SELECTION', "Which gene clusters should be analyzed. You can ask for a single gene cluster,\
                                       or multiple ones listed in a file, or you can use a collection and bin name to list gene clusters\
                                       of interest.")
    groupC.add_argument(*anvio.A('collection-name'), **anvio.K('collection-name'))
    groupC.add_argument(*anvio.A('bin-id'), **anvio.K('bin-id'))

    groupD = parser.add_argument_group('OPTIONAL',"Optional stuff available for you to use")
    groupD.add_argument(*anvio.A('num-threads'), **anvio.K('num-threads'))
    groupD.add_argument(*anvio.A('just-do-it'), **anvio.K('just-do-it'))
    groupD.add_argument( '-B', '--num-bootstraps' , metavar = "NUM_BOOTSTRAPS", default = 0, type = int, help = "number of boostraps run on the partitioning to evaluate it's quality"  )

    args = parser.get_args(parser)

    try:
        main(args)
    except ConfigError as e:
        print(e)
        sys.exit(-1)
    except FilesNPathsError as e:
        print(e)
        sys.exit(-1)
