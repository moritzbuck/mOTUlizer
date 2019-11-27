# Bacterial pangenomes are constrained by genome size

Understanding the population structure of wild microbial populations has been hard, most studies have been done on isolates, and more specifically isolates of a medical nature. The population of these will significantly differ from wild aquatic populations, where microbes have to survive in much poorer environment and are subjected to less bottlenecking events.
Recent advances in sequencing and bioinformatics allow us finally to have better insight into a wild aquatic communities through metagenomics. Similarly
In this analysis we have used a set of 15000 metagenomic assembled genomes (MAGs) clustered into around 1500 metagenomic operational taxonomic units (e.g. species) from a study involving 400 samples from 40 different freshwater bodies.

This set of MAGs has been used to compute core genomes and pangenomes for each of these bacterial species. This analysis comes to the conclusion that the number of genes available for a species in it's pangenome is constrained by the size of the genomes, allowing larger genomes to access to even more niches then previously thought.


# Main

Size and diversity of microbial pangenomes is an open question, the data and the tools necessary have not been available for very long. we have used a newly developed statistical method for the computation of core genomes to estimate core genomes and auxiliary genomes for ~1200 (meta-)genomic Operational Taxonomic Units (mOTUs) from public repositories and ~300 mOTUs from our own dataset of freshwater metagenomes.

For the selected ~1500 mOTUs we have computed core genomes, as in cluster of orthologous genes (COGs) that are likely to be present in all members of this mOTU. Conversely all genes found in these genomes that are not in the core must be part of the auxiliary genomes

Computing core-genomes for has been limited in taxonomic scope, as most methods require a number of high-quality genomes of closely related organisms to have an accurate estimate of presence abscence in the population of interest. We present here a novel bayesian method for the computation of core genomes relying on the presence-abscence of genes/COGs/annotations in sets of draft (or complete) genomes. The method is wrapped in a tool available at this webpage: www.github.com/moritzbuck/0039_mOTUlizer.

Each genome (we will use genome as a shorthand for any set of nucleotide sequences belonging to the same biological entity, e.g. draft genome, complete genome, or Metagenome Assembled Genome), is described as a set of traits, in the case of this analysis a set of COGs, but mOTUlizer is agnostic to the type of traits. And each genome it self is part of a set of genomes which we will call an mOTU (metagenomic Operational Taxonomic Unit, due to the nature of the data we analyse here). mOTUlizer uses an interative approach to classify each trait of the genome in an mOTU as a "core"-trait or "auxiliary"-trait based on a likelihood ratio. For each of the two hypotheses (core-trait or auxiliary-trait) a probability is computed assuming a certain completeness values for the each genome. Whicever of these is more likely is picked as class for that trait. Using this new classification we update the completeness estimate and recalculate the likelyhood ratios and repeat the process until convergence.

(@) $$p_{\textrm{core,cog}}=\prod_{\substack{\textrm{g in mOTU}\\ \textrm{if cog in g}}} (1-\bar{p}_{\textrm{core,cog,g}}) \prod_{\substack{\textrm{g in mOTU}\\ \textrm{if cog not in g}}} \bar{p}_{\textrm{core,cog,g}}$$

(@) $$\bar{p}_{\textrm{core,cog,g}}=(1-\frac{|cog|}{|G|})^{len(g)-(lc_{g}c_{g})}$$





 To account for possible errors both in gene-prediction and metagenomic binning, all COGs appearing only once have been removed.

The selected ~1500 mOTUs have more than 6 genomes (from isolates or as metagenome assembled genome, e.g. MAG),
