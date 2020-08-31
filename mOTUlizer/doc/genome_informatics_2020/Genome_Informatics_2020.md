---
title: "mOTUlizer"
author: Moritz Buck
date: December, 2019
geometry: margin=2cm
output: pdf_document
---

# mOTUlizer: Bayesian approach to leverage metagenomic bins for pan-genome analysis

Understanding the pan-genomic structure of environmental microbe populations is challenging. Most of the tools available for pan-genome analyses rely on high quality genomes that are just not available for environmentally-relevant microbes. Also these environmental pan-genomes are significantly different from the more studied pan-genomes of host-associated microbes as they tend to be less prone to bottlenecking events, and more prone to complex selection pressures.
The recent advances in sequencing and bioinformatics finally give us access to many more genomic bins from environmental clades. However, these often suffer from incompleteness, and most methods for analyzing pan-genomes are very sensitive to this issue.
We propose a new method of computing core and accessory genes explicitly relying on the incompleteness of metagenomic bins for pan-genome analysis. By default, we compute clusters of orthologous genes (COGs) using the diamond and silix tools for the target set of bins.
The likelihood of each COG to belong to the core-genome or the accessory-genome is estimated by computing the probability of its presence/absence pattern in the target set of bins, assuming the estimated incompleteness of each bin. COGs are then tagged accordingly, completeness reestimated, and the processed iterated until convergence.
This procedure is agnostic to the traits it classifies, and can be used for any set of genomically encoded trait not only COGs it computes by default. If traits are pre-computed, the pan-genome analysis set (i.e., membership of the trait in the core or in the accessory genome) is computationally efficient and can be scaled up to large sets of bins. The method has been tested on simulated data to show high recovery and specificity.  
Additionally, we have applied mOTUlizer to all species of the GTDB, including amongst others the 9443 Staphylococcus aureus genomes, 9072 Escherichia flexneri genomes and 8679 Salmonella enterica genome. We also ran this tool on around 100Gb of bins from a large freshwater metagenomic dataset to obtain around 3600 microbial species. These analyses come to the conclusion that the number of genes available for a species in its pan-genome is constrained by the size of the genomes. This approach is implemented in the mOTUlizer python package (available at github.com/moritzbuck/mOTUlizer).
