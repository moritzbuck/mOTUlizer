---
title: "mOTUlizer"
author: Moritz Buck
date: August, 2020
geometry: margin=2cm
output: pdf_document
---

# Genome Informatics Lightning talk --- Script #

## mOTUlizer: Bayesian approach to leverage metagenomic bins for pan-genome analysis ##

Typicaly, in a large scale genome resolved metagenomic study, we make many different assemblies, from which we cluster many microbial genomes for each clade. These are incomplete and contaminated, so only the best genome, as in most complete, and least contaminated, of a specific population will be picked, losing a lot of information. However, not every genome is "equaly" incomplete, the presence/abscence pattern of genes in a set of genomes can give us information. With mOTUpan, we leverage this information to decide if any gene is either in every genome of this clade, so part of the core-genome of that clade, or only some, hence in the accessory genome. This method is fast and flexible and can be used on large data easily!
