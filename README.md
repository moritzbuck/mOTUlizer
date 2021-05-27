# mOTUlizer

Utility to analyse a group of closely related MAGs/Genomes/bins/SUBs of more or less dubious origin. Right now it is composed of three programs:

* `mOTUlize.py` takes a set of genomes (I will use the term genome as a short hand for set of nucleotide sequences that presumably come from the same organism/population, can be incomplete, redundant or contaminated) and cluster them in to metagenomic Operational Taxonomic Units (mOTUs). Using similarity scores (by default ANI as computed by fastANI, but user can provide other similarities) a network is built based on (user defined) better quality genomes (for historical reasons called MAGs) by thresholding the similarities at a specific value (95% by default). The connected components of this graph are the mOTUs. Additionally lower quality genomes (SUBs, ) are recruited to the mOTU of whichever MAG they are most similar too if the similarity is above the threshold.

* `mOTUpan.py` computes the likelihood of gene-encoded traits to be expected in all of a set of genomes, e.g. of a trait to be in the core genome of a set of genomes (of possibly varying quality). Basically you provide to `mOTUpan` the set of proteoms of your genomes of interest (for example from the same mOTU or Genus) as well as a completeness prior of these genomes (for example [`checkm`](https://ecogenomics.github.io/CheckM/) output or a fixed value) and it computes gene clusters using [`mmseqs2`](https://github.com/soedinglab/MMseqs2), you can also provide your own genome encoded traits either as a `JSON`-file, or `TAB`-separated file (see example files). For each of these gene-clusters it will then compute the likelihood of it being in the core vs the likelihood of it not being, the ratio of these likelihoods will determine if a trait is considered core or not. This new partitioning can be used to update our completeness prior, and recomputed iteratively until convergence.

* `mOTUconvert.py` converts the output of diverse programs into input files for `mOTUpan.py`, currently includes methods for [`mmseqs2`](https://github.com/soedinglab/MMseqs2), [`roary`](https://sanger-pathogens.github.io/Roary/), [`PPanGGOLiN`](https://github.com/labgem/PPanGGOLiN), [`eggNOGmapper`](https://github.com/eggnogdb/eggnog-mapper).

## INSTALL

With pip (right now not the right version, try again in a few days for hopefully version `0.2.0`)

```
pip install mOTUlizer
```

manually:

```
git clone https://github.com/moritzbuck/mOTUlizer.git
cd mOTUlizer
python setup.py install
```


## USAGE

### mOTUlize

To make OTUs and get some stats, needs fastANI in the `PATH`, and output of checkm

```
mOTUlize.py -k checkm_output.txt --fnas myfolderwithgenomes/*.fna
```

Loads of little options if you do : `mOTUlize.py -h`


### mOTUpan

An intro video [here](https://www.youtube.com/watch?v=VIeV1Gg5NS4):

[![mOTUpan for beginners](https://img.youtube.com/vi/VIeV1Gg5NS4/0.jpg)](https://www.youtube.com/watch?v=VIeV1Gg5NS4)


```
mOTUpan.py -h
```

Simplest command to run (needs mmseqs2 installed), but many options:

```
mOTUpan.py --faas *.faa
```

### mOTUconvert
