# mOTUlizer
Utility to analyse a group of closely related MAgs/Genomes/bins/SUBs of more or less dubious origin

## INSTALL

In the future (outdated versino right now)

```
pip install mOTUlizer
```

Now:

```
python setup.py install
```


## USAGE

### EASY

To make OTUs and get some stats, needs fastANI in the `PATH`, and output of checkm

```
mOTUlize.py -k checkm_output.txt  --output a_messy_json-file_with_the_output.json
```

Loads of little options if you do : `mOTUlize.py -h`

Also there is `mOTUpan.py` that can compute core genomes and pangenomes, likelihood of a gene to be in all genomes of a set or to only be in some. An intro video here:

[![mOTUpan for beginners](https://img.youtube.com/vi/VIeV1Gg5NS4/0.jpg)](https://www.youtube.com/watch?v=VIeV1Gg5NS4)

Needs to be more debugged so try  out at own risk:

```
mOTUpan.py -h
```

Simplest command to run (needs mmseqs2 installed), but many options:

```
python mOTUlizer/bin/mOTUpan.py --faas *.faa
```
