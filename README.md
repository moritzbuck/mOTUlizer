# mOTUlizer
Utility to analyse a group of closely related MAgs/Genomes/bins/SUBs of more or less dubious origin

## INSTALL

```
pip install mOTUlizer
```


## USAGE

### EASY

To make OTUs and get some stats, needs fastANI in the `PATH`, and output of checkm

```
mOTUlize.py -k checkm_output.txt  --output a_messy_json-file_with_the_output.json
```

Loads of little options if you do : `mOTUlize.py -h`

Also there is `mOTUpan.py` that can compute core genomes and pangenomes. Needs to be more debugged so try  out at own risk:

```
mOTUpan.py -h
```
