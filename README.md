# FluPipe

[![Twitter Follow](https://img.shields.io/twitter/follow/rki_de.svg?style=social)](https://twitter.com/rki_de)
[![DOI](https://zenodo.org/badge/652663468.svg)](https://zenodo.org/doi/10.5281/zenodo.13684139)



## 1. Introduction

FluPipe provides a fully automated, flexible and reproducible workflow for reconstructing genome sequences from Illumina NGS data. The pipeline is optimized for Influenza data.

## 2. Setup 

The most convenient way to install the pipeline is by using git and [conda](https://docs.conda.io/en/latest/miniconda.html):

```bash
# installing the pipeline using git
cd designated/path
git clone https://github.com/rki-mf1/FluPipe.git/
cd FluPipe
conda env create -f flupipe.yml -n FluPipe
conda activate FluPipe
```

# 3. Usage

As a minimum the pipeline needs the following input:
- folder containing gz-compressed FASTQ files (`-d`)
- output folder (`-o`), in which a subfolder named results is automatically created to store all results
- a reference sequence (`--ref`) or an influenza segment database, containing representative fasta files for each genome segment (`--segmentdb`)

```bash
# activate conda environment once before using the pipeline 
conda activate FluPipe 


flupipe.py    -d path/to/myInputFolder \
              -o path/to/myOutputFolder \
              --segmentdb path/to/segmentdb
```

# 4. Options to customize the workflow

The manual page provides information on all options available.

```bash
flupipe.py --help
```


## 4.1 Adjusting Read Filtering

### 4.1.1 Read Length

Per default the minimum read length filter is set to 50. 

(-l 50) 

### 4.1.2 Read Quality

Qualitative read quality used by fastp to filter reads.
By default the "--read_filter_qual" option uses a phredscore of 20 as a cutoff.

## 4.2 Taxonomic Read Filtering

If necessary, reads not derived from the Orthomyxoviridae family can be excluded. 
Read classification is based on corresponding k-mer frequencies using a defined kraken2 database (`--kraken`). 
A database containing Influenza A , Influenza B and human genome sequences is recommended. 

## 4.3 Find a reference for each segment

For each segment, a multifasta file with any number of reference sequences can be provided. The pipeline compares the sequencing reads to the given references and determines the optimal reference sequence per segment for the given data based on read coverage, read depth, and uniformity of mapping. 

## 4.4 Adapting variant calling

Sites considered for variant calling can be restricted based on the following parameters at the respective position.

- the minimum sequencing depth (`--vvar_mincov`; default: 20) 
- the minimum number of reads supporting a variant (`--var_call_count`; default: 10)
- the relative number of reads supporting a variant (`--var_call_frac`; default: 0.1) 

CHECK PARAMETERS


## 4.5 Consensus generation

When generating the consensus sequence, all positions whose read coverage is below a defined threshold can be hard-masked by N (`--cns_min_cov`; default: 20). 
In addtion, genotypes can be adjusted meaning that variants supported by a given fraction of all reads covering the respective site are called explicitely (`--cns_gt_adjust`; default: 0.9). 
This means that a variant that shows a read fraction of 0.94 would be set to full alternate allele and variants showing only 0.03 readfraction are changed to reference.

