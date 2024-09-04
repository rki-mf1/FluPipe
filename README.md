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
conda create -f flupipe.yml -n FluPipe
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
+++++++++++++++++++++
from here README NEEDS WORK
Add description on how to build the database. + TaxID information

## 4.3 Adapting variant calling

Sites considered for variant calling can be restricted based on the following parameters at the respective position.

- the minimum sequencing depth (`--vvar_mincov`; default: 20) 
- the minimum number of reads supporting a variant (`--var_call_count`; default: 10)
- the relative number of reads supporting a variant (`--var_call_frac`; default: 0.1) 

CHECK PARAMETERS


## 4.4 Adapting variant filtering

Variants can be excluded, if the mean mapping quality of observed alternate alleles does not meet a given threshold (`--var_filter_mqm`; default: 40).
The mapping quality measures how good reads align to the respective reference genome region. Good mapping qualities are around MQ 60. 
The filter can be deactivated by setting it to 0.
Variants exceeding the value will be retained.

Additionally, variants can be filtered if the strand balance probability for the alternate allele exceeds a given threshold (`--var_filter_sap`; default: -1, i.e. inactive).
The strand balance probability is a Phred-scaled measure of strand bias. A value near 0 indicates little or no strand bias. 
For amplicon data sets this value can be very high due to unequal performance of neighbouring amplicons. 
Since most of the data processed with covpipe is amplicon data, the default behaviour is to *not* filter variants with high strand bias values.
A recommended value for whole genome approaches, e.g. given bey GATK, is 60.
Values smaller than the set value will be retained.


Finally, variants can be filtered by the site quality (`--var_filter_qual`; default: 10) indicating the pobability for a polymorphism at the respective site.
To deactivate thsi filter, please set to 0.
Variants exceeding this value will be retained.

## 4.5 Refernce options
--ref vs segment db 


## 4.6 Adapting consensus generation

When generating the consensus sequence, all positions whose read coverage is below a defined threshold can be hard-masked by N (`--cns_min_cov`; default: 20). 
In addtion, genotypes can be adjusted meaning that variants supported by a given fraction of all reads covering the respective site are called explicitely (`--cns_gt_adjust`; default: 0.9). 
This means that a variant that shows a read fraction of 0.94 would be set to full alternate allele and variants showing only 0.03 readfraction are changed to reference.

```bash
# hard-masking sites covered by less than 10 reads and explicitely call variants
# supported by at least 95% of all reads at the respective site
ncov_minipipe --cns_min_cov 10 \
              --gt_adjust 0.95 \
              --reference path/to/myReference.fasta \
              --input path/to/myInputFolder \
              -o path/to/myOutputFolder
```




# 5 Workflow

![workflow](docs/workflow.png "workflow")

The illustration shows a slightly simplified workflow that can be divided into the following sections: 


## 5.1 Pre-processing

To ensure the quality of sequence data used for reconstruction, the workflow provides: 

steps: 
- 5' clipping of amplification primer sequences (optional)
- 3' removal of
  - Illumina adapter sequences (optional)
  - low complexity stretches
  - low quality stretches
- excluding of non-SARS-CoV-2 reads (optional)
  - to our experience this step improves analyses quality a lot and hence we recommend it

used tools:
- bamclipper
- fastp
- kraken2 

## 5.2 Mapping

Pre-processed reads are aligned to the reference sequence using bwa-mem (default parameters).

## 5.3 Variant calling, filtering, genotype adjustment and annotation

Variants are called using freebayes allowing the following customizations:
  - minimal total coverage (`--var_call_cov`; default: 20) 
  - absolute number of variant supporting reads (`--var_call_count`; default: 10)
  - relative number of variant supporting reads (`--var_call_frac`; default: 0.1) 

The called variants can be further filtered using bcftools. For this custom thresholds can be applied to  
  - mean mapping quality of observed alternate alleles (MQM, `--var_filter_mqm`; default: 40)
  - strand balance probability for the alternate allele (SAP, `--var_filter_sap`; default: 60)
  - polymorphism probability (QUAL, `--var_filter_qual`; default: 10)

As an unique feature, the worklow allows to adjust the genotype of certain sites by explicitely calling variants supported by a given fraction of reads covering the respective site.

By default, all variants are inspected regarding their impact on coding sequences using SNPeff. 
This only works, if  your reference genome is listed in the SNPeff supported genomes (see Appendix below). Otherwise please deactivate this feature (`--no-var-annotation`).

## 5.4 Consensus generation

Using bcftools the consensus is created allowing to mask lowly covered regions by N (user-defined threshold).
The default consensus (suffix iupac) follows iupac nomenclature at sites mixed variants have been detected 
(and were not adjusted by genotype as described before).
In addition, a second consensus is created for each dataset (suffix masked) where these ambiguous sites are
masked by N (additionally to the lowly covered regions).

## 5.5 QC-report

In addition to the consensus sequences, a HTML-based report is generated summarizing different quality measures and mapping statistics for each dataset such as:

- runID (optional))
- conditional table that warns the user, if samples that were identified as negative controls show high reference genome coverage
- table of read properties:
  - number of bases (before / after trimming)
    - if amplification primer clipping was done
  - length of reads (before / after trimming)
    - if amplification primer clipping was done
  - number of bases mapped (Q = 30)
- optional table listing the species filtering results emitted by Kraken2
- table of mapping properties:
  - reads mapped to reference genome (number & fraction of input)
  - median / sd  of fragment size
- genome wide plot of coverage
- histogram of fragment sizes

- table of reference coverage characteristics
- table of lineage assigments using pangolin (optional)

Samples showing at least 20X sequencing depth at more than 95% of the reference genome are designated a successful genome sequencing.

