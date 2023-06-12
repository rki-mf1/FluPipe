# IMPORT MODULES
import re
import os
import pprint
import sys
import yaml


# DEBUGGING
## use as DEBUG(variable) for debugging at runtime
pp = pprint.PrettyPrinter(indent=4)
DEBUG = pp.pprint


def default_if_not(key, _dict, default):
    try:
        return _dict[key]
    except KeyError:
        return default

# INPUT CHECK
## preparing checks
def checkConfigKey(key, config):
    if key in config and config[key] not in ["", False, None]:
        return True
    return False

def isFileConfig(key, config):
    if not checkConfigKey(key, config) or not os.path.isfile(config[key]):
        return False
    return True

err = []

## check if config was provided
if not config:
    err.append("input error: no config data provided.")

## checking sample yaml file [mandatory file]
if not isFileConfig('samples', config):
    err.append("input error: missing or invalid sample file definition.")

## checking reference file [mandatory file]
# if 'reference' not in config or not config['reference']:
#     config['reference'] = os.path.join(workflow.basedir, "data", "KM517573.1.fasta")
# if not isFileConfig('reference', config):
#     err.append("input error: reference file does not exists")

## checking primer files [optional file]
if not checkConfigKey('primer', config):
    sys.stderr.write("note: amplicon primer clipping turned off\n")

elif not os.path.isfile(config['primer']):
    err.append("input error: primer file does not exists")

## checking adapter file [optional file]
if not checkConfigKey('adapter', config):
    sys.stderr.write("note: adapter clipping turned off\n")
elif not os.path.isfile(config['adapter']):
    err.append("input error: adapter file does not exists.")

## checking kraken database [optional folder]
if not checkConfigKey('krakenDb', config):
    sys.stderr.write("note: taxonomic read filtering turned off\n")
elif not os.path.isdir(config["krakenDb"]):
    err.append("input error: kraken database does not exists.")

## checking segment database [optional folder]
if not checkConfigKey('segmentdb', config):
    sys.stderr.write("note: automatic reference detection turned off - no segment database\n")
elif not os.path.isdir(config["segmentdb"]):
    err.append("input error: segment database does not exists.")

## checking annotation file [optional file]
if not checkConfigKey('var_annotation', config):
    sys.stderr.write("note: variation inspection turned off\n")

if not checkConfigKey('cns_annotation', config):
	sys.stderr.write("note: consensus genome annotation turned off\n")
elif not os.path.isfile(config['cns_annotation']):
    err.append("input error: gff file does not exist.")

## input error reporting
if err:
    sys.stderr.write("\n".join(err) + "\n")
    sys.exit(1)

# CONSTANT DEFINITIONS
## sample files
SAMPLES = dict()
with open(config["samples"], 'r') as handle:
    SAMPLES = yaml.safe_load(handle)
    SAMPLES = {str(x[0]): x[1] for x in SAMPLES.items()}
#segment names
ALL_SEGMENTS= ["HA", "MP", "NA", "NP", "NS", "PA", "PB1", "PB2"]

## additional files or folders
VAR_ANNOT = default_if_not("var_annotation", config, None)
RUN_ANNOT = isFileConfig('cns_annotation', config)
CNS_ANNOT = default_if_not("cns_annotation", config, None)
KRAKEN_DB = default_if_not("krakenDb", config, None)
SEGMENT_DB = default_if_not("segmentdb", config, None)
#KRAKEN_TAX_ID = default_if_not("krakenTaxID", config, None)

## read clipping parameters
PRIMER = default_if_not("primer", config, None)

## variant calling
VAR_CALL_COV = config['var_call_cov']
VAR_CALL_COUNT = config['var_call_count']
VAR_CALL_FRAC = config['var_call_frac']

## variant hard filtering
VAR_FILTER_MQM = config['var_filter_mqm']
VAR_FILTER_SAP = config['var_filter_sap']
VAR_FILTER_QUAL = config['var_filter_qual']

#reference sequence provided
REFERENCE = default_if_not("reference", config, None)


## read filtering
PCR_DEDUP = default_if_not("pcr_dedup", config, None)
READ_FILTER_QUAL = default_if_not("read_filter_qual", config, 20)
READ_FILTER_LEN = default_if_not("read_filter_len", config, 50)

## consensus generation
CNS_MIN_COV = config['cns_min_cov']
CNS_GT_ADJUST = default_if_not("cns_gt_adjust", config, None)

## reporting parameters
REPORT_RUNID = default_if_not("run_id", config, "")

## output folders
PROJFOLDER = os.path.join(config["output"], "results")
IUPAC_CNS_FOLDER = os.path.join(PROJFOLDER, "consensuses_iupac")
#MASKED_CNS_FOLDER = os.path.join(PROJFOLDER, "consensuses_masked")
DATAFOLDER = ["logs", "trimmed"]
if KRAKEN_DB:
    DATAFOLDER.extend(["classified", "filtered"])
DATAFOLDER.extend(["mapping"])
if PCR_DEDUP:
	DATAFOLDER.extend(["dedup"])
DATAFOLDER.extend(["mapping_stats", "variant_calling", "masking", "reporting"])
DATAFOLDER = { x[1]: os.path.join(PROJFOLDER, "intermediate_data", str(x[0]).zfill(2) + "_" + x[1]) for x in enumerate(DATAFOLDER) }

if not KRAKEN_DB:
    DATAFOLDER["classified"] = PROJFOLDER
    DATAFOLDER["filtered"] = PROJFOLDER


## files
#REFERENCE = os.path.join(DATAFOLDER["mapping"], "reference.fasta")
ADAPTERS = default_if_not("adapter", config, None) # the adpater file cannot be provided since it is copyright protected ILLUMINA!

## ref indexes
#PICARD_INDEX = os.path.splitext(REFERENCE)[0] + '.dict'
#SAMTOOLS_INDEX = REFERENCE + ".fai"
#BWA_INDEX = REFERENCE + ".bwt"

# SANITY CHECKS FOR SAMPLES & THRESHOLDS
## kraken input test
#if KRAKEN_DB and not KRAKEN_TAX_ID:
#    err.append("input error: kraken database defined but no TaxID")

## variant and cns threshold test
if VAR_CALL_FRAC < 0 or VAR_CALL_FRAC > 1:
    err.append("input error: the value of var_call_frac cannot be lower than 0 or greater than 1")
if CNS_MIN_COV < VAR_CALL_COV:
    err.append("input error: var_call_cov cannot be smaller than cns_min_cov.\nThey are {varcall} and {cns}".format(varcall=V, cns=CNS_MIN_COV))
if CNS_GT_ADJUST and (CNS_GT_ADJUST <= 0.5 or CNS_GT_ADJUST > 1):
    err.append("input error: the value of cns_gt_adjust has to be greater than 0.5 and not greater than 1")

## sanity error report
if err:
    sys.stderr.write("\n".join(err) + "\n")
    sys.exit(1)


# GENERAL FUNCTIONS
def getFastq(wildcards):
    return SAMPLES[str(wildcards.sample)]["read1"], SAMPLES[str(wildcards.sample)]["read2"]


# RULE ALL
def input_all(wildcards):
    files = []

    ## consensus
    for sample in SAMPLES:
        files.append(os.path.join(IUPAC_CNS_FOLDER, sample + ".iupac_consensus.fasta"))

    ## masked consensus
    #for sample in SAMPLES:
    #    files.append(os.path.join(MASKED_CNS_FOLDER, sample + ".masked_consensus.fasta"))

    ## variant annotation
    #if VAR_ANNOT:
    #    for sample in SAMPLES:
   #         files.append(os.path.join(DATAFOLDER["variant_calling"], sample, sample + ".annotation.html"))
    #if RUN_ANNOT:
        # annotate iupac consensus
     #   for sample in SAMPLES:
      #      files.append(os.path.join(IUPAC_CNS_FOLDER, sample + ".iupac_consensus.gff"))

        # annotate merged consensus
        #for sample in SAMPLES:
            #files.append(os.path.join(MASKED_CNS_FOLDER, sample + ".masked_consensus.gff"))

    ## report
    files.append(os.path.join(PROJFOLDER, "qc_report.html"))

    return files

rule all:
    input:
        input_all


# RULE IMPORT
## general rules
include: "rules/get_version.smk"
####
#reference detection
#minimap against segmentDBs
include: "rules/map_reads_minimapp2.smk"
#get refs fastas from segment DB
include: "rules/r_stats_minimapp2.smk"
include: "rules/seqkit_extract_ref.smk"
#map again + stats
include: "rules/best_Hits_map_reads_minimapp2.smk"
include: "rules/bam_stats_refDetec.smk"
#get final reference sequence for flupipe input
include: "rules/r_stats_final_ref.smk"
include: "rules/seqkit_extract_final_ref.smk"
include: "rules/combine_final_refs.smk"

#final perpare ref script
include: "rules/prepare_reference.smk"
####


## indexing
include: "rules/index_samtools.smk"
include: "rules/index_picard.smk"
include: "rules/index_bwa.smk"
include: "rules/index_bam.smk"
include: "rules/index_tabix.smk"

## amplicon primer clipping
include: "rules/trim_reads.smk"

## taxonomic read classification
include: "rules/classify_reads.smk"
include: "rules/filter_reads.smk"

## read mapping
include: "rules/map_reads.smk"
# PCR_DEDUP:
include: "rules/dedup_reads.smk"

include: "rules/sort_bam.smk"
include: "rules/get_read_cov.smk"

## variant calling
include: "rules/call_vars_lofreq.smk"
#include: "rules/bgzip_vcf.smk"

## consensus generation
include: "rules/create_consensus.smk"

## statistics
include: "rules/get_bamstats.smk"
include: "rules/get_insert_size.smk"
include: "rules/combine_readlength.smk"
include: "rules/getReadLength.smk"
## report
include: "rules/create_report.smk"

##extra rules to cover a specific variant case, good quality variant fails due to strand bias variant filter, most likely at the beginning or end of the sequence
#instead of the reference base an N is placed at this position
include: "rules/lofreq_secial_variant_case.smk"
include: "rules/R_filter_variants_special_variant_case.smk"
