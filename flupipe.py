#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF-1, fuchss@rki.de)
#adapted winterk


__version__ = "0.0.10"

import os
import argparse
import re
import sys
import yaml
from pathlib import Path
import subprocess

default_regex = "([^" + os.path.sep + "]+)_R[12]_.+"

def parse_args():
	parser = argparse.ArgumentParser(prog="rsvpipe.py", description="an automated pipeline for reference-based RSV genome assembly")

	# data input
	group_name = "data input"
	group_data = parser.add_argument_group(group_name + "\n" + ":" * (len(group_name)-1))
	group_data.add_argument('-f', '--file', metavar="FILE", help="fastq files with paired-end read data (forward and reverse reads are expected in individual files)", type=str, nargs='+', default=None)
	group_data.add_argument('-d', '--dir', metavar="DIR", help="directory that contains fatsq files with paired-end read data", type=str, default=None)
	group_data.add_argument('-c', '--configfile', metavar="FILE", help="use existing contig file (no other option allowed)", type=str, default=None)
	group_data.add_argument('-t', '--tag', metavar="STR", help="file name tag labelling forward and reverse read files (per default _R1_ and _R2_ )", type=str, nargs="+", default=['_R1_', '_R2_'])
	group_data.add_argument('-r', '--recursive', help="include files in subfolders of the directory ()", action="store_true")
	group_data.add_argument('-e', '--ext', metavar="STR", help="allowed file extension (only effective when --dir is used; per default: .fq  .fastq  .fq.gz .fastq.gz )", type=str, nargs="+", default=['.fq', '.fastq', '.fq.gz', '.fastq.gz'])


	# read trimming and filtering
	group_name = "read trimming and filtering"
	group_trimming = parser.add_argument_group(group_name + "\n" + ":" * (len(group_name)-1))
	group_trimming.add_argument('-p', '--phred', metavar="INT", help="minimal averaged PHRED (default: 20)", type=int, default=20)
	group_trimming.add_argument('-l', '--length', metavar="INT", help="minimal read length (default: 50)", type=int, default=50)
	group_trimming.add_argument('-a', '--adapter', metavar="FILE", help="fasta file with adapter sequences to clip (default: no clipping)", type=str, default=None)
	group_trimming.add_argument('-k', '--krakendb', metavar="DIR", help="kraken2 database for taxonomic read filtering (default: no filtering)", type=str, default=None)
	group_trimming.add_argument('-s', '--segmentdb', metavar="DIR", help="folder containing fasta files with the influenza reference sequences per segment - needed for automatic reference detection", type=str, default=None)
	#group_trimming.add_argument('-x', '--taxid', metavar="STR", help="target taxonomic identifier (default: 12814)", type=str, default="12814")
	group_trimming.add_argument('-m', '--primer', metavar="FILE", help="bedpe file with information about amplicon primers to be clipped (default: no clipping)", type=str, default=None)
	group_trimming.add_argument('-u', '--dedup', help="remove read duplicates (default: no deduplication)", action="store_true")

	# mapping
	group_name = "read mapping"
	group_mapping = parser.add_argument_group(group_name + "\n" + ":" * (len(group_name)-1))
	group_mapping.add_argument('--ref', metavar="FILE", help="fasta file with reference sequence (default reference: automatic reference selection from GISAID DB)", type=str, default=None)

	# variant calling
	group_name = "variant calling"
	group_calling = parser.add_argument_group(group_name + "\n" + ":" * (len(group_name)-1))
	group_calling.add_argument('--var_mincov', metavar="INT", help="minimum number of reads required for variant calling (default: 50)", type=int, default=50)
	group_calling.add_argument('--mincount', metavar="INT", help="minimum number of supporting reads required to call a variant (default: 10)", type=int, default=10)
	group_calling.add_argument('--minfrac', metavar="FLOAT", help="minimum percentage of supporting reads required to call a variant (default: 0.1)", type=int, default=0.1)

	# variant masking
	group_name = "variant filtering"
	group_filter = parser.add_argument_group(group_name + "\n" + ":" * (len(group_name)-1))
	group_filter.add_argument('--mqm', metavar="INT", help="minimal mean mapping quality of observed alternate alleles (default: 40)", type=float, default=40)
	group_filter.add_argument('--sap', metavar="INT", help="minimal strand balance probability for the alternate allele (default: 60)", type=int, default=60)
	group_filter.add_argument('--qual', metavar="INT", help="minimal variant call quality (default: 10)", type=int, default=10)

	# consensus building
	group_name = "consensus building"
	group_masker = parser.add_argument_group(group_name + "\n" + ":" * (len(group_name)-1))
	group_masker.add_argument('--cns_mincov', metavar="INT", help="minimum number of reads to call reference or alternate alleles (default: 50)", type=float, default=50)
	group_masker.add_argument('--gt', metavar="FLOAT", help="minimum percentage of reads to avoid ambiguous genotypes (genotype adjustment). (default: 0.9)", type=float, default=0.9)

	# general
	group_name = "general"
	group_out = parser.add_argument_group(group_name + "\n" + ":" * (len(group_name)-1))
	group_out.add_argument('-o', '--out', metavar="STR", help="output directory (default: current working directory)", type=str, default=".")
	group_out.add_argument('--threads', metavar="INT", help="number of threads to use (default: 1)", type=int, default=1)
	group_out.add_argument('--run_id', metavar="STR", help="run ID (default: None)", type=str, default=None)
	group_out.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)

	return parser.parse_args()


def main():
	config = {}
	samples = {}
	args = parse_args()

	# data file input checks
	i = sum(x is not None for x in [args.file, args.dir, args.configfile])
	if i > 1:
		sys.exit("input error: options --file, --dir, and --configfile are mutually exclusive.")
	elif i == 0:
		sys.exit("input error: no file, directory or existing config file has been set.")

	if len(args.tag) != 2:
		sys.exit("input error: there must be exactly two data file name tags defined.")

	if args.file:
		for fname in args.file:
			if not os.path.isfile(fname):
				sys.exit("input error: data file does not exist (" + os.path.basename(fname) + ").")
		fnames = args.file
	elif args.dir:
		if not os.path.isdir(args.dir):
			sys.exit("input error: data directory does not exist (" + args.dir + ").")
		p = Path(args.dir).rglob('*') if args.recursive else Path(args.dir).glob('*')
		fnames = [str(x.resolve()) for x in p if x.name.endswith(tuple(args.ext)) ]
	else:
		with open(args.configfile, "r") as handle:
			if not os.path.isfile(args.configfile):
				sys.exit("input error: config file does not exist (" + os.path.basename(args.configfile) + ").")
			try:
				config = yaml.safe_load(handle)
			except:
				sys.exit("input error: config file cannot processed.")

			if not "samples" in config:
				sys.exit("input error: missing samples definition in contig file")

			with open(config['samples'], "r") as handle:
				samples = yaml.safe_load(handle)

	if len(fnames)%2 != 0:
		sys.exit("input error: paired-end data expected (even number of input files mandatory)")
	elif len(fnames) == 0:
		sys.exit("input error: no data file found")

	if args.dir or args.file:
		fnames.sort()
		args.ext.sort(key = len, reverse=True)
		for i in range(0, len(fnames), 2):
			f1 = os.path.basename(fnames[i])
			f2 = os.path.basename(fnames[i+1])
			if f1.replace(args.tag[0], args.tag[1]) != f2 or f1 != f2.replace(args.tag[1], args.tag[0]):
				sys.exit("input error: file name pairing failed for data file pair " + os.path.basename(args.file[i]) + " and " + os.path.basename(args.file[i+1]))
			sample_id = f1.replace(args.tag[0], args.tag[0][-1])
			for ext in args.ext:
				if sample_id.endswith(ext):
					sample_id = sample_id[:-len(ext)]
					break
			if sample_id in samples:
				sys.exit("input error: sample id collision for different file pairs (id: " + sample_id + ").")
			samples[sample_id] = {}
			samples[sample_id]['altid'] = 'Sample_' + str(int(i/2))
			samples[sample_id]['read1'] = fnames[i]
			samples[sample_id]['read2'] = fnames[i+1]


		# trimming
		if args.phred < 0:
			 sys.exit("Input error: value of --phred cannot be smaller than 0.")

		if args.length < 0:
			 sys.exit("Input error: value of --lengths cannot be smaller than 0.")

		if args.adapter and not os.path.isfile(args.adapter):
			 sys.exit("Input error: value of --adapter is not a valid file.")

		if args.krakendb and not os.path.isdir(args.krakendb):
			 sys.exit("Input error: value of --krakendb is not a valid directory.")

		if args.primer and not os.path.isfile(args.primer):
			 sys.exit("Input error: value of --primer is not a valid file.")

		# mapping
		if args.ref and not os.path.isfile(args.ref):
			 sys.exit("Input error: value of --ref is not a valid file.")

		# variant calling
		if args.var_mincov < 0:
			 sys.exit("Input error: value of --var_mincov cannot be smaller than 0.")

		if args.mincount < 0:
			 sys.exit("Input error: value of --mincount cannot be smaller than 0.")

		if args.minfrac < 0 or args.minfrac > 1:
			 sys.exit("Input error: value of --minfrac has to be between 0 and 1.")

		# variant filter
		if args.mqm < 0:
			 sys.exit("Input error: value of --mqm cannot be smaller than 0.")
		if args.sap < 0:
			 sys.exit("Input error: value of --sap cannot be smaller than 0.")

		# consensus building
		if args.cns_mincov < 0:
			 sys.exit("Input error: value of --cns_mincov cannot be smaller than 0.")

		if args.gt < 0 or args.gt > 1:
			 sys.exit("Input error: value of --gt has to be between 0 and 1.")

		# assigning config keys
		args.out = os.path.abspath(args.out)
		config['output'] = args.out
		config['samples'] = os.path.join(args.out, "results", "samples.yaml")
		config['reference'] = os.path.abspath(args.ref) if args.ref else None
		config['krakenDb'] = os.path.abspath(args.krakendb) if args.krakendb else None
		config['segmentdb'] = os.path.abspath(args.segmentdb) if args.segmentdb else None
		#config['krakenTaxID'] = args.taxid
		config['adapter'] = os.path.abspath(args.adapter) if args.adapter else None
		config['primer'] = os.path.abspath(args.primer) if args.primer else None
		config['pcr_dedup'] = args.dedup
		config['var_call_cov'] = args.var_mincov
		config['var_call_count'] = args.mincount
		config['var_call_frac'] = args.minfrac
		config['var_filter_mqm'] = args.mqm
		config['var_filter_sap'] = args.sap
		config['var_filter_qual'] = args.qual
		config['read_filter_qual'] = args.phred
		config['read_filter_len'] = args.length
		config['cns_min_cov'] = args.cns_mincov
		config['run_id'] = args.run_id

	configfile = os.path.join(args.out, "results", "config.yaml")
	os.makedirs(os.path.join(args.out, "results"), exist_ok=True)

	with open(configfile, "w") as handle:
		yaml.dump(config, handle)

	with open(config['samples'], "w") as handle:
		yaml.dump(samples, handle)

	#run snakemake
	path = os.path.dirname(os.path.abspath(__file__))
	snakefile = os.path.join(path, "flupipe.Snakefile")
	cmd = ["snakemake", "-s", snakefile, "--configfile", configfile, "--cores",  str(args.threads), "--use-conda" , "--conda-frontend", "mamba", "--rerun-incomplete",  "--nolock"]
	#cmd = ["snakemake", "--dag -s", snakefile, "--configfile", configfile, "--cores",  str(args.threads), " | dot -Tsvg > dag.svg"]
	print("executing:", " ".join(cmd))
	subprocess.Popen(cmd).wait()

if __name__ == "__main__":
		main()
