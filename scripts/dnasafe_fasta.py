#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, FG13, fuchss@rki.de)

VERSION = "0.0.9"
import os
import sys
import argparse
import textwrap
from Bio import SeqIO

def parse_args(CMD=None):
	parser = argparse.ArgumentParser(prog="dna_safe.py", description="writes a dna-safe fasta file", )
	parser.add_argument('-f', '--fasta', metavar='FASTA_FILE', help="fasta file to check", type=str, required=True)
	parser.add_argument('-o', '--out', metavar='FASTA_FILE', help="dna-safe fasta file to write (exisiting file will be overwritten!)", type=str, required=True)
	parser.add_argument('-l', '--len', metavar='INT', help="sequence line length (set 0 to avoid new lines, default=60)", type=int, default=60)
	parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
	return parser.parse_args(CMD)


def get_cmd_from_snake(snakemake):
	cmd = ["-o", snakemake.output[0],
		  "--fasta", snakemake.input[0]]
	cmd = [str(arg) for arg in cmd]
	return cmd

def main():
	valid_letters = "ATGCURYSWKMBDHVN.-"

	#args
	CMD = None
	if "snakemake" in globals():
		CMD = get_cmd_from_snake(snakemake)
	args = parse_args(CMD=CMD)

	#check file
	if not os.path.isfile(args.fasta):
		sys.exit("error: the input fasta file does not exist")

	# process file
	with open(args.out, "w") as handle:
		for record in SeqIO.parse(args.fasta, "fasta"):
			header = str(record.description)
			#remove colon and spaces in header for subsequent tools to work
			header = header.replace(":", "_").replace(" ", "_").replace("(", "_").replace(")", "_").replace("'", "_").replace(".", "_").replace("[", "_").replace("]", "_")
			seq = str(record.seq).upper()
			if not all(i in valid_letters for i in seq):
				sys.exit("error: sequence '" + header + "' in " + args.fasta + " conatins non-IUPAC characters.")
			#replace any non standard characters - otherwise problems with lofreq and bcltools consensus
			seq = seq.replace(".", "").replace("-", "").replace("U", "T").replace("W", "A").replace("S", "C").replace("M", "A").replace("K", "G").replace("R", "A").replace("Y", "C").replace("B", "C").replace("D", "A").replace("H", "A").replace("V", "A").replace("N", "A")
			if args.len > 0:
				seq = textwrap.fill(seq, width=args.len)
			handle.write(">" + header + "\n" + seq + "\n")

if __name__ == "__main__":
	main()
