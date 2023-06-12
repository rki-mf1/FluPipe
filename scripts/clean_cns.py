#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF-1, fuchss@rki.de)

VERSION = "0.0.9"
import os
import argparse
from Bio import SeqIO

def parse_args(CMD=None):
    parser = argparse.ArgumentParser(prog="clean_cns.py", description="clean consensus sequences", )
    parser.add_argument('--fasta', metavar="FILE", help="fasta file", type=str, required=True)
    parser.add_argument('-o', metavar="FILE", help="output file (will be overwritten!)", type=str, required=True)
    parser.add_argument('--mask', help="mask ambiguities using N *(default n)", action="store_true")
    return parser.parse_args(CMD)

def get_cmd_from_snake(snakemake):
    cmd = ["-o", snakemake.output[0],
	       "--fasta", snakemake.input[0]]
    if snakemake.params['mask']:
        cmd.append("--mask")
    return cmd

# open file handles considering compression state
def process(in_fname, out_fname, mask=False):
	entries = []
	for record in SeqIO.parse(in_fname, "fasta"):
		seq = str(record.seq).replace("?", "")
		if mask:
			for char in "RYSWKMBDHV":
				seq = "".join([x if x not in set("RYSWKMBDHV") else "N" for x in seq])
		n = 80
		seq = "".join([seq[x : x + n] for x in range(0, len(seq), n)])
		entries.append(">" + record.description + "\n" + seq)

	with open(out_fname, "w") as handle:
		handle.write("\n".join(entries))

def main(CMD=None):
        args = parse_args(CMD)
        process(args.fasta, args.o, args.mask)

if __name__ == "__main__":
        CMD = None
        if "snakemake" in globals():
            CMD = get_cmd_from_snake(snakemake)
        main(CMD=CMD)
