#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF-1, fuchss@rki.de)

VERSION = "0.0.9"
import os
import argparse
import re
import sys
import vcf

def parse_args(CMD=None):
    parser = argparse.ArgumentParser(prog="rename_in_gff3.py", description="changes genotype in VCFs", )
    parser.add_argument('--vcf', metavar="FILE", help="vcf file", type=str, required=True)
    parser.add_argument('-o', help="output gz file (will be overwritten!)", type=str, required=True)
    return parser.parse_args(CMD)

def get_cmd_from_snake(snakemake):
    cmd = ["-o", snakemake.output[0],
	       "--vcf", snakemake.input[0]]
    return cmd

# open file handles considering compression state
def process(in_fname, out_fname):
    vcf_reader = vcf.Reader(filename=in_fname)
    vcf_writer = vcf.Writer(open(out_fname, 'w'), vcf_reader)
    for record in vcf_reader:
         if len(record.REF) > len(record.ALT):
              vcf_writer.write_record(record)
         #if len(record.REF) != 1:
         #     record.ALT = ",".join([x + "?"*(len(record.REF)-len(x)) for x in record.ALT.split(",")])
         #vcf_writer.write_record(record)

def main(CMD=None):
        args = parse_args(CMD)
        process(args.vcf, args.o)

if __name__ == "__main__":
        CMD = None
        if "snakemake" in globals():
            CMD = get_cmd_from_snake(snakemake)
        main(CMD=CMD)
