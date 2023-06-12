#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, FG13, fuchss@rki.de)

VERSION = "0.0.9"
import os
import argparse
import difflib
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(prog="compare_seqs.py", description="compares sequences in two FASTA files", )
    parser.add_argument('file1', metavar='FASTA_FILE', help="FASTA file 1", type=argparse.FileType('r'))
    parser.add_argument('file2', metavar='FASTA_FILE', help="FASTA file 2", type=argparse.FileType('r'))
    parser.add_argument('-c', help='case-sensitive comparison', action="store_true")
    parser.add_argument('-d', help='detailed ouput', action="store_true")
    parser.add_argument('--cutoff', metavar="INT", help='ignore given numer of chars at both sequence ends [0]', type=int, default=0)
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    return parser.parse_args()

def readFasta(fname):
    return str(SeqIO.read(fname, "fasta").seq)

def main():
    args = parse_args()
    seq1 = readFasta(args.file1)
    seq2 = readFasta(args.file2)
    if not args.c:
        seq1 = seq1.upper()
        seq2 = seq2.upper()
    if args.cutoff > 0:
        args.cutoff = abs(args.cutoff)
        seq1 = seq1[args.cutoff:-1*args.cutoff]
        seq2 = seq2[args.cutoff:-1*args.cutoff]
    if seq1 == seq2:
        print("sequences identical")
    else:
        if args.d:
            seq1 = [seq1[i:i+60] for i in range(0, len(seq1), 60)]
            seq2 = [seq2[i:i+60] for i in range(0, len(seq2), 60)]
            d = difflib.Differ()
            diff = d.compare(seq1, seq2)
            print('\n'.join(diff))
            print()
            print()
        print("Sequences are different!")
        exit(1)            

if __name__ == "__main__":
	main()
