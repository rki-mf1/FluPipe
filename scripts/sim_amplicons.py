#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, FG13, fuchss@rki.de)

VERSION = "0.0.9"
import os
import argparse
import re
from Bio import SeqIO
import sys
import random
import time
import gzip

def parse_args():
    parser = argparse.ArgumentParser(prog="sim_amplicon.py", description="generates simple and ideal amplicon data (paired-end). Six files are generated or overwritten: forward reads (<prefix>_R1.fq.gz), reverse reads (<prefix>_R2.fq.gz), primer clipped forward reads (<prefix>_clipped_R1.fq.gz), primer clipped reverse reads (<prefix>_clipped_R2.fq.gz), primer list in ptrimmer format (<prefix>.tsv), reference sequence with ambiguous call replaced by random explicit base call (<prefix>_reference.fasta)", )
    parser.add_argument('ref', metavar='FASTA_FILE', help="reference genome sequence in FASTA format", type=argparse.FileType('r'))
    parser.add_argument('--imin', metavar='INT', help="minimal insert length (default: 150)", type=int, default=150)
    parser.add_argument('--imax', metavar='INT', help="maximal insert length (default: 1200)", type=int, default=1200)
    parser.add_argument('--pmin', metavar='INT', help="maximal primer length (default: 18)", type=int, default=18)
    parser.add_argument('--pmax', metavar='INT', help="minimal primer length (default: 22)", type=int, default=23)
    parser.add_argument('--cmin', metavar='INT', help="minimal copies per amplicon (default: 25)", type=int, default=25)
    parser.add_argument('--cmax', metavar='INT', help="maximal copies per amplicon (default: 50)", type=int, default=50)
    parser.add_argument('--smin', metavar='INT', help="minimal sliding window (default: 10)", type=int, default=10)
    parser.add_argument('--smax', metavar='INT', help="maximal sliding window (default: 100)", type=int, default=100)
    parser.add_argument('--qmin', metavar='INT', help="PHRED score (default: 28)", type=int, default=28)
    parser.add_argument('--qmax', metavar='INT', help="PHRED score (default: 38)", type=int, default=38)
    parser.add_argument('--len', metavar='INT', help="read length (default: 150)", type=int, default=150)
    parser.add_argument('--debug', help="debug mode", action="store_true")
    parser.add_argument('--out', metavar='STR', help="file basenames (default: out). Files will be overwritten!", type=str, default="out")
    parser.add_argument('-a', help="allow non-unique R2 primer (faster)", action="store_true")
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    return parser.parse_args()

def fasta_to_dict(fname, to_strip=False):
    entries = {}
    for record in SeqIO.parse(fname, "fasta"):
        entries[record.description] = str(record.seq)
    if to_strip:
        entries = {x[0]: x[1].strip(to_strip) for x in entries.items()}
    return entries

def rev_compl(seq):
    return seq.lower().replace("a", "T").replace("t", "A").replace("g", "C").replace("c", "G")[::-1].upper()

def get_phred(score, base=33):
    return chr(base+score)

def get_qualstr(length, minqual, maxqual):
    quals = []
    for i in range(length):
        qual = random.randint(minqual, maxqual)
        quals.append(get_phred(qual))
    return "".join(quals)

def main():
    args = parse_args()
    allowed = {"A", "T", "G", "C"}
    refseq = "".join([x if x in allowed else list(allowed)[random.randint(0, len(allowed)-1)] for x in  list(str(SeqIO.read(args.ref, "fasta").seq).upper()) ])
    with open(args.out + "_reference.fasta", "w") as handle:
        handle.write(">reference\n" + refseq)
    rev_compl_refseq = rev_compl(refseq)
    seqlen = len(refseq)
    with open(args.out + ".tsv", "w") as handle:
        repeat = True
        r1 = []
        r2 = []
        r1_clipped = []
        r2_clipped = []
        primer_pairs = [("N", "N")]
        read = 0
        amplicon_start = 0
        while repeat:
            if seqlen - amplicon_start >= args.imax:
                insert_length = random.randint(args.imin, args.imax)
                primer_lengths = (random.randint(args.pmin, args.pmax), random.randint(args.pmin, args.pmax))
            else:
                primer_lengths = (random.randint(args.pmin, args.pmax), random.randint(args.pmin, args.pmax))
                insert_length = seqlen - amplicon_start - sum(primer_lengths) - 1
                repeat = False
            amplicon_end = amplicon_start + insert_length + sum(primer_lengths)
            full_fragment = refseq[amplicon_start:amplicon_end+1]
            if args.debug and amplicon_end - amplicon_start + 1 != len(full_fragment):
                exit(amplicon_end - amplicon_start + 1, len(full_fragment))

            #fwd read
            r1_read = full_fragment[:args.len]
            r1_primer = r1_read[:primer_lengths[0]]
            r1_read_qual = get_qualstr(args.len, args.qmin, args.qmax)
            r1_read_clipped = r1_read[primer_lengths[0]:]
            r1_read_clipped_qual = r1_read_qual[primer_lengths[0]:]

            #rev read
            r2_read = rev_compl(full_fragment[-args.len:])
            r2_primer = r2_read[:primer_lengths[1]]
            r2_read_qual = get_qualstr(args.len, args.qmin, args.qmax)
            r2_read_clipped = r2_read[primer_lengths[1]:]
            r2_read_clipped_qual = r2_read_qual[primer_lengths[1]:]

            copies = random.randint(args.cmin, args.cmax)
            handle.write(r1_primer + "\t" + r2_primer + "\t" + str(insert_length) + "\tcopies " + str(copies) + "\n")
            read += 1
            for i in range(copies):
                r1.append("@SEQ_ID_" + str(read) + "." + str(i+1))
                r1.append(r1_read)
                r1.append("+")
                r1.append(r1_read_qual)
                r2.append("@SEQ_ID_" + str(read) + "." + str(i+1))
                r2.append(r2_read)
                r2.append("+")
                r2.append(r2_read_qual)

                r1_clipped.append("@SEQ_ID_" + str(read) + "." + str(i+1))
                r1_clipped.append(r1_read_clipped)
                r1_clipped.append("+")
                r1_clipped.append(r1_read_clipped_qual)
                r2_clipped.append("@SEQ_ID_" + str(read) + "." + str(i+1))
                r2_clipped.append(r2_read_clipped)
                r2_clipped.append("+")
                r2_clipped.append(r2_read_clipped_qual)

            if args.debug:
                print("read pair:       " + str(read))
                print("ampl start:      " + str(amplicon_start))
                print("insert len:      " + str(insert_length))
                print("copies:          " + str(copies))
                print("r1 read:         " + str(r1_read))
                print("r1 clipped:      " + " " * primer_lengths[0] + r1_read_clipped)
                print("r1 primer:       " + r1_primer)
                print("r2 read:         " + str(r2_read))
                print("r2 clipped:      " + " " * primer_lengths[1] + r2_read_clipped)
                print("r2 primer:       " + r2_primer)
                print()

            sliding_window = random.randint(args.smin, args.smax)
            amplicon_start += sliding_window

    with gzip.open(args.out + '_R1.fq.gz', 'wt') as handle:
        handle.write("\n".join(r1))
    with gzip.open(args.out + '_R2.fq.gz', 'wt') as handle:
        handle.write("\n".join(r2))
    with gzip.open(args.out + '_clipped_R1.fq.gz', 'wt') as handle:
        handle.write("\n".join(r1_clipped))
    with gzip.open(args.out + '_clipped_R2.fq.gz', 'wt') as handle:
        handle.write("\n".join(r2_clipped))
if __name__ == "__main__":
	main()
