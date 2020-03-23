#!/usr/bin/env python3
""" Print numbers of true and false positives and negatives.
Scales linearly with input file size, be careful."""
__author__ = "Luisa W. Hugerth"
__date__ = "2020"
__version__ = "0.2"

from sys import argv, exit, stderr
from collections import defaultdict
import csv
import argparse
import gzip
from Bio import SeqIO
import re


def parse_args():
    desc = f"{__doc__} Copyright (c) {__date__} {__author__}. Version v{__version__}."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--table", "-t",
        dest="class_tab",
        help="Unoise UC or SAM table mapping to human.",
    )
    parser.add_argument("--unoise",
        dest="unoise",
        action="store_true",
        help="Mapping table is in Unoise UC format.",
    )
    parser.add_argument("--sam",
        dest="unoise",
        action="store_false",
        help="Mapping table is in SAM format (default).",
    )
    parser.set_defaults(unoise=False)
    parser.add_argument("fastqs", 
        nargs="+", 
        help="Fastq files to quantify human and non-human reads."
    )

    if len(argv) < 3:
        parser.print_help()
        exit()

    return parser.parse_args()


def sort_reads(class_tab, unoise):
    human = set()
    nonhuman = set()
    with open(class_tab) as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        if unoise:
            print("Parsing Unoise file...", file=stderr)
            for row in reader:
                if row[0] == "H":
                    human.add(row[8])
                else:
                    nonhuman.add(row[8])
        else:
            print("Parsing SAM file...", file=stderr)
            for row in reader:
                if row[0][0] != "@":
                    if row[5] == "*":
                        nonhuman.add(row[0])
                    else:
                        human.add(row[0])
    # Sanity check in case Usearch did something weird:
    nonhuman = nonhuman.difference(human)
    return (human, nonhuman)


def parse_fastq(infile, humanlist, buglist):
    print("Parsing", infile, file=stderr)
    counts = {}
    with gzip.open(infile, "rt") as fastq_file:
        record_list = list(SeqIO.to_dict(SeqIO.parse(fastq_file, "fastq")))
        record_list = [record.replace("/1", "") for record in record_list]
        record_list = [record.replace("/2", "") for record in record_list]
        record_set = set(record_list)
        counts["tp"] = len(humanlist.difference(record_set))   # True positives: human reads NOT in record
        counts["fp"] = len(buglist.difference(record_set))     # False positives: bacterial reads NOT in record
        counts["tn"] = len(buglist.intersection(record_set))   # True negatives: bacterial reads in record
        counts["fn"] = len(humanlist.intersection(record_set)) # False negatives: human reads in record
    return counts


def parse_fastqs(fastq_files, humanreads, bugreads):
    pos_neg_counts = defaultdict(lambda: defaultdict(list))
    for fastq in fastq_files:
        counts = parse_fastq(fastq, humanreads, bugreads)
        pos_neg_counts[fastq].update(counts)
    return pos_neg_counts


def main(class_tab, unoise, fastq_files):
    human, nonhuman = sort_reads(class_tab, unoise)
    pos_neg_counts = parse_fastqs(fastq_files, human, nonhuman)
    print("File\tTP\tFP\tTN\tFN")
    for fastq, counts in pos_neg_counts.items():
        print(f"{fastq}\t{counts['tp']}\t{counts['fp']}\t{counts['tn']}\t{counts['fn']}")


if __name__ == "__main__":
    args = parse_args()
    main(args.class_tab, args.unoise, args.fastqs)
