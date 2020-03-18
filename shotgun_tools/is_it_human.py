#!/usr/bin/env python3
""" Count human and not human reads """
__author__ = "Luisa W. Hugerth"
__date__ = "2020"

from sys import argv, exit
from collections import defaultdict
import csv
import argparse
import gzip
from Bio import SeqIO


def parse_args():
    """Count human and not human reads; slow, but ok memory usage"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", "-t", 
            dest="class_tab", 
            help="Unoise mapping to human UC output."
    )
    parser.add_argument("fastqs", nargs="+", help="fastq files to quantify human %.")

    if len(argv) < 2:
        parser.print_help()
        exit()

    return parser.parse_args()


def sort_reads(class_tab):
    human_reads = set()
    with open(class_tab) as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        for row in reader:
            if row[0] == "H":
                human_reads.add(row[8])
    return human_reads


def parse_fastq(infile, humanlist):
    counts = defaultdict(int)
    with gzip.open(infile, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if record.id in humanlist:
                counts["human"] += 1
            else:
                counts["nothuman"] += 1
    return counts


def get_total_counts(fastq_files, humanreads):
    total_counts = defaultdict(lambda: defaultdict(list))
    for fastq in fastq_files:
        counts = parse_fastq(fastq, humanreads)
        total_counts[fastq].update(counts)
    return total_counts


def main(class_tab, fastq_files):
    human_reads = sort_reads(class_tab)
    total_counts = get_total_counts(fastq_files, human_reads)
    print("File\tHuman\tNot_human")
    for fastq, counts in total_counts.items():
        print(f"{fastq}\t{counts['human']}\t{counts['nothuman']}")


if __name__ == "__main__":
    args = parse_args()
    main(args.class_tab, args.fastqs)
