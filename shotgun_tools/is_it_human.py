#!/usr/bin/env python3
""" Count human and not human reads """
__author__ = "Luisa W. Hugerth"
__date__ = "2020"

from sys import argv, exit
import csv
import argparse
import gzip
from Bio import SeqIO


def parse_args():
    """Count human and not human reads; slow, but ok memory usage"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", "-t", dest="class_tab", 
            help="Unoise mapping to human UC output ")
    parser.add_argument("fastqs", nargs="+", 
            help="fastq files to quantify human %")
    
    if len(argv) < 2:
        parser.print_help()
        exit()
            
        return parser.parse_args()


def sort_reads(class_tab):
    human = set()
    with open(class_tab) as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        for row in reader:
            if(row[0]=="H"):
                human.add(row[8])
    return(human)


def parse_fastq(infile, humanlist):
    human_counts = 0
    nothuman_counts = 0
    with gzip.open(infile, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if(record.id in humanlist):
                human_counts += 1
            else:
                nothuman_counts += 1
    return(human_counts, nothuman_counts)



def parse_fastqs(fastq_files, humanreads):
    humancounts = [None] * len(fastq_files)
    nothumancounts = [None] * len(fastq_files)
    for i in range(0, len(fastq_files)):
        myfile = fastq_files[i]
        h, n = parse_fastq(myfile, humanreads)
        humancounts[i] = h
        nothumancounts[i] = n
    return(humancounts, nothumancounts)


def main(class_tab, fastq_files):
    human = sort_reads(class_tab)
    Nhuman, Nnot = parse_fastqs(fastq_files, human)
    print("File\t" + "\t".join(fastq_files))
    print("Human_counts\t" + "\t".join([str(elem) for elem in Nhuman]))
    print("Not_human_counts\t" + "\t".join([str(elem) for elem in Nnot]))

if __name__ == "__main__":
    args = parse_args()
    main(args.class_tab, args.fastqs)
