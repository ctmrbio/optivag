#!/usr/bin/env python3
""" Print numbers of true and false positives and negatives """
__author__ = "Luisa W. Hugerth"
__date__ = "2020"

from sys import argv, exit
import csv
import argparse
import gzip
from Bio import SeqIO


def parse_args():
    """Scales nearly with input file size, be careful"""
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
    nonhuman = set()
    with open(class_tab) as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        for row in reader:
            if(row[0]=="H"):
                human.add(row[8])
            else:
                nonhuman.add(row[8])
    # Sanity check in case Usearch did something weird:                
    nonhuman = nonhuman.difference(human)
    return(human, nonhuman)


def parse_fastq(infile, humanlist, buglist):
    with gzip.open(infile, "rt") as handle:
        record_set = set(SeqIO.to_dict(SeqIO.parse(handle, "fastq")))
        #True positives: human reads NOT in record
        tpr = len(humanlist.difference(record_set))
        #False positives: bacterial reads NOT in record
        fpr = len(buglist.difference(record_set))
        #True negatives: bacterial reads in record
        tnr = len(buglist & record_set)
        #False negatives: human reads in record
        fnr = len(humanlist & record_set)
    return(tpr, fpr, tnr, fnr)



def parse_fastqs(fastq_files, humanreads, bugreads):
    tpr = [None] * len(fastq_files)
    fpr = [None] * len(fastq_files)
    tnr = [None] * len(fastq_files)
    fnr = [None] * len(fastq_files)
    for i in range(0, len(fastq_files)):
        myfile = fastq_files[i]
        tp, fp, tn, fn = parse_fastq(myfile, humanreads, bugreads)
        tpr[i] = tp
        fpr[i] = fp
        tnr[i] = tn
        fnr[i] = fn
    return(tpr, fpr, tnr, fnr)


def main(class_tab, fastq_files):
    human, nonhuman = sort_reads(class_tab)
    TPR, FPR, TNR, FNR = parse_fastqs(fastq_files, human, nonhuman)
    print("File\t" + "\t".join(fastq_files))
    print("True positives\t" + "\t".join([str(elem) for elem in TPR]))
    print("False positives\t" + "\t".join([str(elem) for elem in FPR]))
    print("True negatives\t" + "\t".join([str(elem) for elem in TNR]))
    print("False negatives\t" + "\t".join([str(elem) for elem in FNR]))

if __name__ == "__main__":
    args = parse_args()
    main(args.class_tab, args.fastqs)
