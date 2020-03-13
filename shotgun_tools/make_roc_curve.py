#!/usr/bin/env python3
""" Print numbers of true and false positives and negatives """
__author__ = "Luisa W. Hugerth"
__date__ = "2020"

from sys import argv, exit
import csv
import argparse
import gzip
from Bio import SeqIO
import re


def parse_args():
    """Scales linearly with input file size, be careful"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", "-t", dest="class_tab", 
            help="Unoise UC or SAM table mapping to human")
    parser.add_argument('--unoise', dest='unoise', action='store_true',
            help="mapping table is in Unoise UC format")
    parser.add_argument('--sam', dest='unoise', action='store_false',
            help="mapping table is in SAM format (default)")
    parser.set_defaults(unoise=False)
    parser.add_argument("fastqs", nargs="+", 
            help="fastq files to quantify human and non-human reads")
    
    if len(argv) < 3:
        parser.print_help()
        exit()
        
    return parser.parse_args()


def sort_reads(class_tab, unoise):
    human = set()
    nonhuman = set()
    with open(class_tab) as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        if(unoise == True):
            print("Parsing Unoise file...")
            for row in reader:
                if(row[0]=="H"):
                    human.add(row[8])
                else:
                    nonhuman.add(row[8])
        else:
            print("Parsing SAM file...")
            for row in reader:
                if(row[0][0]!="@"):
                    if(row[5]=='*'):
                        nonhuman.add(row[0])
                    else:
                        human.add(row[0])
    # Sanity check in case Usearch did something weird:                
    nonhuman = nonhuman.difference(human)
    return(human, nonhuman)


def parse_fastq(infile, humanlist, buglist):
    with gzip.open(infile, "rt") as handle:
        record_list = list(SeqIO.to_dict(SeqIO.parse(handle, "fastq")))
        record_list = [record.replace('/1', '') for record in record_list]
        record_list = [record.replace('/2', '') for record in record_list]
        record_set = set(record_list)
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
        print("Parsing", myfile)
        tp, fp, tn, fn = parse_fastq(myfile, humanreads, bugreads)
        tpr[i] = tp
        fpr[i] = fp
        tnr[i] = tn
        fnr[i] = fn
    return(tpr, fpr, tnr, fnr)


def main(class_tab, unoise, fastq_files):
    human, nonhuman = sort_reads(class_tab, unoise)
    TPR, FPR, TNR, FNR = parse_fastqs(fastq_files, human, nonhuman)
    print("File\t" + "\t".join(fastq_files))
    print("True positives\t" + "\t".join([str(elem) for elem in TPR]))
    print("False positives\t" + "\t".join([str(elem) for elem in FPR]))
    print("True negatives\t" + "\t".join([str(elem) for elem in TNR]))
    print("False negatives\t" + "\t".join([str(elem) for elem in FNR]))

if __name__ == "__main__":
    args = parse_args()
    print(args)
    main(args.class_tab, args.unoise, args.fastqs)
