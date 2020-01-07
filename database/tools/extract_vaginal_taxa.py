#!/usr/bin/env python3
""" Goes through RefSeq's "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt" 
to extract sequence identiers for OptiVagDB. Taxa not found are flagged in a separate file for manual download.
"""

from collections import defaultdict
import argparse
import csv
import re

def read_bac_list(taxfile):
    wantedtaxa = set()
    wantedsamples = set()
    allwanted = set()
    with open(taxfile) as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        next (reader, None)
        for row in reader:
            allwanted.add(row[1])
            if(row[0] == "BioSample"):
                wantedsamples.add(row[1])
            else:
                wantedtaxa.add(row[1])
    return(wantedtaxa, wantedsamples, allwanted)


def sort_taxa(wantedtaxa, wantedsamples, allwanted, taxfile, foundfile, notfoundfile):
    bacfound = open(foundfile, "w+")
    notfound = open(notfoundfile, "w+")
    with open(taxfile) as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        for row in reader:
            if(row[0][0] != "#"):
                sample = row[2]
                taxon = row[7]
                parts = taxon.split(" ")
                subspecies = re.search("subsp", taxon)
                if subspecies:
                    binom = " ".join(parts[0:4])
                else:
                    binom = " ".join(parts[0:2])
                if(binom in wantedtaxa or sample in wantedsamples):
                    bacfound.write("\t".join(row) + "\n")
                    allwanted.discard(binom)
                    allwanted.discard(sample)
    for element in allwanted:
        notfound.write(element + "\n")
            

def main(assem_sum, taxlist, foundfile, notfoundfile):
    wantedtaxa, wantedsamples, allwanted = read_bac_list(taxlist)
    sort_taxa(wantedtaxa, wantedsamples, allwanted, assem_sum, foundfile, notfoundfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract the relevant rows from RefSeq:s assembly_summary.txt')
    parser.add_argument('--assem', help='path to assembly_summary.txt')
    parser.add_argument('--touse', help='TSV with taxa to use in the second column, comment on the first column')
    parser.add_argument('--found', help='subset of assembly_summary with taxa present in touse')
    parser.add_argument('--notfound', help='taxa in touse not found in assembly_summary')
    args = parser.parse_args()

    main(args.assem, args.touse, args.found, args.notfound)

