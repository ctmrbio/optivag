#!/usr/bin/env python3
""" Goes through RefSeq's "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt" 
to extract sequence identiers for OptiVagDB. Taxa not found are flagged in a separate file for manual download.
"""
__author__ = "Luisa Hugerth"
__date__ = "2020"
__version__ = "0.2"

from collections import defaultdict
from sys import argv, exit
import argparse
import csv
import re


def read_bac_list(taxfile):
    wanted_taxa = set()
    wanted_samples = set()
    all_wanted = set()
    with open(taxfile) as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        next(reader, None)
        for row in reader:
            all_wanted.add(row[1])
            if row[0] == "BioSample":
                wanted_samples.add(row[1])
            else:
                wanted_taxa.add(row[1])
    return (wanted_taxa, wanted_samples, all_wanted)


def sort_taxa(wanted_taxa, wanted_samples, all_wanted, taxfile, foundfile, notfoundfile):
    with open(taxfile) as csvfile, open(foundfile, "w+") as bacfound, open(notfoundfile, "w+") as notfound:
        reader = csv.reader(csvfile, delimiter="\t")
        for row in reader:
            if row[0][0] != "#":
                sample = row[2]
                taxon = row[7]
                parts = taxon.split(" ")
                subspecies = re.search("subsp", taxon)
                genomospecies = re.search("genomosp", taxon)
                draft = re.search(" bacterium ", taxon)
                if subspecies:
                    binom = " ".join(parts[0:4])
                elif genomospecies or draft:
                    binom = " ".join(parts[0:3])
                else:
                    binom = " ".join(parts[0:2])
                if binom in wanted_taxa or sample in wanted_samples:
                    bacfound.write("\t".join(row) + "\n")
                    all_wanted.discard(binom)
                    all_wanted.discard(sample)
        for element in all_wanted:
            notfound.write("Not found\t" + element + "\n")


def main(assem_sum, taxlist, foundfile, notfoundfile):
    wanted_taxa, wanted_samples, allwanted_ = read_bac_list(taxlist)
    sort_taxa(wanted_taxa, wanted_samples, allwanted_, assem_sum, foundfile, notfoundfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=f"{__doc__} Copyright (c) {__date__} {__author__}."
    )
    parser.add_argument("--assem", help="path to assembly_summary.txt")
    parser.add_argument("--to-use",
        help="TSV with taxa to use in the second column, comment on the first column.",
    )
    parser.add_argument("--found", 
        help="Subset of assembly_summary with taxa present in touse.",
    )
    parser.add_argument("--not-found", 
        help="Taxa in to-use not found in assembly_summary.",
    )
    
    if len(argv) < 2:
        parser.print_help()
        exit()

    args = parser.parse_args()

    main(args.assem, args.to_use, args.found, args.not_found)
