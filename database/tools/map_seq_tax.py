#!/usr/bin/env python3
""" Uses the assembly summary and the list of downloaded genomes to create
a taxonomy file compatible with Bracken.
"""

from sys import argv, exit
import argparse
import csv
import glob
import os


def read_taxa(taxa):
    id_dict = dict()
    with open(taxa) as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        for row in reader:
            id_dict[row[0]] = row[1]
    return id_dict


def make_list(id_dict, genomes, outfile):
    with open(outfile, "w+") as o:
        for filename in glob.glob(os.path.join(genomes, "*.fna")):
            genome = os.path.basename(filename)
            tax_id = id_dict[genome]
            with open(filename, "r") as csvfile:
                reader = csv.reader(csvfile, delimiter=" ")
                for row in reader:
                    if row[0][0] == ">":
                        seq_id = row[0][1:]
                        o.write(seq_id + "\t" + tax_id + "\n")


def main(taxa, genomes, outfile):
    id_dict = read_taxa(taxa)
    make_list(id_dict, genomes, outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Combines sequence headers with the appropriate taxonomy code"
    )
    parser.add_argument("-t", "--taxa", 
        help="TSV file connecting each genome to its NCBI tax_id."
    )
    parser.add_argument("-g", "--genomes", 
        help="Folder containing each genome in the database."
    )
    parser.add_argument("-o", "--out",
        nargs="?",
        default="seqid2taxid.map",
        help="Path to outfile (default seqid2taxid.map).",
    )
    
    if len(argv) < 2:
        parser.print_help()
        exit()

    args = parser.parse_args()

    main(args.taxa, args.genomes, args.out)
