#!/usr/bin/env python3
""" Filters a subset of RefSeq's "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt" 
to keep a maximum number of genomes per species, privileging reference genomes.
"""

from sys import argv, exit
from collections import defaultdict
import argparse
import csv
import re


def parse_infile(infile):
    """Parse a file with the following structure:
	GCF_000067165.1 PRJNA224116 SAMEA3138271		representative genome   448385  56  Sorangium cellulosum So ce56	strain=So ce 56	 latest  Complete Genome Major   Full	2007/11/27  ASM6716v1   Bielefeld Univ  GCA_000067165.1 identical   ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/067/165/GCF_000067165.1_ASM6716v1
    """
    with open(infile) as csvfile:
        eachspecies = set()
        ref_genomes = defaultdict(lambda: defaultdict(str))
        repr_genomes = defaultdict(lambda: defaultdict(str))
        complete_genomes = defaultdict(lambda: defaultdict(str))
        chromosomes = defaultdict(lambda: defaultdict(str))
        scaffolds = defaultdict(lambda: defaultdict(str))
        tax_dict = defaultdict(lambda: defaultdict(str))

        reader = csv.reader(csvfile, delimiter="\t")
        for row in reader:
            if row[0][0] != "#":
                genome_type = row[4]
                tax_id = row[5]
                tax_name = row[7]
                strain = row[8]
                if "=" in strain:
                    strain = strain.split("=")[1]
                elif strain == "":
                    strain = row[9]
                seq_type = row[11]
                url = row[19]
                taxa = tax_name.split(" ")

                if len(taxa) > 2 and taxa[2] == "subsp.":
                    taxon = " ".join(taxa[0:4])
                elif len(taxa) > 2 and taxa[1] == "genomosp.":
                    taxon = " ".join(taxa[0:3])
                else:
                    taxon = " ".join(taxa[0:2])
                eachspecies.add(taxon)
                strain = " ".join([taxon, strain])

                tax_dict[taxon][strain] = tax_id

                if genome_type == "reference genome":
                    ref_genomes[taxon][strain] = url
                elif genome_type == "representative genome":
                    repr_genomes[taxon][strain] = url
                elif seq_type == "Complete genome":
                    complete_genomes[taxon][strain] = url
                elif seq_type == "Chromosome":
                    chromosomes[taxon][strain] = url
                elif seq_type == "Scaffold":
                    scaffolds[taxon][strain] = url
        return (
            eachspecies,
            ref_genomes,
            repr_genomes,
            chromosomes,
            scaffolds,
            tax_dict,
        )


def print_reference(outgenomes, outproteins, outtaxa, strain_name, url, tax_id):
    with open(outgenomes, "a+") as g, open(outproteins, "a+") as p, open(
        outtaxa, "a+"
    ) as t:
        parts = url.split("/")
        identifier = parts[-1]
        g.write("wget " + url + "/" + identifier + "_genomic.fna.gz\n")
        g.write(
            "mv " + identifier + "_genomic.fna.gz " + strain_name + "_genomic.fna.gz\n"
        )
        p.write("wget " + url + "/" + identifier + "_protein.faa.gz\n")
        p.write(
            "mv " + identifier + "_protein.faa.gz " + strain_name + "_protein.faa.gz\n"
        )
        t.write(strain_name + "_genomic.fna\t" + tax_id + "\n")


def select_genomes(
    eachspecies,
    ref_genomes,
    repr_genomes,
    chromosomes,
    scaffolds,
    taxa_id,
    maxgenomes,
    outgenomes,
    outproteins,
    outtaxa,
):
    for species in eachspecies:
        genomes_found = 0
        if species in ref_genomes:
            strains = ref_genomes[species].keys()
            for strain in strains:
                if genomes_found < maxgenomes:
                    genomes_found += 1
                    tax_id = taxa_id[species][strain]
                    strain_name = strain.replace(" ", "_")
                    url = ref_genomes[species][strain]
                    print_reference(
                        outgenomes, outproteins, outtaxa, strain_name, url, tax_id
                    )
        if species in repr_genomes and genomes_found < maxgenomes:
            strains = repr_genomes[species].keys()
            for strain in strains:
                if genomes_found < maxgenomes:
                    genomes_found += 1
                    tax_id = taxa_id[species][strain]
                    strain_name = strain.replace(" ", "_")
                    url = repr_genomes[species][strain]
                    print_reference(
                        outgenomes, outproteins, outtaxa, strain_name, url, tax_id
                    )
        if species in chromosomes and genomes_found < maxgenomes:
            strains = chromosomes[species].keys()
            for strain in strains:
                if genomes_found < maxgenomes:
                    genomes_found += 1
                    tax_id = taxa_id[species][strain]
                    strain_name = strain.replace(" ", "_")
                    url = chromosomes[species][strain]
                    print_reference(
                        outgenomes, outproteins, outtaxa, strain_name, url, tax_id
                    )
        if species in scaffolds and genomes_found < maxgenomes:
            strains = scaffolds[species].keys()
            for strain in strains:
                if genomes_found < maxgenomes:
                    genomes_found += 1
                    tax_id = taxa_id[species][strain]
                    strain_name = strain.replace(" ", "_")
                    url = scaffolds[species][strain]
                    print_reference(
                        outgenomes, outproteins, outtaxa, strain_name, url, tax_id
                    )


def main(infile, outgenomes, outproteins, outtaxa, maxgenomes):
    eachspecies, ref_genomes, repr_genomes, chromosomes, scaffolds, taxa_id = parse_infile(
        infile
    )
    select_genomes(
        eachspecies,
        ref_genomes,
        repr_genomes,
        chromosomes,
        scaffolds,
        taxa_id,
        maxgenomes,
        outgenomes,
        outproteins,
        outtaxa,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract the relevant rows from RefSeq:s assembly_summary.txt"
    )
    parser.add_argument("-i", "--infile", 
        help="Subset of assembly_summary with selected taxa.",
    )
    parser.add_argument("-g", "--outgenomes",
        help="sh file to be created for getting genome information.",
    )
    parser.add_argument("-p", "--outproteins",
        help="sh file to be created for getting protein information.",
    )
    parser.add_argument("-t", "--outtaxa", 
        help="tsv file to be created with the taxon id for each file."
    )
    parser.add_argument("-n", "--maxgenomes",
        type=int,
        default=10,
        help="Maximum number of genomes per sub/species",
    )

    if len(argv) < 2:
        parser.print_help()
        exit()

    args = parser.parse_args()

    main(args.infile, args.outgenomes, args.outproteins, args.outtaxa, args.maxgenomes)
