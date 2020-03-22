# OptiVag:
# Tools and databases for annotating vaginal communities

## Citation
If using any part of this repo, please refer to _Hugerth et al, submitted_

## Contents
### Database
**db**

_16S:_ 

* optivag_db.aln.fasta.gz: aligned file with all 16S sequences used to simulate amplicons

*	optivag_db.fasta.gz: unaligned file with all 16S sequences used to simulate amplicons

*	optivag_seqinfo.csv: information on each of these sequences, including accession ID and taxonomy

_genome_info:_

* bacteria_list.tsv: list of bacteria, needed for creating a database locally

* updated_taxonomy.tsv: taxon names which changed since the inclusion in the database

**tools**

3 scripts, required for recreating the shotgun database from the files in _genome_info_

### Amplicon simulation
A single script, extracts amplicons and reads of a given length, given forward and reverse primer sequences 

### Shotgun tools
Two scripts:

* is_it_human.py: classifies reads in a fasta file as mapped or unmapped, given a reference file in UC format

*	make_roc_curve.py: classifies reads in one or more fastas files as correctly mapped, incorrectly mapped, correclty unmapped or incorrectly unmapped, given a reference file in UC or SAM format
