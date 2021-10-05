# OptiVag:
# Tools and databases for annotating vaginal communities

## Citation
If using any part of this repo, please refer to:

Luisa W. Hugerth, Marcela Pereira, Yinghua Zha, Maike Seifert, Vilde Kaldhusdal, Fredrik Boulund, Maria C. Krog, , Zahra Bashir, Marica Hamsten, Emma Fransson, Henriette Svarre Nielsen, Ina Schuppe-Koistinen, and Lars Engstrand (2018) [_Assessment of In Vitro and In Silico Protocols for Sequence-Based Characterization of the Human Vaginal Microbiome_](https://journals.asm.org/doi/full/10.1128/mSphere.00448-20) mSphere, 5(6): e00448-20


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

For instructions on how to create your local database, look [here](https://github.com/ctmrbio/optivag/tree/master/database/tools)

### Amplicon simulation
A single script, extracts amplicons and reads of a given length, given forward and reverse primer sequences 

### Shotgun tools
Two scripts:

* is_it_human.py: classifies reads in a fasta file as mapped or unmapped, given a reference file in UC format

*	make_roc_curve.py: classifies reads in one or more fastas files as correctly mapped, incorrectly mapped, correclty unmapped or incorrectly unmapped, given a reference file in UC or SAM format
