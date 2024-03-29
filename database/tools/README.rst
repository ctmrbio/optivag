======
How to use the tools in this repo to create your own reference database
======


Part 1: Creating a genome database
------------
1. Download the latest **assembly summaries** from  `RefSeq <ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt>`_
and `GenBank <ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt>`_ 

* RefSeq: ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
* Genbank: ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
  
2. Use the script **extract_vaginal_taxa.py** and the file
   **bacteria_list.tsv** (found under `database/db/genome_info`) to extract the
   relevant genome headers from RefSeq. Adjust the version number and output folder as appropriate::

    python optivag/database/tools/extract_vaginal_taxa.py \
     --assem assembly_summary_refseq.txt \
     --to-use optivag/database/db/genome_info/vXX/bacteria_list.tsv \
     --found outdir/refseq_genomes_found.tsv \
     --not-found outdir/refseq_missing_genomes.tsv

3. Go through your missing genomes and assess whether there are taxa with
   recently updated taxonomy or any typos in your bacteria_list.tsv. If you
   have any updates to our current list of vaginal bacteria, please issue a
   pull request.

4. If necessary, repeat steps 2 and 3 until you're satisfied that you cannot
   find any more relevant genomes in RefSeq.

5. Use the script **filter_found_genomes.py** to create bash scripts to
   automatically download the relevant genomes for you (in our example, no more
   than 25 per species)::

    python filter_found_genomes.py \
     -i refseq_genomes_found.tsv \
     -g get_refseq_genomes.sh \
     -p get_refseq_proteins.sh \
     -t refseq_taxa.tsv \
     -n 25
  
6. Run the produced bash scripts. These scripts may crash at unexpected symbols
   such as unpaired parentheses in the output genome name. As of right now
   these need to be manually removed.

7. Repeat step 2 using the input files **assembly_summary_genbank.txt** and
   **refseq_missing_genomes.tsv** to expand your database with genomes still
   not available in RefSeq.

8. Repeat steps 3-6.


Part 2: Creating a searchable Kraken2 database
---------------
While you can do anything you want with your brand new custom-made database, we
tend to run `Kraken2 <https://ccb.jhu.edu/software/kraken2/>`_ followed by
`_Bracken <https://ccb.jhu.edu/software/bracken/>`_. So we have some extra
tools to enable this.

9. Unzip the genomes you've downloaded. This will be needed for Kraken anyway.

10. Concatenate your taxon id lists::

    cat refseq_taxa.tsv genbank_taxa.tsv > all_taxa.tsv

11. Use the script **map_seq_tax.py** to associate each sequence in your
    download to its NCBI tax_id::

    python map_seq_tax.py -t all_taxa.tsv -g path/to/mygenomes
    
    You might at this point discover other special characters, such as "/", and yes, this has to be fixed manually.
  
12. You're now ready to follow `Kraken2's <https://ccb.jhu.edu/software/kraken2/>`_ 
    and `_Bracken's <https://ccb.jhu.edu/software/bracken/>`_ manual pages! You 
    will need to use Kraken2 to add each downloaded genome to your new
    database, and provide the seqid2taxid.map as input for Bracken.
