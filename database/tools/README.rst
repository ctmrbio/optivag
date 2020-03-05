======
How to use the tools in this repo to create your own reference database
======

1. Download the latest **assembly summaries** from  `RefSeq <ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt>`_
and `GenBank <ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt>`_ 

RefSeq: ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt

Genbank: ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
  
2. Use the script **extract_vaginal_taxa.py** and the file **bacteria_list.tsv** to extract the relevant genome headers from RefSeq:

  python extract_vaginal_taxa.py --assem assembly_summary_refseq.txt --touse bacteria_list.tsv --found refseq_genomes_found.tsv --notfound refseq_missing_genomes.tsv

3. Go through your missing genomes and assess whether there are taxa with recently updated taxonomy or any typos in your bacteria_list.tsv. If you have any updates to our current list of vaginal bacteria, please issue a pull request.

4. If necessary, repeat steps 2 and 3 until you're satisfied that you cannot find any more relevant genomes in RefSeq.

5. Use the script **filter_found_genomes.py** to create bash scripts to automatically download the relevant genomes for you (in our example, no more than 25 per species):

  python filter_found_genomes.py -i refseq_genomes_found.tsv -g get_refseq_genomes.sh -p get_refseq_proteins.sh -n 25
  
6. Run the produced bash scripts.

7. Repeat step 2 using the input files **assembly_summary_genbank.txt** and **refseq_missing_genomes.tsv** to expand your database with genomes still not available in RefSeq.

8. Repeat steps 3-6.
