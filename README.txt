MetaPhinder

MetaPhinder classifies metagenomic contigs as of phage origin or not based on a
comparison to a phage database.
The script relies on BLAST which must be installed on you machine prior to running
MetaPhinder.
The input to MetaPhinder is a FASTA file of metagenomic contigs.

Here is an example of how to run MetaPhinder:

python MetaPinder.py -i infile.fsa -o output.txt -d database -b path_to_blast

OPTIONS:

-i inputfile in FASTA format
-o output file (text file)
-d location of MetaPhinders database
-b location of your local BLAST installation 
