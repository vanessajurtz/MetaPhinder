# MetaPhinder

MetaPhinder classifies metagenomic contigs as of phage origin or not based on a comparison to a phage database. This tool is described in https://doi.org/10.1371/journal.pone.0163111.

# Installation

First, clone this repository:

```
git clone git@github.com:vanessajurtz/MetaPhinder.git
```

Change into the repository directory:

```
cd MetaPhinder
```

The following dependecies are required to run MetaPhinder:

* BLAST
* Python 2.7
* numpy

The easiest way of obtaining the necessary dependencies is through Anaconda.

To create a conda environment in the current directory (*preferred*):

```
conda create --prefix ./env -c bioconda python=2.7 blast numpy
conda activate ./env
```

Alternatively, you can create the environment in your `$HOME`:

```
conda create --name MetaPhinder -c bioconda python=2.7 blast numpy
conda activate MetaPhinder
```

# Usage

The following options are provided:

* `-i/--infile`: FASTA file of metagenomic contigs
* `-o/--outpath`: Path to which output files will be written (default: current directory)
* `-d/--database`: Path to MetaPhinder database (*e.g.* `database/ALL_140821_hr`)
* `-b/--blast`: Path to BLAST installation (*e.g.* `./env/bin`)

A help menu can be accessed with the `-h/--help` flag:

```
$ ./MetaPhinder.py --help
Usage: MetaPhinder.py [options]

Options:
  -h, --help            show this help message and exit
  -i INFILE, --infile=INFILE
                        contigs FASTA file format
  -o OUTPATH, --outpath=OUTPATH
                        path to output file(s)
  -d DATABASE, --database=DATABASE
                        MetaPhinder database
  -b BLAST, --blast=BLAST
                        path to BLAST installation
```

Running MetaPhinder on example data:

```
$ ./MetaPhinder.py -i phage_genomes/phages_1_hr.fsa -d database/ALL_140821_hr -b ./env/bin
parsing commandline options...
running BLAST...
calculating ANI...
preparing output...
DONE!
```

# Output

Two files are generated, `blast.out` and `output.txt`.

`output.txt` is a tab (`\t`) separated file containing the classifications by MetaPhinder:

```
$ head output.txt
#contigID       classification  ANI [%] merged coverage [%]     number of hits  size[bp]
NC_003714       phage   98.991  99.966  4       2948
AB091475        phage   91.533  99.969  9       3254
12021.1 phage   87.461  99.971  10      3486
NC_001628       phage   99.272  99.972  3       3588
NC_004175       phage   88.066  99.976  3       4100
NC_001890       phage   81.272  99.976  9       4215
FJ539134        phage   80.908  99.977  7       4273
NC_003438       phage   99.54   99.977  5       4421
NC_001741       phage   94.39   99.979  3       4877
```

The columns of the output file are:

* `#contigID`: Query sequence ID from FASTA input file
* `classification`: Query sequence classification; either *phage* or *negative*
* `ANI [%]`: Average nucleotide identity, relative to phage database
* `merged coverage [%]`: Coverage of the query sequence by all hits after overlapping hit regions have been merged
* `number of hits`: Number of blastn hits between one query sequence and one subject
* `size[bp]`: Query sequence length

*Note*: only contigs of length >500 are classified
