#!/usr/bin/env python

from __future__ import division

from optparse import OptionParser
import subprocess
import sys
import os
import numpy as np
import re
import math


##########################################################################
#   SUBROUTINES
##########################################################################

def get_contig_size(contigfile):

    # open inputfile (contigs in FASTA format):
    infile = open(contigfile, "r")

    contigID = []
    size = {}
    s = -1

    for l in infile:
        l = l.strip()
        if l[0] == ">":

            # save size of previous contig:
            if s != -1:
                size[contigID[-1]] = s
                s = 0
            # save ID of new contig:
            l = l.split(" ")
            contigID.append(l[0].strip(">"))
        else:
            # count bases
            if s == -1:
                s = 0
            s = s + len(l)
    infile.close()

    # save size of last contig:
    if s == -1:
        sys.stderr.write("No contigs found! Problem with FASTA file format\n")
        sys.exit(2)
    else:
        size[contigID[-1]] = s

    return contigID, size


# average %ID---------------------------------------
def calc_a_id(p_id, aln_l):

    tot_aln_l = 0
    a_id = 0
    for i in range(0, len(p_id)):
        a_id = a_id + (p_id[i] * aln_l[i])
        tot_aln_l = tot_aln_l + aln_l[i]

    if (tot_aln_l > 0 and a_id > 0):
        a_id = (a_id / tot_aln_l) / 100
    return a_id

# genomewide %ID----------------------------------


def calc_rel_mcov(positions, gsize):
    # calculate merged coverage:
    mcov = 0
    rel_mcov = 0

    if len(positions) > 1:
        spos = sorted(positions, key=lambda firstp: firstp[0])

        start = spos[0][0]
        end = spos[0][1]
        for i in range(0, (len(spos) - 1)):
            if spos[i + 1][0] > end:
                mcov += end - start
                start = spos[i + 1][0]
                end = spos[i + 1][1]
            else:
                if spos[i + 1][1] > end:
                    end = spos[i + 1][1]

        mcov += end - start
    # only one hit:
    elif len(positions) == 1:
        mcov = positions[0][1] - positions[0][0]

    rel_mcov = float(mcov) / gsize

    return rel_mcov


################################################################################
#       PARSE COMMANDLINE OPTIONS
################################################################################



parser = OptionParser()
parser.add_option('-i','--infile', dest="infile", help="contigs FASTA file format")
parser.add_option('-o','--outpath', dest="outpath", help="path to output file(s)")
parser.add_option('-d','--database', dest="database", help="MetaPhinder database")
parser.add_option('-b','--blast', dest="blast", help="path to BLAST installation")
(options, args) = parser.parse_args()

print "parsing commandline options..."

# open input file:
if options.infile != None:
    contigfile = options.infile
else:
    sys.stderr.write("Please specify inputfile!\n")
    sys.exit(2)

# save outfile:
if options.outpath != None:
    if options.outpath[-1] != "/":
        outPath = options.outpath + "/"
    else:
        outPath = options.outpath

else:
    outPath = ""

# save databasepath:
if options.database != None:
    blastDB = options.database
else:
    sys.stderr.write("Please specify path to database!\n")
    sys.exit(2)

# save path to nnlinplayer:
if options.blast != None:
    if options.blast[-1] != "/":
        blastPath = options.blast + "/"
    else:
        blastPath = options.blast
else:
    blastPath=""


################################################################################
#   GET CONTIG LENGTH
################################################################################


contigID,size = get_contig_size(contigfile)

################################################################################
#   RUN BLAST
################################################################################

print "running BLAST..."

os.system(blastPath + "blastn -query " + contigfile + " -task blastn -evalue 0.05 -outfmt 7  -num_threads 4 -db "
    + blastDB + " -out " + outPath + "blast.out")



################################################################################
#   PARSE BLAST OUTPUT
################################################################################

print "calculating ANI..."

res={} #results
p_id = []  # percent identity
aln_l = []  # alignment length
positions = []  # start and stop positions in alignment

old_id = ''
s_id='' # subject ID
n_s_id=0 # count hits in DB
count = 0

evalue=0.05

infile=open(outPath + "blast.out","r")

for l in infile:
    l = l.strip()

    if l[0] != "#":
        l = l.split("\t")

        if (old_id != str(l[0])) and (old_id != ""):

            # calc average %ID, relative merged coverage and genomewide %ID:
            a_id = calc_a_id(p_id, aln_l)
            rel_mcov = calc_rel_mcov(positions, size[old_id])

            g_id = a_id * rel_mcov

            #save result:
            res[old_id] = str(round(g_id*100,3)) + "\t" + str(round(rel_mcov*100,3)) + "\t"+ str(n_s_id)

            # reset variables:
            p_id = []
            aln_l = []
            positions = []
            old_id = str(l[0])
            count = count + 1
            s_id=''
            n_s_id=0

        #check for evalue:
        if float(l[10]) <= evalue:
            # save output:
            if s_id != l[1]:
                s_id = l[1]
                n_s_id = n_s_id +1

            p_id.append(float(l[2]))
            aln_l.append(int(l[3]))
            if int(l[6]) < int(l[7]):
                positions.append((int(l[6]), int(l[7])))
            else:
                positions.append((int(l[7]), int(l[6])))
            old_id = str(l[0])

if (old_id != str(l[0])) and (old_id != ""):

    # calc average %ID, relative merged coverage and genomewide %ID:
    a_id = calc_a_id(p_id, aln_l)
    rel_mcov = calc_rel_mcov(positions, size[old_id])

    g_id = a_id * rel_mcov

    #save result:
    res[old_id] = str(round(g_id*100,3)) + "\t" + str(round(rel_mcov*100,3)) + "\t"+ str(n_s_id)



################################################################################
#   PRINT RESULTS
################################################################################

print "preparing output..."

outfile=open(outPath + "output.txt","w")

outfile.write("#contigID\tclassification\tANI [%]\tmerged coverage [%]\tnumber of hits\tsize[bp]\n")

threshold=1.7

for i in contigID:

    if int(size[i])<500:
        outfile.write(i + "\tnot processed\tnot processed\tnot processed\tnot processed\t" + str(size[i]) + "\n")
    elif i in res:
        ani=float(res[i].split("\t")[0])
        if ani>threshold:
            outfile.write(i + "\tphage\t" + res[i] + "\t" + str(size[i]) +  "\n")
        else:
            outfile.write(i + "\tnegative\t" + res[i] + "\t" + str(size[i]) +  "\n")
    else:
        outfile.write(i + "\tnegative\t0\t0\t0\t" + str(size[i]) + "\n")
outfile.close()

# for wrapper:
sys.stderr.write("DONE!")
