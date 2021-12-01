#!/usr/bin/env python3
"""
MetaPhinder
Purpose: Classify metagenomic contigs as phage or not
"""

from __future__ import division
import argparse
import sys
import os

# pylint: disable=unspecified-encoding,consider-using-enumerate
# pylint: disable=too-many-branches


# ----------------------------------------------------------------------------
def get_args():
    """ Parse command line arguments"""

    parser = argparse.ArgumentParser(
        description='Classify metagenomic contigs as phage or not',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i',
                        '--infile',
                        metavar='FILE',
                        type=argparse.FileType('rt'),
                        help="Input FASTA file",
                        required=True)
    parser.add_argument('-d',
                        '--database',
                        metavar='DB',
                        type=str,
                        help="MetaPhinder database",
                        required=True)
    parser.add_argument('-b',
                        '--blast',
                        metavar='BLAST',
                        type=str,
                        help="Path to BLAST installation",
                        required=True)
    parser.add_argument('-o',
                        '--outpath',
                        metavar='DIR',
                        type=str,
                        help="Path to output file(s)",
                        default='.')

    return parser.parse_args()


# ----------------------------------------------------------------------------
def get_contig_size(contigfile):
    """ Calculate contig lengths """

    contigID = []
    size = {}
    s = -1

    for line in contigfile:
        line = line.strip()
        if line[0] == ">":

            # save size of previous contig:
            if s != -1:
                size[contigID[-1]] = s
                s = 0
            # save ID of new contig:
            line = line.split(" ")
            contigID.append(line[0].strip(">"))
        else:
            # count bases
            if s == -1:
                s = 0
            s = s + len(line)
    contigfile.close()

    # save size of last contig:
    if s == -1:
        sys.stderr.write("No contigs found! Problem with FASTA file format\n")
        sys.exit(2)
    else:
        size[contigID[-1]] = s

    return contigID, size


# ----------------------------------------------------------------------------
def calc_a_id(p_id, aln_l):
    """ Calculate ANI """

    tot_aln_l = 0
    a_id = 0
    for i in range(0, len(p_id)):
        a_id = a_id + (p_id[i] * aln_l[i])
        tot_aln_l = tot_aln_l + aln_l[i]

    if (tot_aln_l > 0 and a_id > 0):
        a_id = (a_id / tot_aln_l) / 100
    return a_id


# ----------------------------------------------------------------------------
def calc_rel_mcov(positions, gsize):
    """ Genome-wide % Identity """

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


# ----------------------------------------------------------------------------
def main():
    """ Main function """

    args = get_args()
    contigfile = args.infile
    blastDB = args.database
    blastPath = args.blast
    outPath = args.outpath

    ##########################################################################
    #   GET CONTIG LENGTH
    ##########################################################################

    contigID, size = get_contig_size(contigfile)

    ##########################################################################
    #   RUN BLAST
    ##########################################################################

    print("running BLAST...")

    blast = os.path.join(blastPath, 'blastn')
    blast_out = os.path.join(outPath, 'blast.out')

    os.system(blast + " -query " + contigfile.name +
              " -task blastn -evalue 0.05 -outfmt 7  -num_threads 4 -db " +
              blastDB + " -out " + blast_out)

    ##########################################################################
    #   PARSE BLAST OUTPUT
    ##########################################################################

    print("calculating ANI...")

    res = {}  # results
    p_id = []  # percent identity
    aln_l = []  # alignment length
    positions = []  # start and stop positions in alignment

    old_id = ''
    s_id = ''  # subject ID
    n_s_id = 0  # count hits in DB
    count = 0

    evalue = 0.05

    infile = open(blast_out, "r")

    for line in infile:
        line = line.strip()

        if line[0] != "#":
            line = line.split("\t")

            if (old_id != str(line[0])) and (old_id != ""):

                # calc average %ID, relative merged coverage
                # and genomewide %ID:
                a_id = calc_a_id(p_id, aln_l)
                rel_mcov = calc_rel_mcov(positions, size[old_id])

                g_id = a_id * rel_mcov

                # save result:
                res[old_id] = str(round(g_id * 100, 3)) + "\t" + str(
                    round(rel_mcov * 100, 3)) + "\t" + str(n_s_id)

                # reset variables:
                p_id = []
                aln_l = []
                positions = []
                old_id = str(line[0])
                count = count + 1
                s_id = ''
                n_s_id = 0

            # check for evalue:
            if float(line[10]) <= evalue:
                # save output:
                if s_id != line[1]:
                    s_id = line[1]
                    n_s_id = n_s_id + 1

                p_id.append(float(line[2]))
                aln_l.append(int(line[3]))
                if int(line[6]) < int(line[7]):
                    positions.append((int(line[6]), int(line[7])))
                else:
                    positions.append((int(line[7]), int(line[6])))
                old_id = str(line[0])

    if (old_id != str(line[0])) and (old_id != ""):

        # calc average %ID, relative merged coverage and genomewide %ID:
        a_id = calc_a_id(p_id, aln_l)
        rel_mcov = calc_rel_mcov(positions, size[old_id])

        g_id = a_id * rel_mcov

        # save result:
        res[old_id] = str(round(g_id * 100, 3)) + "\t" + str(
            round(rel_mcov * 100, 3)) + "\t" + str(n_s_id)

    ##########################################################################
    #   PRINT RESULTS
    ##########################################################################

    print("preparing output...")

    out_file_name = os.path.join(outPath, 'output.txt')
    outfile = open(out_file_name, "w")

    outfile.write(
        "#contigID\tclassification\tANI [%]\t"
        "merged coverage [%]\tnumber of hits\tsize[bp]\n"
    )

    threshold = 1.7

    for i in contigID:

        if int(size[i]) < 500:
            outfile.write(
                i +
                "\tnot processed\tnot processed\t" +
                "not processed\tnot processed\t"
                + str(size[i]) + "\n")
        elif i in res:
            ani = float(res[i].split("\t")[0])
            if ani > threshold:
                outfile.write(i + "\tphage\t" + res[i] + "\t" + str(size[i]) +
                              "\n")
            else:
                outfile.write(i + "\tnegative\t" + res[i] + "\t" +
                              str(size[i]) + "\n")
        else:
            outfile.write(i + "\tnegative\t0\t0\t0\t" + str(size[i]) + "\n")
    outfile.close()

    # for wrapper:
    sys.stderr.write("DONE!")


# ----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
