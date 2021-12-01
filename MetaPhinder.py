#!/usr/bin/env python3
"""
MetaPhinder
Purpose: Classify metagenomic contigs as phage or not
"""

from __future__ import division
from Bio import SeqIO
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
def get_contig_size(contig_file):
    """ Calculate contig lengths """

    contig_ids = []
    size = {}

    for rec in SeqIO.parse(contig_file, 'fasta'):
        contig_ids.append(rec.id)
        size[rec.id] = len(rec.seq)

    if len(contig_ids) == 0:
        sys.stderr.write("No contigs found! Problem with FASTA file format\n")
        sys.exit(2)

    return contig_ids, size


# ----------------------------------------------------------------------------
def calc_a_id(p_id, aln_l):
    """ Calculate AID """

    total_alignment = sum(aln_l)

    alignment_id = sum([id * align for id, align in zip(p_id, aln_l)])

    if total_alignment > 0 and alignment_id > 0:
        alignment_id = (alignment_id / total_alignment) / 100

    return alignment_id


# ----------------------------------------------------------------------------
def test_calc_a_id() -> None:
    """ Test AID calculation function """

    # Zero alignment or identity result in 0% AID
    assert calc_a_id([0], [0]) == 0.
    assert calc_a_id([100], [0]) == 0.
    assert calc_a_id([0], [100]) == 0.
    assert calc_a_id([0, 0, 0], [0, 0, 0]) == 0.

    # 100% identity results in 100% AID
    assert calc_a_id([100], [65]) == 1.
    assert calc_a_id([100, 100, 100], [57, 75, 87]) == 1.

    # Spot check equation
    assert calc_a_id([100, 75, 25], [100, 50, 50]) == 0.75
    assert calc_a_id([100, 25, 25], [1000, 1000, 1000]) == 0.5


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

    rel_mcov = float(mcov + 1) / gsize

    return rel_mcov


# ----------------------------------------------------------------------------
def test_calc_rel_mcov() -> None:
    """ Test function for calcualting relative merged coverage """

    # No coverage returns 0
    assert calc_rel_mcov([], 1000) == 0.
    
    # 100% coverage returns 1.0
    assert calc_rel_mcov([(1, 1000)], 1000) == 1.
    assert calc_rel_mcov([(1, 500), (501, 1000)], 1000) == 1.
    
    # Other spot checks
    assert calc_rel_mcov([(1, 500)], 1000) == 0.5

# ----------------------------------------------------------------------------
def main():
    """ Main function """

    args = get_args()
    contig_file = args.infile
    blast_db = args.database
    blast_path = args.blast
    out_path = args.outpath

    contig_ids, size = get_contig_size(contig_file)

    print("running BLAST...")

    blast = os.path.join(blast_path, 'blastn')
    blast_out = os.path.join(out_path, 'blast.out')

    os.system(f"{blast} -query {contig_file.name} -task blastn " +
              f"-evalue 0.05 -outfmt 7 -num_threads 4 -db {blast_db} " +
              f"-out {blast_out}")

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

    in_file = open(blast_out, "r")

    for line in in_file:
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

    print("preparing output...")

    out_file_name = os.path.join(out_path, 'output.txt')
    out_file = open(out_file_name, "w")

    out_file.write("#contig_ids\tclassification\tANI [%]\t"
                   "merged coverage [%]\tnumber of hits\tsize[bp]\n")

    threshold = 1.7

    for i in contig_ids:

        if int(size[i]) < 500:
            out_file.write(i + "\tnot processed\tnot processed\t" +
                           "not processed\tnot processed\t" + str(size[i]) +
                           "\n")
        elif i in res:
            ani = float(res[i].split("\t")[0])
            if ani > threshold:
                out_file.write(i + "\tphage\t" + res[i] + "\t" + str(size[i]) +
                               "\n")
            else:
                out_file.write(i + "\tnegative\t" + res[i] + "\t" +
                               str(size[i]) + "\n")
        else:
            out_file.write(i + "\tnegative\t0\t0\t0\t" + str(size[i]) + "\n")
    out_file.close()

    # for wrapper:
    sys.stderr.write("DONE!")


# ----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
