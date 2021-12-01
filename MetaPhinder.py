#!/usr/bin/env python3
"""
MetaPhinder
Purpose: Classify metagenomic contigs as phage or not
"""

from __future__ import division
from Bio import SeqIO
from Bio.Blast import NCBIXML
import argparse
import numpy as np
import os
import sys

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
def calc_align_id(percent_ids, align_lengths):
    """ Calculate AID """

    total_alignment = sum(align_lengths)

    alignment_id = sum(
        [id * align for id, align in zip(percent_ids, align_lengths)])

    if total_alignment > 0 and alignment_id > 0:
        alignment_id = (alignment_id / total_alignment) / 100

    return alignment_id


# ----------------------------------------------------------------------------
def test_calc_align_id() -> None:
    """ Test AID calculation function """

    # Zero alignment or identity result in 0% AID
    assert calc_align_id([0], [0]) == 0.
    assert calc_align_id([100], [0]) == 0.
    assert calc_align_id([0], [100]) == 0.
    assert calc_align_id([0, 0, 0], [0, 0, 0]) == 0.

    # 100% identity results in 100% AID
    assert calc_align_id([100], [65]) == 1.
    assert calc_align_id([100, 100, 100], [57, 75, 87]) == 1.

    # Spot check equation
    assert calc_align_id([100, 75, 25], [100, 50, 50]) == 0.75
    assert calc_align_id([100, 25, 25], [1000, 1000, 1000]) == 0.5


# ----------------------------------------------------------------------------
def calc_rel_mcov(positions, gsize):
    """ Genome-wide % Identity """

    rel_mcov = 0
    coverages = [(0, 0)]

    for (start, end) in positions:
        overlapping = [
            cov_start <= start <= cov_end for (cov_start, cov_end) in coverages
        ]
        if any(overlapping):
            i = int(np.where(overlapping)[0])
            (cov_start, cov_end) = coverages[i]
            coverages[i] = (cov_start, max(end, cov_end))
        else:
            coverages.append((start, end))

    rel_mcov = sum([end - (start - 1) for start, end in coverages[1:]]) / gsize

    return rel_mcov


# ----------------------------------------------------------------------------
def test_calc_rel_mcov() -> None:
    """ Test function for calculating relative merged coverage """

    # No coverage returns 0
    assert calc_rel_mcov([], 1000) == 0.

    # 100% coverage returns 1.0
    assert calc_rel_mcov([(1, 1000)], 1000) == 1.
    assert calc_rel_mcov([(1, 500), (501, 1000)], 1000) == 1.

    # Other spot checks
    assert calc_rel_mcov([(1, 500)], 1000) == 0.5
    assert calc_rel_mcov([(501, 1000)], 1000) == 0.5
    assert calc_rel_mcov([(1, 500), (250, 750)], 1000) == 0.75
    assert calc_rel_mcov([(1, 500), (501, 750)], 1000) == 0.75
    assert calc_rel_mcov([(501, 750), (1, 500)], 1000) == 0.75
    assert calc_rel_mcov([(1, 500), (250, 525), (501, 750)], 1000) == 0.75


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
              f"-evalue 0.05 -outfmt 5 -num_threads 4 -db {blast_db} " +
              f"-out {blast_out}")

    print("calculating ANI...")

    res = {}  # results
    align_num = 0

    e_thresh = 0.05
    ani_thresh = 1.7

    # Iterate through each query in BLAST output
    for query in NCBIXML.parse(open(blast_out, 'rt')):

        percent_ids = []
        align_lengths = []
        positions = []

        # Iterate through each alignment in query
        for align_num, alignment in enumerate(query.alignments, start=1):

            # Only first alignment retains original record ID
            if align_num == 1:
                query_id = alignment.title.split(" ")[0]

            # Check that alignment is below e-value threshold, and gather stats
            for hsp in alignment.hsps:
                if hsp.expect <= e_thresh:
                    align_length = hsp.align_length
                    align_lengths.append(align_length)

                    identity = (hsp.positives / align_length) * 100
                    percent_ids.append(identity)

                    start = hsp.query_start
                    end = hsp.query_end
                    positions.append((min(start, end), max(start, end)))

        align_id = calc_align_id(percent_ids, align_lengths)
        rel_mcov = calc_rel_mcov(positions, size[query_id])

        g_id = align_id * rel_mcov

        # save result:
        res[query_id] = str(round(g_id * 100, 3)) + "\t" + str(
            round(rel_mcov * 100, 3)) + "\t" + str(align_num)

    out_file_name = os.path.join(out_path, 'output.txt')
    out_file = open(out_file_name, "w")

    out_file.write("#contig_ids\tclassification\tANI [%]\t"
                   "merged coverage [%]\tnumber of hits\tsize[bp]\n")

    for i in contig_ids:

        if int(size[i]) < 500:
            out_file.write(i + "\tnot processed\tnot processed\t" +
                           "not processed\tnot processed\t" + str(size[i]) +
                           "\n")
        elif i in res:
            ani = float(res[i].split("\t")[0])
            if ani > ani_thresh:
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
