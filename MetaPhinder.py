#!/usr/bin/env python3
"""
MetaPhinder
Purpose: Classify metagenomic contigs as phage or not
"""

from __future__ import division
from Bio import SeqIO
from Bio.Blast import NCBIXML
from subprocess import getstatusoutput
from typing import Dict, TextIO, Tuple, List, NamedTuple
import argparse
import numpy as np
import os
import sys

# pylint: disable=unspecified-encoding,consider-using-enumerate
# pylint: disable=too-many-branches


# ----------------------------------------------------------------------------
class ANI_results(NamedTuple):
    """ Results from classification by ANI% """
    id: str
    ani: float
    coverage: float
    hits: int
    classification: str
    size: int


# ----------------------------------------------------------------------------
class Args(NamedTuple):
    """ Command line arguments """
    infile: TextIO
    database: str
    blast: str
    outpath: str


# ----------------------------------------------------------------------------
def get_args() -> Args:
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

    args = parser.parse_args()

    rv, blast_path = getstatusoutput(f'which {args.blast}/blastn')

    if rv != 0:
        parser.error(f'BLAST not found on --blast path "{args.blast}".')

    args.blast = blast_path

    return Args(args.infile, args.database, args.blast, args.outpath)


# ----------------------------------------------------------------------------
def get_contig_size(contig_file: TextIO) -> Dict[str, int]:
    """ Calculate contig lengths """

    sizes: Dict[str, int] = {}

    for rec in SeqIO.parse(contig_file, 'fasta'):
        sizes[rec.id] = len(rec.seq)

    if len(sizes) == 0:
        sys.stderr.write("No contigs found! Problem with FASTA file format\n")
        sys.exit(2)

    return sizes


# ----------------------------------------------------------------------------
def calc_align_id(percent_ids: List[int], align_lengths: List[int]) -> float:
    """ Calculate AID """

    alignment_id: float = 0.
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
def calc_rel_mcov(positions: List[Tuple[int, int]], gsize: int) -> float:
    """ Genome-wide % Identity """

    rel_mcov = 0.
    coverages = [(0, 0)]

    for (start, end) in sorted(positions):
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
    assert calc_rel_mcov([(283, 359), (130, 359), (293, 347)], 1000) == 0.23


# ----------------------------------------------------------------------------
def calc_ani(blast_out: str, sizes: Dict[str, int]):
    """ Calculate ANI from BLAST results"""

    results: Dict[str, ANI_results] = {}  # results
    align_num = 0

    e_thresh = 0.05

    # Iterate through each query in BLAST output
    for query in NCBIXML.parse(open(blast_out, 'rt')):

        query_id = query.query.split(" ")[0]
        percent_ids = []
        align_lengths = []
        positions = []

        # Iterate through each alignment in query
        for align_num, alignment in enumerate(query.alignments, start=1):

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
        rel_mcov = calc_rel_mcov(positions, sizes[query_id])

        results[query_id] = summarize_results(query_id, align_id, rel_mcov,
                                              align_num, sizes[query_id])

    return results


# ----------------------------------------------------------------------------
def summarize_results(q_id: str, align_id: float, rel_mcov: float, hits: int,
                      size: int) -> ANI_results:
    """ Summarize query results """

    # Experimentally optimized threshold
    ANI_THRESH = 1.7

    ani = align_id * rel_mcov * 100

    classification = "phage" if ani > ANI_THRESH else "negative"

    return ANI_results(q_id, round(ani, 3), round(rel_mcov * 100, 3), hits,
                       classification, size)


# ----------------------------------------------------------------------------
def test_summarize_results() -> None:
    """ Test summarizing function """

    # Far over ANI threshold
    result1 = ANI_results('CZMN', 25., 50., 4, "phage", 500)
    assert summarize_results('CZMN', 0.5, 0.5, 4, 500) == result1

    # Just under ANI threshold
    result2 = ANI_results('PDFB', 1.6, 10., 1, "negative", 500)
    assert summarize_results('PDFB', 0.16, 0.1, 1, 500) == result2

    # Just at ANI threshold
    result3 = ANI_results('ARFS', 1.7, 10., 1, "phage", 500)
    assert summarize_results('ARFS', 0.17, 0.1, 1, 500) == result3

    # Check rounding
    result4 = ANI_results('ADFB', 15.665, 47., 3, "phage", 500)
    assert summarize_results('ADFB', 0.3333, 0.47, 3, 500) == result4


# ----------------------------------------------------------------------------
def get_out_str(result: ANI_results, contig_id: str, size: int) -> str:
    """ Determine output string for contig """

    out_str = ""

    if int(size) < 500:
        out_str = (f'{contig_id}\tnot processed\tnot processed\t'
                   f'not processed\tnot processed\t{str(size)}\n')
    elif result:
        ani = result.ani
        coverage = result.coverage
        hits = result.hits
        classification = result.classification

        out_str = (f'{contig_id}\t{classification}\t{ani}\t' +
                   f'{coverage}\t{hits}\t{str(size)}\n')
    else:
        out_str = (f'{contig_id}\tnegative\t0\t0\t0\t{str(size)}\n')

    return out_str


# ----------------------------------------------------------------------------
def test_get_out_str():
    """ Test get_out_str() """

    # No hits from BLAST, not present in results dict
    assert get_out_str(None, 'FOO', 1000) == "FOO\tnegative\t0\t0\t0\t1000\n"

    # Length < 500
    shorty = ANI_results('SMALL', 18.012, 14.012, 3, "phage", 499)
    exp = ("SMALL\tnot processed\tnot processed\t"
           "not processed\tnot processed\t499\n")
    assert get_out_str(shorty, 'SMALL', 499) == exp

    # Predicted phage
    phage = ANI_results('PHAGE', 99.024, 100.0, 4, "phage", 2948)
    exp = "PHAGE\tphage\t99.024\t100.0\t4\t2948\n"
    assert get_out_str(phage, 'PHAGE', 2948) == exp

    # Predicted negative
    bacteria = ANI_results('BACTERIA', 1.699, 22.524, 2, "negative", 3926)
    exp = "BACTERIA\tnegative\t1.699\t22.524\t2\t3926\n"
    assert get_out_str(bacteria, 'BACTERIA', 3926) == exp


# ----------------------------------------------------------------------------
def main():
    """ Main function """

    args = get_args()
    contig_file = args.infile
    blast_db = args.database
    blast = args.blast
    out_path = args.outpath

    if not os.path.isdir(out_path):
        os.makedirs(out_path)

    sizes = get_contig_size(contig_file)

    print("running BLAST...")

    blast_out = os.path.join(out_path, 'blast.out')

    rv, out = getstatusoutput(
        f"{blast} -query {contig_file.name} -task blastn " +
        f"-evalue 0.05 -outfmt 5 -num_threads 4 -db {blast_db} " +
        f"-out {blast_out}")

    if rv != 0:
        if os.path.isfile(blast_out):
            os.remove(blast_out)
        sys.exit(out)

    print("calculating ANI...")

    res = calc_ani(blast_out, sizes)

    out_file_name = os.path.join(out_path, 'output.txt')
    out_file = open(out_file_name, "w")

    out_file.write("#contig_ids\tclassification\tANI [%]\t"
                   "merged coverage [%]\tnumber of hits\tsize[bp]\n")

    for contig_id, size in sizes.items():

        result = res.get(contig_id)

        out_str = get_out_str(result, contig_id, size)

        out_file.write(out_str)
    out_file.close()

    # for wrapper:
    sys.stderr.write("DONE!\n")


# ----------------------------------------------------------------------------
if __name__ == "__main__":
    main()
