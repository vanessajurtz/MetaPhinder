""" tests """

from subprocess import getstatusoutput
import platform
import os
import random
import re
import shutil
import string

PRG = './MetaPhinder.py'
RUN = 'python3 MetaPhinder.py' if platform.system() == 'Windows' else PRG
INPUT1 = './tests/inputs/input1.fa'
INPUT2 = './tests/inputs/input3.fa'
DB = 'database/ALL_140821_hr'
BLAST = './env/bin'

# pylint: disable=unspecified-encoding


# --------------------------------------------------
def test_exists():
    """ Program exists """

    assert os.path.isfile(PRG)


# --------------------------------------------------
def test_usage():
    """ Usage """

    for flag in ['-h', '--help']:
        retval, out = getstatusoutput(f'{PRG} {flag}')
        assert retval == 0
        assert re.search('usage', out)


# --------------------------------------------------
def test_no_args() -> None:
    """ Dies on no args """

    rv, out = getstatusoutput(RUN)
    assert rv != 0
    assert re.match("usage", out, re.IGNORECASE)


# --------------------------------------------------
def test_missing_args() -> None:
    """ Dies on missing args """

    # Missing FASTA file argument
    rv, out = getstatusoutput(f'{RUN} -d {DB} -b {BLAST}')
    assert rv != 0
    assert re.match("usage", out, re.IGNORECASE)
    assert re.search('required: -i/--infile', out)

    # Missing database argument
    rv, out = getstatusoutput(f'{RUN} -i {INPUT1} -b {BLAST}')
    assert rv != 0
    assert re.match("usage", out, re.IGNORECASE)
    assert re.search('required: -d/--database', out)

    # Missing blast argument
    rv, out = getstatusoutput(f'{RUN} -i {INPUT1} -d {DB}')
    assert rv != 0
    assert re.match("usage", out, re.IGNORECASE)
    assert re.search('required: -b/--blast', out)


# --------------------------------------------------
def test_bad_infile() -> None:
    """ Dies on bad input file name """

    bad = random_string()

    rv, out = getstatusoutput(f'{RUN} -i {bad} -b {BLAST} -d {DB}')
    assert rv != 0
    assert re.search(f"No such file or directory: '{bad}'", out)


# --------------------------------------------------
def test_bad_blast() -> None:
    """ Dies when no BLAST found """

    bad = random_string()

    rv, out = getstatusoutput(f'{RUN} -i {INPUT1} -b {bad} -d {DB}')
    assert rv != 0
    assert re.search(f"{bad}/blastn: not found", out)


# --------------------------------------------------
def test_bad_db() -> None:
    """ Dies when bad database is given """

    bad = random_string()

    rv, out = getstatusoutput(f'{RUN} -i {INPUT1} -b {BLAST} -d {bad}')
    assert rv != 0
    assert re.search(f'"{bad}" is not a valid database', out)


# --------------------------------------------------
def test_runs_ok() -> None:
    """ Runs with good inputs """

    out_dir = "out_test"

    for in_file in [INPUT1, INPUT2]:
        try:
            if os.path.isdir(out_dir):
                shutil.rmtree(out_dir)

            os.makedirs(out_dir)

            rv, out = getstatusoutput(
                f'{RUN} -i {in_file} -b {BLAST} -d {DB} -o {out_dir}')

            assert rv == 0
            assert re.search("DONE!", out)
            assert os.path.isdir(out_dir)
            out_file1 = os.path.join(out_dir, 'blast.out')
            out_file2 = os.path.join(out_dir, 'output.txt')
            assert os.path.isfile(out_file1)
            assert os.path.isfile(out_file2)

            num_contigs = open(in_file).read().count(">")
            assert open(out_file2).read().count('\n') == num_contigs + 1

        finally:
            if os.path.isdir(out_dir):
                shutil.rmtree(out_dir)


# --------------------------------------------------
def random_string() -> str:
    """ generate a random string """

    k = random.randint(5, 10)
    return ''.join(random.choices(string.ascii_letters + string.digits, k=k))
