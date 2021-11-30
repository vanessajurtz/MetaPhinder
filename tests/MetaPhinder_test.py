""" tests """

from subprocess import getstatusoutput
import platform
import os
import re

PRG = './MetaPhinder.py'
RUN = 'python3 MetaPhinder.py' if platform.system() == 'Windows' else PRG
TEST1 = './tests/inputs/input1.fa'


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
