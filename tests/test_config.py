"""
Shared test configuration: paths to the TAE DATA directory and the
various sample molecule files used by the integration tests.

All of these can be overridden via environment variables so the test
suite is portable across machines, while defaulting to the paths used
during this package's own development sandbox. Tests that depend on
files which aren't present are skipped (via HAVE_* flags / each test's
own setUp/setUpClass), not failed - this suite is not expected to be
fully runnable without the original TAE DATA directory and several
non-redistributable sample molecule files.

To run these tests against your own data, set:
    RECON6_DATA_DIR        - directory of TAE .dat files (contains 'bond')
    RECON6_GAC_SDF         - GAC_withH.sdf test file (974 molecules)
    RECON6_GAC_FF          - matching Fortran recon.ff reference
    RECON6_TOXX_SDF        - toxx2.sdf test file (278 molecules)
    RECON6_TOXX_FF         - matching Fortran recon.ff reference
    RECON6_TOXX_MOL2       - toxx.mol2 (MOL2 conversion of the above)
    RECON6_BENZOX_TXT      - Benzoxazines.txt SMILES list (106 lines)
    RECON6_BENZOX_FF       - matching Fortran recon.ff reference
    RECON6_PDB_NNA         - NNA-7jyc.pdb sample (no CONECT, no H)
    RECON6_PDB_SV6         - SV6-7lb7.pdb sample (no CONECT, partial D)
"""
import os

DATA_DIR = os.environ.get('RECON6_DATA_DIR', '/home/claude/data/DATA')
BOND_FILE = os.path.join(DATA_DIR, 'bond')

GAC_SDF = os.environ.get('RECON6_GAC_SDF', '/home/claude/gac/GAC_withH.sdf')
GAC_FF = os.environ.get('RECON6_GAC_FF', '/home/claude/gac/recon.ff')

TOXX_SDF = os.environ.get('RECON6_TOXX_SDF', '/home/claude/toxx/toxx2.sdf')
TOXX_FF = os.environ.get('RECON6_TOXX_FF', '/home/claude/toxx/recon.ff')
TOXX_MOL2 = os.environ.get('RECON6_TOXX_MOL2', '/home/claude/pdbtest/toxx.mol2')

BENZOX_TXT = os.environ.get('RECON6_BENZOX_TXT', '/home/claude/benzox/Benzoxazines.txt')
BENZOX_FF = os.environ.get('RECON6_BENZOX_FF', '/home/claude/benzox/recon.ff')

PDB_NNA = os.environ.get('RECON6_PDB_NNA', '/home/claude/pdbtest/NNA-7jyc.pdb')
PDB_SV6 = os.environ.get('RECON6_PDB_SV6', '/home/claude/pdbtest/SV6-7lb7.pdb')

HAVE_DATA = os.path.exists(DATA_DIR)
