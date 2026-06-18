"""Tests for the ORCA input reader (Cartesian '* xyz' block format
only)."""
import unittest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.dirname(__file__))

from recon6.readers.orca import read_orca_xyz
from recon6.bonds import load_bond_table

from test_config import DATA_DIR, HAVE_DATA


def _write(path, text):
    with open(path, 'w') as f:
        f.write(text)


_METHANE_SPACED = """! B3LYP def2-SVP Opt

* xyz 0 1
C        0.000000    0.000000    0.000000
H        0.629312    0.629312    0.629312
H        0.629312   -0.629312   -0.629312
H       -0.629312    0.629312   -0.629312
H       -0.629312   -0.629312    0.629312
*
"""

_METHYL_NO_SPACE_WITH_BLOCK = """! HF def2-SVP

%scf
  MaxIter 200
end

*xyz -1 2
C   0.000000    0.000000    0.000000
H   0.629312    0.629312    0.629312
H   0.629312   -0.629312   -0.629312
H  -0.629312    0.629312   -0.629312
*
"""

_FORMALDEHYDE = """! B3LYP def2-SVP

* xyz 0 1
C    0.000000    0.000000    0.000000
O    0.000000    0.000000    1.210000
*
"""

_BY_ATOMIC_NUMBER = """! HF def2-SVP

* xyz 0 1
6   0.000000    0.000000    0.000000
1   0.629312    0.629312    0.629312
1   0.629312   -0.629312   -0.629312
1  -0.629312    0.629312   -0.629312
1  -0.629312   -0.629312    0.629312
*
"""


class TestOrcaXyzParsing(unittest.TestCase):
    def setUp(self):
        if not HAVE_DATA:
            self.skipTest("Test data not available")
        self.bondl = load_bond_table(os.path.join(DATA_DIR, 'bond'))

    def test_methane_spaced_marker(self):
        path = '/tmp/test_orca_methane_spaced.inp'
        _write(path, _METHANE_SPACED)
        mol = read_orca_xyz(path, bondl=self.bondl)
        self.assertEqual(mol['natom'], 5)
        self.assertEqual(mol['mol_charge'], 0)
        self.assertEqual(mol['mol_multiplicity'], 1)
        self.assertEqual(mol['nuc'][1], 6)
        self.assertEqual(mol['ival'][1], 4)

    def test_no_space_marker_and_skipped_block(self):
        """'*xyz' (no space) after a '%scf ... end' block should
        parse identically to the spaced form, with the intervening
        block correctly skipped."""
        path = '/tmp/test_orca_noSpace.inp'
        _write(path, _METHYL_NO_SPACE_WITH_BLOCK)
        mol = read_orca_xyz(path, bondl=self.bondl)
        self.assertEqual(mol['natom'], 4)
        self.assertEqual(mol['mol_charge'], -1)
        self.assertEqual(mol['mol_multiplicity'], 2)
        self.assertEqual(mol['ival'][1], 3)

    def test_atomic_number_input(self):
        path = '/tmp/test_orca_atomicnum.inp'
        _write(path, _BY_ATOMIC_NUMBER)
        mol = read_orca_xyz(path, bondl=self.bondl)
        self.assertEqual(mol['natom'], 5)
        self.assertEqual(mol['nuc'][1], 6)
        for i in range(2, 6):
            self.assertEqual(mol['nuc'][i], 1)

    def test_carbonyl_bond_order_inference(self):
        path = '/tmp/test_orca_formaldehyde.inp'
        _write(path, _FORMALDEHYDE)
        mol = read_orca_xyz(path, bondl=self.bondl)
        from recon6.hydrogenate import _bond_order_sum, _best_standard_valence
        used_o = _bond_order_sum(mol, 2)
        target_o = _best_standard_valence(mol['nuc'][2], used_o)
        self.assertEqual(used_o, target_o, "carbonyl O should already satisfy its valence")

    def test_no_hydrogen_addition_capability_at_all(self):
        """ORCA input always specifies all atoms including hydrogens
        explicitly, so this reader does not offer hydrogen addition at
        all (no auto_add_h parameter)."""
        import inspect
        sig = inspect.signature(read_orca_xyz)
        self.assertNotIn('auto_add_h', sig.parameters)

    def test_heavy_atom_only_input_left_unsaturated(self):
        path = '/tmp/test_orca_formaldehyde.inp'
        _write(path, _FORMALDEHYDE)
        mol = read_orca_xyz(path, bondl=self.bondl)
        self.assertEqual(mol['natom'], 2)
        num_h = sum(1 for i in range(1, mol['natom'] + 1) if mol['nuc'][i] == 1)
        self.assertEqual(num_h, 0)


class TestOrcaMalformedInput(unittest.TestCase):
    def setUp(self):
        if not HAVE_DATA:
            self.skipTest("Test data not available")
        self.bondl = load_bond_table(os.path.join(DATA_DIR, 'bond'))

    def test_missing_xyz_block_raises(self):
        path = '/tmp/test_orca_bad1.inp'
        _write(path, "! HF def2-SVP\n\nno xyz block here\n")
        with self.assertRaises(ValueError):
            read_orca_xyz(path, bondl=self.bondl)

    def test_empty_atom_block_raises(self):
        path = '/tmp/test_orca_bad2.inp'
        _write(path, "! HF\n\n* xyz 0 1\n*\n")
        with self.assertRaises(ValueError):
            read_orca_xyz(path, bondl=self.bondl)

    def test_unrecognized_element_raises(self):
        path = '/tmp/test_orca_bad3.inp'
        _write(path, "! HF\n\n* xyz 0 1\nXx 0.0 0.0 0.0\n*\n")
        with self.assertRaises(ValueError):
            read_orca_xyz(path, bondl=self.bondl)

    def test_missing_terminator_does_not_crash(self):
        """A file with no closing '*' after the atom block should
        simply read to EOF rather than raising or hanging."""
        path = '/tmp/test_orca_noterm.inp'
        _write(path, "! HF\n\n* xyz 0 1\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
        mol = read_orca_xyz(path, bondl=self.bondl)
        self.assertEqual(mol['natom'], 2)


class TestOrcaEndToEnd(unittest.TestCase):
    def setUp(self):
        if not HAVE_DATA:
            self.skipTest("Test data not available")

    def test_run_recon_processes_orca_input(self):
        from recon6.recon import ReconConfig, run_recon
        path = '/tmp/test_orca_e2e.inp'
        _write(path, _METHANE_SPACED)
        config = ReconConfig(data_dir=DATA_DIR, input_files=[path], fmt='orca')
        results = run_recon(config)
        self.assertEqual(len(results), 1)
        self.assertIn('Energy', results[0])

    def test_auto_format_detection_inp_extension(self):
        from recon6.recon import ReconConfig, run_recon
        path = '/tmp/test_orca_autofmt.inp'
        _write(path, _METHANE_SPACED)
        config = ReconConfig(data_dir=DATA_DIR, input_files=[path], fmt='auto')
        results = run_recon(config)
        self.assertEqual(len(results), 1)

    def test_gaussian_and_orca_agree_on_same_molecule(self):
        """Methane specified via Gaussian .com and ORCA '* xyz' with
        identical coordinates should produce identical descriptors."""
        from recon6.recon import ReconConfig, run_recon
        gpath = '/tmp/test_cross_methane.com'
        opath = '/tmp/test_cross_methane.inp'
        _write(gpath,
               "%chk=t.chk\n# HF/3-21G\n\ntitle\n\n0 1\n"
               "C        0.000000    0.000000    0.000000\n"
               "H        0.629312    0.629312    0.629312\n"
               "H        0.629312   -0.629312   -0.629312\n"
               "H       -0.629312    0.629312   -0.629312\n"
               "H       -0.629312   -0.629312    0.629312\n\n")
        _write(opath, _METHANE_SPACED)

        config_g = ReconConfig(data_dir=DATA_DIR, input_files=[gpath], fmt='gaussian')
        config_o = ReconConfig(data_dir=DATA_DIR, input_files=[opath], fmt='orca')
        results_g = run_recon(config_g)
        results_o = run_recon(config_o)
        self.assertAlmostEqual(results_g[0]['Energy'], results_o[0]['Energy'], places=4)


if __name__ == '__main__':
    unittest.main()
