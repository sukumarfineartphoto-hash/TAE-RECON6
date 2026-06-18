"""Tests for the Gaussian .com/.gjf reader (Cartesian coordinates
only, no z-matrix support)."""
import unittest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.dirname(__file__))

from recon6.readers.gaussian import read_gaussian_com, _parse_atom_symbol_or_number
from recon6.bonds import load_bond_table

from test_config import DATA_DIR, HAVE_DATA


def _write(path, text):
    with open(path, 'w') as f:
        f.write(text)


_METHANE = """%chk=test.chk
%mem=4GB
# B3LYP/6-31G(d) Opt Freq

Title line describing the molecule

0 1
C        0.000000    0.000000    0.000000
H        0.629312    0.629312    0.629312
H        0.629312   -0.629312   -0.629312
H       -0.629312    0.629312   -0.629312
H       -0.629312   -0.629312    0.629312

"""

_METHANE_BY_ATOMIC_NUMBER = """# HF/3-21G

Title line one
Title line two, still part of the title

-1 2
6        0.000000    0.000000    0.000000
1        0.629312    0.629312    0.629312
1        0.629312   -0.629312   -0.629312
1       -0.629312    0.629312   -0.629312

"""

_FORMALDEHYDE = """%chk=hcho.chk
# B3LYP/6-31G(d)

formaldehyde, heavy atoms only

0 1
C    0.000000    0.000000    0.000000
O    0.000000    0.000000    1.210000

"""


class TestAtomTokenParsing(unittest.TestCase):
    def test_element_symbol_carbon(self):
        self.assertEqual(_parse_atom_symbol_or_number('C'), 'C ')

    def test_element_symbol_chlorine_mixed_case(self):
        self.assertEqual(_parse_atom_symbol_or_number('Cl'), 'CL')

    def test_atomic_number_carbon(self):
        self.assertEqual(_parse_atom_symbol_or_number('6'), 'C ')

    def test_atomic_number_oxygen(self):
        self.assertEqual(_parse_atom_symbol_or_number('8'), 'O ')

    def test_unknown_token_returns_none(self):
        self.assertIsNone(_parse_atom_symbol_or_number('Xx'))

    def test_out_of_range_atomic_number_returns_none(self):
        self.assertIsNone(_parse_atom_symbol_or_number('99'))


class TestGaussianComParsing(unittest.TestCase):
    def setUp(self):
        if not HAVE_DATA:
            self.skipTest("Test data not available")
        self.bondl = load_bond_table(os.path.join(DATA_DIR, 'bond'))

    def test_methane_atom_count_and_symbols(self):
        path = '/tmp/test_gaussian_methane.com'
        _write(path, _METHANE)
        mol = read_gaussian_com(path, bondl=self.bondl)
        self.assertEqual(mol['natom'], 5)
        self.assertEqual(mol['atom'][1].strip(), 'C')
        self.assertEqual(mol['nuc'][1], 6)
        for i in range(2, 6):
            self.assertEqual(mol['nuc'][i], 1)

    def test_methane_charge_and_multiplicity(self):
        path = '/tmp/test_gaussian_methane.com'
        _write(path, _METHANE)
        mol = read_gaussian_com(path, bondl=self.bondl)
        self.assertEqual(mol['mol_charge'], 0)
        self.assertEqual(mol['mol_multiplicity'], 1)

    def test_methane_tetrahedral_connectivity(self):
        path = '/tmp/test_gaussian_methane.com'
        _write(path, _METHANE)
        mol = read_gaussian_com(path, bondl=self.bondl)
        self.assertEqual(mol['ival'][1], 4)
        for i in range(2, 6):
            self.assertEqual(mol['ival'][i], 1)
            self.assertEqual(mol['idcon'][i][1], 1)

    def test_atomic_number_input_and_negative_charge(self):
        path = '/tmp/test_gaussian_atomicnum.com'
        _write(path, _METHANE_BY_ATOMIC_NUMBER)
        mol = read_gaussian_com(path, bondl=self.bondl)
        self.assertEqual(mol['natom'], 4)
        self.assertEqual(mol['nuc'][1], 6)
        self.assertEqual(mol['mol_charge'], -1)
        self.assertEqual(mol['mol_multiplicity'], 2)

    def test_no_chk_line_still_parses(self):
        path = '/tmp/test_gaussian_nochk.com'
        _write(path, _METHANE_BY_ATOMIC_NUMBER)
        mol = read_gaussian_com(path, bondl=self.bondl)
        self.assertEqual(mol['natom'], 4)

    def test_carbonyl_bond_order_inference(self):
        """Same bond-order inference used for PDB applies here: a
        short C=O bond should satisfy oxygen's valence without an
        extra hydrogen being implied."""
        path = '/tmp/test_gaussian_formaldehyde.com'
        _write(path, _FORMALDEHYDE)
        mol = read_gaussian_com(path, bondl=self.bondl)
        from recon6.hydrogenate import _bond_order_sum, _best_standard_valence
        used_c = _bond_order_sum(mol, 1)
        used_o = _bond_order_sum(mol, 2)
        target_o = _best_standard_valence(mol['nuc'][2], used_o)
        self.assertEqual(used_o, target_o, "carbonyl O should already satisfy its valence")
        self.assertEqual(used_c, 2, "carbon should show bond order 2 from the C=O")

    def test_no_hydrogen_addition_capability_at_all(self):
        """Gaussian input always specifies all atoms including
        hydrogens explicitly, so this reader does not offer hydrogen
        addition at all (no auto_add_h parameter) - confirmed by
        checking the function signature directly."""
        import inspect
        sig = inspect.signature(read_gaussian_com)
        self.assertNotIn('auto_add_h', sig.parameters)

    def test_heavy_atom_only_input_left_unsaturated(self):
        """A heavy-atom-only Gaussian file (e.g. hand-edited or from
        an external tool) must be read as-is, with no hydrogens
        silently added."""
        path = '/tmp/test_gaussian_formaldehyde.com'
        _write(path, _FORMALDEHYDE)
        mol = read_gaussian_com(path, bondl=self.bondl)
        self.assertEqual(mol['natom'], 2)
        num_h = sum(1 for i in range(1, mol['natom'] + 1) if mol['nuc'][i] == 1)
        self.assertEqual(num_h, 0)


class TestGaussianComMalformedInput(unittest.TestCase):
    def setUp(self):
        if not HAVE_DATA:
            self.skipTest("Test data not available")
        self.bondl = load_bond_table(os.path.join(DATA_DIR, 'bond'))

    def test_missing_charge_multiplicity_line_raises(self):
        path = '/tmp/test_gaussian_bad1.com'
        _write(path, "# route\n\ntitle\n\n")
        with self.assertRaises(ValueError):
            read_gaussian_com(path, bondl=self.bondl)

    def test_unrecognized_element_raises(self):
        path = '/tmp/test_gaussian_bad2.com'
        _write(path, "# route\n\ntitle\n\n0 1\nXx 0.0 0.0 0.0\n\n")
        with self.assertRaises(ValueError):
            read_gaussian_com(path, bondl=self.bondl)

    def test_no_atoms_raises(self):
        path = '/tmp/test_gaussian_bad3.com'
        _write(path, "# route\n\ntitle\n\n0 1\n\n")
        with self.assertRaises(ValueError):
            read_gaussian_com(path, bondl=self.bondl)


class TestGaussianEndToEnd(unittest.TestCase):
    def setUp(self):
        if not HAVE_DATA:
            self.skipTest("Test data not available")

    def test_run_recon_processes_gaussian_com(self):
        from recon6.recon import ReconConfig, run_recon
        path = '/tmp/test_gaussian_e2e.com'
        _write(path, _METHANE)
        config = ReconConfig(data_dir=DATA_DIR, input_files=[path], fmt='gaussian')
        results = run_recon(config)
        self.assertEqual(len(results), 1)
        self.assertIn('Energy', results[0])

    def test_auto_format_detection_com_extension(self):
        from recon6.recon import ReconConfig, run_recon
        path = '/tmp/test_gaussian_autofmt.com'
        _write(path, _METHANE)
        config = ReconConfig(data_dir=DATA_DIR, input_files=[path], fmt='auto')
        results = run_recon(config)
        self.assertEqual(len(results), 1)

    def test_auto_format_detection_gjf_extension(self):
        from recon6.recon import ReconConfig, run_recon
        path = '/tmp/test_gaussian_autofmt.gjf'
        _write(path, _METHANE)
        config = ReconConfig(data_dir=DATA_DIR, input_files=[path], fmt='auto')
        results = run_recon(config)
        self.assertEqual(len(results), 1)


if __name__ == '__main__':
    unittest.main()
