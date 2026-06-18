"""Tests for molecular file format readers."""
import unittest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.dirname(__file__))

from recon6.readers.smiles import parse_smiles
from recon6.readers.sdf import read_sdf
from recon6.element import ele
from recon6.periodic import atomic_number
from test_config import GAC_SDF as _GAC_SDF


class TestEle(unittest.TestCase):
    def test_carbon(self):     self.assertEqual(ele('C'), 'C ')
    def test_chlorine(self):   self.assertEqual(ele('Cl'), 'CL')
    def test_chlorine_cap(self): self.assertEqual(ele('CL'), 'CL')
    def test_bromine(self):    self.assertEqual(ele('Br'), 'BR')
    def test_silicon(self):    self.assertEqual(ele('Si'), 'SI')
    def test_nitrogen(self):   self.assertEqual(ele('N'), 'N ')
    def test_unknown(self):    self.assertEqual(ele('Zn'), 'Zn')

class TestAtomicNumber(unittest.TestCase):
    def test_H(self): self.assertEqual(atomic_number('H '), 1)
    def test_C(self): self.assertEqual(atomic_number('C '), 6)
    def test_Cl(self): self.assertEqual(atomic_number('CL'), 17)
    def test_unknown(self): self.assertIsNone(atomic_number('XX'))

class TestSmiles(unittest.TestCase):
    def test_methane(self):
        mol = parse_smiles('C')
        self.assertEqual(mol['natom'], 5)
        self.assertEqual(mol['nuc'][1], 6)

    def test_ethane(self):
        mol = parse_smiles('CC')
        self.assertEqual(mol['natom'], 8)

    def test_water(self):
        mol = parse_smiles('O')
        self.assertEqual(mol['natom'], 3)

    def test_ethanol(self):
        mol = parse_smiles('CCO')
        self.assertEqual(mol['natom'], 9)

    def test_benzene(self):
        mol = parse_smiles('c1ccccc1')
        self.assertEqual(mol['natom'], 12)

    def test_cyclopropane_ring_connectivity(self):
        mol = parse_smiles('C1CC1')
        heavy = [i for i in range(1, mol['natom']+1) if mol['nuc'][i] == 6]
        self.assertEqual(len(heavy), 3)
        for c in heavy:
            heavy_nbrs = [mol['idcon'][c][j] for j in range(1, mol['ival'][c]+1)
                          if mol['nuc'][mol['idcon'][c][j]] == 6]
            self.assertEqual(len(heavy_nbrs), 2)

    def test_rejects_leading_paren(self):
        with self.assertRaises(ValueError):
            parse_smiles('(C)CC')

    def test_rejects_stereo_brackets(self):
        with self.assertRaises(ValueError):
            parse_smiles('C[C@H](N)O')

    def test_rejects_charge_brackets(self):
        with self.assertRaises(ValueError):
            parse_smiles('O=[N+]([O-])c1ccccc1')

    def test_rejects_bond_direction(self):
        with self.assertRaises(ValueError):
            parse_smiles('C(/C=C/C)C')

class TestSdf(unittest.TestCase):
    GAC_SDF = _GAC_SDF

    def test_gac_natom(self):
        if not os.path.exists(self.GAC_SDF):
            self.skipTest("GAC SDF not available")
        with open(self.GAC_SDF) as fh:
            mol = read_sdf(fh)
        # GAC_withH.sdf is a 974-molecule batch file; the first molecule
        # (C8H12N8S2) has 18 heavy+H atoms.
        self.assertEqual(mol['natom'], 18)
        self.assertEqual(mol['nuc'][1], 7)

if __name__ == '__main__':
    unittest.main()
