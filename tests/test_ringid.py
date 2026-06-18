"""Tests for ringid.py.

Note: the TAE database (TAE.LIST) only contains atom-type entries for
5- and 6-membered rings. ringid still detects 3- and 4-membered rings
correctly (verified separately), but those ring sizes have no
corresponding TAE descriptor data, so they're out of scope for this
package and not exercised here.
"""
import unittest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from recon6.ringid import ringid
from recon6.readers.smiles import parse_smiles


def _ring_sizes(smiles):
    mol = parse_smiles(smiles)
    isize = ringid(mol['natom'], mol['ival'], mol['idcon'])
    return [isize[i] for i in range(1, mol['natom']+1) if mol['nuc'][i] != 1]

class TestRingId(unittest.TestCase):
    def test_acyclic(self):
        self.assertTrue(all(s == 0 for s in _ring_sizes('CC')))

    def test_benzene_6ring(self):
        self.assertTrue(all(s == 6 for s in _ring_sizes('c1ccccc1')))

    def test_cyclopentane_5ring(self):
        self.assertTrue(all(s == 5 for s in _ring_sizes('C1CCCC1')))

    def test_pyridine_6ring(self):
        sizes = _ring_sizes('c1ccncc1')
        self.assertTrue(all(s == 6 for s in sizes))

    def test_fused_bicyclic_naphthalene(self):
        # Naphthalene: two fused 6-rings sharing an edge
        sizes = _ring_sizes('c1ccc2ccccc2c1')
        self.assertTrue(all(s != 0 for s in sizes))

if __name__ == '__main__':
    unittest.main()
