"""Tests for the hydrogenate.py H-addition utility."""
import unittest
import sys, os, math, itertools
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.dirname(__file__))

from recon6.hydrogenate import (
    add_missing_hydrogens, needs_hydrogens, _new_h_directions, _unit,
)
from recon6.sparse import SparseMatrix
from test_config import GAC_SDF as _GAC_SDF


def _make_mol(natom, atoms, nuc, coords, bonds, charge=None):
    MAX = natom + 1
    idcon = [[0] * 5 for _ in range(MAX + 1)]
    ival = [0] * (MAX + 1)
    nbo = SparseMatrix()
    for (i, j, bo) in bonds:
        ival[i] += 1; idcon[i][ival[i]] = j
        ival[j] += 1; idcon[j][ival[j]] = i
        nbo[i][j] = bo; nbo[j][i] = bo
    icon = [[0] * 5 for _ in range(MAX + 1)]
    isum = [0] * (MAX + 1)
    for i in range(1, natom + 1):
        s = 0
        for k in range(1, ival[i] + 1):
            nb = nuc[idcon[i][k]]
            icon[i][k] = nb
            s += nb
        isum[i] = s
    return dict(natom=natom, numbond=len(bonds), atom=atoms, nuc=nuc,
                coords=coords, idcon=idcon, ival=ival, nbo=nbo,
                icon=icon, isum=isum, charge=charge or {})


def _angle(a, b, c):
    v1 = [a[k] - b[k] for k in range(3)]
    v2 = [c[k] - b[k] for k in range(3)]
    dot = sum(v1[k] * v2[k] for k in range(3))
    n1 = math.sqrt(sum(x * x for x in v1))
    n2 = math.sqrt(sum(x * x for x in v2))
    return math.degrees(math.acos(dot / (n1 * n2)))


def _bond_len(a, b):
    return math.sqrt(sum((a[k] - b[k]) ** 2 for k in range(3)))


class TestHydrogenateGeometry(unittest.TestCase):
    def test_isolated_carbon_gets_4_tetrahedral_h(self):
        # Atom 1 = C at origin, no bonds at all
        mol = _make_mol(1, [None, 'C '], [0, 6], [[0, 0, 0], [0, 0, 0]], [])
        result = add_missing_hydrogens(mol)
        self.assertEqual(result['hydrogens_added'], 4)
        self.assertEqual(result['natom'], 5)
        c = result['coords'][1]
        hs = [result['coords'][i] for i in range(2, 6)]
        for h in hs:
            self.assertAlmostEqual(_bond_len(c, h), 1.09, places=3)
        for a, b in itertools.combinations(hs, 2):
            self.assertAlmostEqual(_angle(a, c, b), 109.47, places=1)

    def test_carbon_with_one_neighbor_gets_3_h_tetrahedral(self):
        # C1-C2, both need 3 H's each for a proper tetrahedral ethane
        mol = _make_mol(2, [None, 'C ', 'C '], [0, 6, 6],
                         [[0, 0, 0], [0, 0, 0], [1.54, 0, 0]], [(1, 2, 1)])
        result = add_missing_hydrogens(mol)
        self.assertEqual(result['hydrogens_added'], 6)
        c1, c2 = result['coords'][1], result['coords'][2]
        h_on_c1 = [result['coords'][i] for i in range(1, result['natom'] + 1)
                   if result['nuc'][i] == 1 and result['idcon'][i][1] == 1]
        self.assertEqual(len(h_on_c1), 3)
        for h in h_on_c1:
            self.assertAlmostEqual(_angle(h, c1, c2), 109.5, delta=0.5)
        for a, b in itertools.combinations(h_on_c1, 2):
            self.assertAlmostEqual(_angle(a, c1, b), 109.47, delta=0.5)

    def test_sp3_two_neighbors_gets_2_h_tetrahedral(self):
        # Proper tetrahedral C1-C2-C3 chain (C2 needs exactly 2 H's)
        def norm(v):
            n = math.sqrt(sum(x * x for x in v))
            return [x / n for x in v]
        d_c1 = norm([1, 1, 1])
        d_c3 = norm([1, -1, -1])
        c2 = [0, 0, 0]
        c1 = [c2[k] + d_c1[k] * 1.54 for k in range(3)]
        c3 = [c2[k] + d_c3[k] * 1.54 for k in range(3)]
        mol = _make_mol(3, [None, 'C ', 'C ', 'C '], [0, 6, 6, 6],
                         [[0, 0, 0], c1, c2, c3], [(1, 2, 1), (2, 3, 1)])
        result = add_missing_hydrogens(mol)
        c1r, c2r, c3r = result['coords'][1], result['coords'][2], result['coords'][3]
        h_on_c2 = [result['coords'][i] for i in range(1, result['natom'] + 1)
                   if result['nuc'][i] == 1 and result['idcon'][i][1] == 2]
        self.assertEqual(len(h_on_c2), 2)
        self.assertAlmostEqual(_angle(h_on_c2[0], c2r, c1r), 109.47, delta=1.0)
        self.assertAlmostEqual(_angle(h_on_c2[0], c2r, c3r), 109.47, delta=1.0)
        self.assertAlmostEqual(_angle(h_on_c2[0], c2r, h_on_c2[1]), 109.5, delta=0.5)

    def test_sp2_aromatic_carbon_gets_1_h_in_plane(self):
        # Aromatic-like ring carbon (atom 2) with 2 ring neighbors
        # (atoms 1 and 3) at 120 deg apart. Index 0 is unused padding.
        mol = _make_mol(3, [None, 'C ', 'C ', 'C '], [0, 6, 6, 6],
                         [[0, 0, 0], [0, 0, 0], [1.4, 0, 0], [2.1, 1.2, 0]],
                         [(1, 2, 2), (2, 3, 1)])
        result = add_missing_hydrogens(mol)
        c1, c2, c3 = result['coords'][1], result['coords'][2], result['coords'][3]
        h_atoms = [result['coords'][i] for i in range(1, result['natom'] + 1)
                   if result['nuc'][i] == 1]
        # Only atom 2 needs an H (1 missing valence); atom 1 and atom 3
        # are themselves only 1-coordinate carbons in this toy fragment
        # and would also want H's - so check specifically the H attached
        # to atom 2, identified by its idcon back-reference.
        h_on_c2 = [result['coords'][i] for i in range(1, result['natom'] + 1)
                   if result['nuc'][i] == 1 and result['idcon'][i][1] == 2]
        self.assertEqual(len(h_on_c2), 1)
        self.assertAlmostEqual(_angle(h_on_c2[0], c2, c1), 120, delta=1.0)
        self.assertAlmostEqual(_angle(h_on_c2[0], c2, c3), 120, delta=1.0)

    def test_methane_carbon_count_matches_formula(self):
        # Single N with no bonds should get 3 H's (standard valence 3)
        mol = _make_mol(1, [None, 'N '], [0, 7], [[0, 0, 0], [0, 0, 0]], [])
        result = add_missing_hydrogens(mol)
        self.assertEqual(result['hydrogens_added'], 3)

    def test_oxygen_gets_2_h(self):
        mol = _make_mol(1, [None, 'O '], [0, 8], [[0, 0, 0], [0, 0, 0]], [])
        result = add_missing_hydrogens(mol)
        self.assertEqual(result['hydrogens_added'], 2)

    def test_double_bonded_carbon_uses_remaining_valence(self):
        # C=C, each carbon already has bond order 2 used, needs 2 more (sp2 alkene)
        mol = _make_mol(2, [None, 'C ', 'C '], [0, 6, 6],
                         [[0, 0, 0], [0, 0, 0], [1.34, 0, 0]], [(1, 2, 2)])
        result = add_missing_hydrogens(mol)
        # Each carbon: used=2 (double bond), target=4, needs 2 H's each = 4 total
        self.assertEqual(result['hydrogens_added'], 4)

    def test_chlorine_gets_no_extra_h_if_already_bonded(self):
        mol = _make_mol(2, [None, 'C ', 'Cl'], [0, 6, 17],
                         [[0, 0, 0], [0, 0, 0], [1.7, 0, 0]], [(1, 2, 1)])
        result = add_missing_hydrogens(mol)
        # Cl already has its single bond satisfied; only C needs H's (3 more)
        self.assertEqual(result['hydrogens_added'], 3)

    def test_no_bond_length_distortion(self):
        mol = _make_mol(1, [None, 'C '], [0, 6], [[0, 0, 0], [0, 0, 0]], [])
        result = add_missing_hydrogens(mol, bond_length=1.09)
        c = result['coords'][1]
        for i in range(2, result['natom'] + 1):
            self.assertAlmostEqual(_bond_len(c, result['coords'][i]), 1.09, places=4)


class TestNeedsHydrogens(unittest.TestCase):
    def test_small_molecule_no_warning(self):
        # 5 atoms, no H's - below the >10 atom threshold, should NOT trigger
        mol = _make_mol(5, [None] + ['C '] * 5, [0] + [6] * 5,
                         [[0, 0, 0]] * 6, [(1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1)])
        self.assertFalse(needs_hydrogens(mol))

    def test_large_molecule_no_h_triggers(self):
        n = 12
        atoms = [None] + ['C '] * n
        nuc = [0] + [6] * n
        coords = [[0, 0, 0]] * (n + 1)
        bonds = [(i, i + 1, 1) for i in range(1, n)]
        mol = _make_mol(n, atoms, nuc, coords, bonds)
        self.assertTrue(needs_hydrogens(mol))

    def test_large_molecule_with_h_does_not_trigger(self):
        n = 12
        atoms = [None] + ['C '] * (n - 1) + ['H ']
        nuc = [0] + [6] * (n - 1) + [1]
        coords = [[0, 0, 0]] * (n + 1)
        bonds = [(i, i + 1, 1) for i in range(1, n)]
        mol = _make_mol(n, atoms, nuc, coords, bonds)
        self.assertFalse(needs_hydrogens(mol))


class TestRealMoleculeFormulaConsistency(unittest.TestCase):
    """Validate H-adder against real GAC molecules with known formulas."""
    GAC_SDF = _GAC_SDF

    def _strip_hydrogens(self, mol):
        heavy_idx = [i for i in range(1, mol['natom'] + 1) if mol['nuc'][i] != 1]
        remap = {old: new + 1 for new, old in enumerate(heavy_idx)}
        n2 = len(heavy_idx)
        atom2 = [None] * (n2 + 1)
        nuc2 = [0] * (n2 + 1)
        coords2 = [[0, 0, 0]] * (n2 + 1)
        idcon2 = [[0] * 5 for _ in range(n2 + 1)]
        ival2 = [0] * (n2 + 1)
        nbo2 = SparseMatrix()
        for old in heavy_idx:
            new = remap[old]
            atom2[new] = mol['atom'][old]
            nuc2[new] = mol['nuc'][old]
            coords2[new] = mol['coords'][old]
            cnt = 0
            for j in range(1, mol['ival'][old] + 1):
                nbr_old = mol['idcon'][old][j]
                if nbr_old in remap:
                    cnt += 1
                    idcon2[new][cnt] = remap[nbr_old]
                    bo = mol['nbo'][old][nbr_old]
                    nbo2[new][remap[nbr_old]] = bo if bo else 1
            ival2[new] = cnt
        icon2 = [[0] * 5 for _ in range(n2 + 1)]
        isum2 = [0] * (n2 + 1)
        for i in range(1, n2 + 1):
            s = 0
            for j in range(1, ival2[i] + 1):
                nb = nuc2[idcon2[i][j]]
                icon2[i][j] = nb
                s += nb
            isum2[i] = s
        return dict(natom=n2, numbond=0, atom=atom2, nuc=nuc2, coords=coords2,
                    idcon=idcon2, ival=ival2, nbo=nbo2, icon=icon2, isum=isum2,
                    charge={})

    def test_h_count_matches_molecular_formula(self):
        if not os.path.exists(self.GAC_SDF):
            self.skipTest("GAC SDF not available")
        import re
        with open(self.GAC_SDF) as fh:
            from recon6.readers.sdf import read_sdf
            checked = 0
            idx = 0
            while checked < 30 and idx < 200:
                try:
                    mol = read_sdf(fh)
                except StopIteration:
                    break
                idx += 1
                name = mol.get('mol_name', '')
                m = re.search(r'H(\d+)', name)
                if not m:
                    continue
                expected_h = int(m.group(1))
                stripped = self._strip_hydrogens(mol)
                if not needs_hydrogens(stripped, threshold_atoms=5):
                    continue
                result = add_missing_hydrogens(stripped)
                self.assertEqual(
                    result['hydrogens_added'], expected_h,
                    "Formula %s implies %d H but got %d" % (
                        name, expected_h, result['hydrogens_added'])
                )
                checked += 1
            self.assertGreater(checked, 0, "No suitable test molecules found")


if __name__ == '__main__':
    unittest.main()
