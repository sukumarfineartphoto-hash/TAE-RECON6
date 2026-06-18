"""Tests for the PDB reader: element-column priority, deuterium
handling, distance-based connectivity fallback, bond-order inference
for H-counting, and integration with the H-adder."""
import unittest
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.dirname(__file__))

from recon6.readers.pdb import read_pdb, _element_from_pdb_line
from recon6.bonds import load_bond_table
from test_config import DATA_DIR, BOND_FILE, PDB_NNA as NNA_PDB, PDB_SV6 as SV6_PDB, HAVE_DATA

HAVE_PDB_SAMPLES = os.path.exists(NNA_PDB) and os.path.exists(SV6_PDB)


class TestElementColumnExtraction(unittest.TestCase):
    def test_element_column_used_when_present(self):
        line = "HETATM 2388  C1  NNA A 401     125.180   4.744  22.565  1.00 28.56           C  \n"
        self.assertEqual(_element_from_pdb_line(line, 'C1'), 'C ')

    def test_deuterium_mapped_to_hydrogen(self):
        line = "HETATM 6449 DNAD SV6 A 401      11.405   6.193  26.179  1.00 39.74           D  \n"
        # Atom name 'DNAD' would NOT be recognized correctly by the
        # name-based heuristic (ele.f has no rule for names starting
        # with 'D' followed by other letters), but the element column
        # correctly identifies it as deuterium, which we fold into H.
        self.assertEqual(_element_from_pdb_line(line, 'DNAD'), 'H ')

    def test_falls_back_to_name_heuristic_when_no_element_column(self):
        short_line = "HETATM 2388  C1  NNA A 401     125.180   4.744  22.565\n"
        result = _element_from_pdb_line(short_line, 'C1')
        self.assertEqual(result, 'C ')

    def test_sulfur_element_column(self):
        line = "HETATM 2398  S11 NNA A 401     126.012   8.398  23.175  1.00 37.65           S  \n"
        self.assertEqual(_element_from_pdb_line(line, 'S11'), 'S ')


@unittest.skipUnless(HAVE_DATA and HAVE_PDB_SAMPLES, "PDB sample files or DATA dir not available")
class TestRealPdbConnectivity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bondl = load_bond_table(BOND_FILE)

    def test_nna_no_conect_falls_back_to_distance(self):
        """NNA-7jyc.pdb has zero CONECT records; connectivity must
        still be derived (via distance fallback) rather than left
        empty."""
        mol = read_pdb(NNA_PDB, bondl=self.bondl, iovr=-1, auto_add_h=False)
        self.assertEqual(mol['natom'], 49)
        total_bonds_per_atom = sum(mol['ival'][i] for i in range(1, mol['natom'] + 1))
        self.assertGreater(total_bonds_per_atom, 0)

    def test_nna_no_valence_violations_before_h(self):
        mol = read_pdb(NNA_PDB, bondl=self.bondl, iovr=-1, auto_add_h=False)
        maxval = {6: 4, 7: 3, 8: 2, 16: 6}
        for i in range(1, mol['natom'] + 1):
            nuc = mol['nuc'][i]
            self.assertLessEqual(mol['ival'][i], maxval.get(nuc, 4),
                                  "atom %d (nuc=%d) exceeds standard valence" % (i, nuc))

    def test_nna_deuterium_or_unknown_atoms_none(self):
        mol = read_pdb(NNA_PDB, bondl=self.bondl, iovr=-1, auto_add_h=False)
        unknown = [i for i in range(1, mol['natom'] + 1) if mol['nuc'][i] == 0]
        self.assertEqual(unknown, [])

    def test_sv6_deuterium_atoms_recognized_as_hydrogen(self):
        """SV6-7lb7.pdb has 5 atoms with element column 'D' and atom
        names like 'DNAD' that the name-based heuristic alone would
        misclassify; all 5 must resolve to nuc=1 (hydrogen), not 0."""
        mol = read_pdb(SV6_PDB, bondl=self.bondl, iovr=-1, auto_add_h=False)
        unknown = [i for i in range(1, mol['natom'] + 1) if mol['nuc'][i] == 0]
        self.assertEqual(unknown, [],
                          "deuterium atoms should resolve to H, not unknown element")


@unittest.skipUnless(HAVE_DATA and HAVE_PDB_SAMPLES, "PDB sample files or DATA dir not available")
class TestBondOrderInference(unittest.TestCase):
    """Carbonyl/imine oxygens and nitrogens are systematically shorter
    than the single-bond reference length; this must be detected and
    used to avoid over-counting hydrogens on those centers."""

    @classmethod
    def setUpClass(cls):
        cls.bondl = load_bond_table(BOND_FILE)

    def test_nna_h_count_matches_formula(self):
        """NNA-7jyc ligand formula is C36H63N5O7S (63 H). Without
        bond-order inference this came out to 73 (10 too many, all on
        carbonyl-type centers); with it, it should match exactly."""
        mol = read_pdb(NNA_PDB, bondl=self.bondl, iovr=-1, auto_add_h=True)
        self.assertEqual(mol['hydrogens_added'], 63)

    def test_carbonyl_oxygen_gets_no_extra_hydrogen(self):
        """The first carbonyl oxygen (atom 2, bonded to atom 1) should
        not gain a hydrogen - its short bond length should be detected
        as a double bond, satisfying its valence without H."""
        mol = read_pdb(NNA_PDB, bondl=self.bondl, iovr=-1, auto_add_h=False)
        from recon6.hydrogenate import _bond_order_sum, _best_standard_valence
        used = _bond_order_sum(mol, 2)
        target = _best_standard_valence(mol['nuc'][2], used)
        self.assertEqual(used, target, "carbonyl O should already satisfy its valence")

    def test_hydroxyl_oxygen_still_gets_hydrogen(self):
        """Atom 42 (O) - C43 is a genuine single-bonded hydroxyl-like
        oxygen (bond length ratio ~0.98, near the single-bond
        reference) and should still correctly receive a hydrogen,
        i.e. the inference shouldn't over-fire on ordinary single
        bonds that happen to be slightly short."""
        mol = read_pdb(NNA_PDB, bondl=self.bondl, iovr=-1, auto_add_h=False)
        from recon6.hydrogenate import _bond_order_sum, _best_standard_valence
        used = _bond_order_sum(mol, 42)
        target = _best_standard_valence(mol['nuc'][42], used)
        self.assertEqual(target - used, 1, "hydroxyl-like O should still need 1 H")

    def test_sulfonyl_oxygens_get_no_extra_hydrogen(self):
        """Atoms 12 and 13 are the two oxygens of a sulfonyl (-SO2-)
        group; S=O bonds are only modestly shortened relative to S-O
        single bonds (smaller contraction than C=O), so this needs the
        dedicated sulfonyl special case, not just the general ratio
        threshold."""
        mol = read_pdb(NNA_PDB, bondl=self.bondl, iovr=-1, auto_add_h=False)
        from recon6.hydrogenate import _bond_order_sum, _best_standard_valence
        for idx in (12, 13):
            used = _bond_order_sum(mol, idx)
            target = _best_standard_valence(mol['nuc'][idx], used)
            self.assertEqual(used, target,
                              "sulfonyl O (atom %d) should already satisfy its valence" % idx)


if __name__ == '__main__':
    unittest.main()
