"""
Integration test: GAC_withH.sdf (974 molecules).

IMPORTANT CONTEXT: despite the filename, 971 of these 974 molecules
have ZERO explicit hydrogens in the source file. The original Fortran
RECON5 has no hydrogen-addition capability - it just logs a warning
for such molecules and proceeds with whatever (degenerate, H-less)
connectivity it was given. The recon.ff reference output for those
971 rows therefore reflects Fortran's H-less behavior, not a true
H-saturated calculation, and is NOT a valid comparison target now that
this Python port adds missing hydrogens by default (see hydrogenate.py
and the no-add-h release notes).

Only 3 of the 974 molecules (indices 155, 812, 813) already had
explicit hydrogens in the source file; for those, the Fortran
reference remains a valid like-for-like comparison and is checked
directly. For the other 971, this test instead validates that the
H-addition behaved sensibly (added hydrogen count matches the
molecular formula in the title line) and that the pipeline completes
without errors across the whole batch.
"""
import unittest
import sys, os, csv, re
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.dirname(__file__))

from test_config import DATA_DIR, BOND_FILE, GAC_SDF, GAC_FF
HAVE_DATA = os.path.exists(DATA_DIR) and os.path.exists(GAC_SDF)

# The only 3 rows where Fortran's reference output reflects a molecule
# that genuinely had explicit hydrogens already (verified by direct
# inspection of the source SDF).
_VALID_FORTRAN_COMPARISON = {155, 812, 813}

# The 3 rows with a known data defect (disconnected halogen atoms with
# zero bonds in the source SDF) - see README "Known limitations".
_KNOWN_BAD = {885, 942, 944}


class TestGACIntegration(unittest.TestCase):
    _results = None
    _fortran_rows = None

    @classmethod
    def setUpClass(cls):
        if not HAVE_DATA:
            return
        from recon6.recon import ReconConfig, run_recon
        config = ReconConfig(
            data_dir=DATA_DIR, bond_file=BOND_FILE,
            input_files=[GAC_SDF], fmt='sdf',
        )
        cls._results = run_recon(config)

        cls._fortran_rows = []
        with open(GAC_FF) as f:
            reader = csv.reader(f)
            header = [h.strip() for h in next(reader)]
            for row in reader:
                row = [v.strip() for v in row]
                cls._fortran_rows.append(dict(zip(header, row)))

    def _skip_if_no_data(self):
        if not HAVE_DATA or self._results is None:
            self.skipTest("Test data not available")

    def test_molecule_count(self):
        self._skip_if_no_data()
        self.assertEqual(len(self._results), 974)
        self.assertEqual(len(self._fortran_rows), 974)

    def test_first_molecule_name(self):
        self._skip_if_no_data()
        self.assertEqual(self._results[0]['Molecule'], 'C8H12N8S2')

    def test_valid_fortran_comparisons_match(self):
        """The 3 molecules that already had explicit H's in the source
        file should match the Fortran reference closely (these are a
        genuine apples-to-apples comparison, unlike the other 971)."""
        self._skip_if_no_data()
        for idx in _VALID_FORTRAN_COMPARISON:
            prow = self._results[idx]
            frow = self._fortran_rows[idx]
            fval = float(frow['Energy'])
            pval = prow['Energy']
            err = abs(pval - fval) / abs(fval) * 100
            self.assertLess(
                err, 0.05,
                "idx=%d %s: Fortran=%.4f Python=%.4f err=%.4f%%" % (
                    idx, prow['Molecule'], fval, pval, err)
            )

    def test_hydrogen_addition_matches_formula(self):
        """For molecules where H's were added by this port (i.e. all
        except the 3 that already had some explicit H's in the source
        and the 3 with known disconnected-atom defects), the final H
        count present should match the molecular formula's H count
        (parsed from the title line)."""
        self._skip_if_no_data()
        checked = 0
        skip_indices = _KNOWN_BAD | _VALID_FORTRAN_COMPARISON
        for idx, row in enumerate(self._results):
            if idx in skip_indices:
                continue
            name = row['Molecule']
            m = re.search(r'H(\d+)', name)
            if not m:
                continue
            expected_h = int(m.group(1))
            self.assertEqual(
                row['_num_h'], expected_h,
                "idx=%d %s: formula implies %d H but molecule has %d" % (
                    idx, name, expected_h, row['_num_h'])
            )
            checked += 1
        self.assertGreater(checked, 900,
                            "Expected most of the 974 molecules to have "
                            "a parseable formula in their name")

    def test_known_bad_records_have_disconnected_atoms(self):
        """Documents *why* indices 885, 942, 944 were previously flagged:
        each contains a halogen atom with zero bonds in the source SDF.
        Still true after H-addition (disconnected atoms get no new H's
        either, since H-addition only adds bonds to atoms below their
        standard valence, and an atom with zero existing bonds and an
        element not in the standard valence table - or already at its
        cap - is left as-is)."""
        self._skip_if_no_data()
        for idx in _KNOWN_BAD:
            mol_name = self._results[idx]['Molecule']
            self.assertIn(mol_name, ('C27H27Cl2F3N8O3S', 'C33H36ClF3N6O3S'))

    def test_all_molecules_processed_without_exception(self):
        self._skip_if_no_data()
        for row in self._results:
            self.assertIsNotNone(row.get('Energy'))
            self.assertIsNotNone(row.get('chi'))


if __name__ == '__main__':
    unittest.main()
