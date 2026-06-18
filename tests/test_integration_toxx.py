"""Integration test: batch validation against toxx2.sdf (278 molecules)."""
import unittest
import sys, os, csv
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.dirname(__file__))

from test_config import DATA_DIR, BOND_FILE, TOXX_SDF, TOXX_FF
HAVE_DATA = os.path.exists(DATA_DIR) and os.path.exists(TOXX_SDF)


class TestToxxBatch(unittest.TestCase):
    # Lists (not dicts) to preserve order and duplicate molecule names.
    _results = None
    _fortran_rows = None

    @classmethod
    def setUpClass(cls):
        if not HAVE_DATA:
            return
        from recon6.recon import ReconConfig, run_recon
        config = ReconConfig(
            data_dir=DATA_DIR, bond_file=BOND_FILE,
            input_files=[TOXX_SDF], fmt='sdf',
        )
        cls._results = run_recon(config)

        cls._fortran_rows = []
        with open(TOXX_FF) as f:
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
        self.assertEqual(len(self._results), 278)
        self.assertEqual(len(self._fortran_rows), 278)

    def test_molecule_names_match_in_order(self):
        self._skip_if_no_data()
        for prow, frow in zip(self._results, self._fortran_rows):
            self.assertEqual(prow['Molecule'], frow['Molecule'])

    def _max_rel_err(self, field):
        max_err = 0.0
        for prow, frow in zip(self._results, self._fortran_rows):
            try:
                fval = float(frow[field])
            except (KeyError, ValueError):
                continue
            pval = prow.get(field)
            if pval is None or fval == 0:
                continue
            err = abs(pval - fval) / abs(fval) * 100
            max_err = max(max_err, err)
        return max_err

    def test_energy_accuracy_batch(self):
        self._skip_if_no_data()
        self.assertLess(self._max_rel_err('Energy'), 0.05)

    def test_sik_accuracy_batch(self):
        self._skip_if_no_data()
        self.assertLess(self._max_rel_err('SIK'), 0.05)

    def test_surfarea_accuracy_batch(self):
        self._skip_if_no_data()
        self.assertLess(self._max_rel_err('SurfArea'), 0.05)


if __name__ == '__main__':
    unittest.main()
