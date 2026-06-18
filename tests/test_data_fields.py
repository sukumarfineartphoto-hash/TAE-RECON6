"""Tests for SDF '> <TAG>' data-field extraction and the positional
alignment guarantee in run_recon (rows preserved 1:1 with input order,
even for unparseable/skipped molecules)."""
import unittest
import sys, os, io, csv
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.dirname(__file__))

from recon6.readers.sdf import read_sdf, _read_molecule_block, _parse_molecule_block
from test_config import DATA_DIR, BOND_FILE, TOXX_SDF, TOXX_FF, TOXX_MOL2, HAVE_DATA


_SAMPLE_SDF = """mol_A
  -ISIS-

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <RESPONSE>
1.23

$$$$
mol_B_broken
  -ISIS-

 BADCOUNT  0  0  0  0  0  0  0  0999 V2000
M  END
> <RESPONSE>
4.56

$$$$
mol_C
  -ISIS-

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <RESPONSE>
7.89

$$$$
"""


class TestDataFieldExtraction(unittest.TestCase):
    def test_single_tag_extracted(self):
        sdf = """mol1
  prog

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <LIC50>
0.085

$$$$
"""
        mol = read_sdf(io.StringIO(sdf))
        self.assertEqual(mol['data_fields'], {'LIC50': '0.085'})

    def test_multiple_tags_extracted_in_order(self):
        sdf = """mol1
  prog

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
M  END
>  <EXTREG>  (some-name)
some-name

> <LIC50>
0.085

$$$$
"""
        mol = read_sdf(io.StringIO(sdf))
        self.assertEqual(mol['data_fields']['EXTREG'], 'some-name')
        self.assertEqual(mol['data_fields']['LIC50'], '0.085')
        self.assertEqual(list(mol['data_fields'].keys()), ['EXTREG', 'LIC50'])

    def test_no_data_fields_when_none_present(self):
        sdf = """mol1
  prog

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""
        mol = read_sdf(io.StringIO(sdf))
        self.assertEqual(mol['data_fields'], {})

    def test_real_toxx_file_field_values(self):
        path = TOXX_SDF
        if not os.path.exists(path):
            self.skipTest("toxx2.sdf not available")
        with open(path) as fh:
            mol1 = read_sdf(fh)
            mol2 = read_sdf(fh)
        self.assertEqual(mol1['data_fields']['EXTREG'], '3-hydroxybenzaldehyde')
        self.assertEqual(mol1['data_fields']['LIC50'], '0.085')
        self.assertEqual(mol2['data_fields']['EXTREG'], 'nonylphenol')
        self.assertEqual(mol2['data_fields']['LIC50'], '2.468')


class TestMalformedRecordRecovery(unittest.TestCase):
    """A malformed molecule block must not desync subsequent reads, and
    its data_fields should still be salvageable even though structural
    parsing fails."""

    def test_malformed_block_does_not_desync_file_position(self):
        fh = io.StringIO(_SAMPLE_SDF)
        block1 = _read_molecule_block(fh)
        mol1 = _parse_molecule_block(block1)
        self.assertEqual(mol1['mol_name'], 'mol_A')

        block2 = _read_molecule_block(fh)
        with self.assertRaises(Exception):
            _parse_molecule_block(block2)

        # Despite the failure above, the third block should still be
        # readable correctly - i.e. file position wasn't corrupted.
        block3 = _read_molecule_block(fh)
        mol3 = _parse_molecule_block(block3)
        self.assertEqual(mol3['mol_name'], 'mol_C')

    def test_data_fields_salvaged_from_broken_block(self):
        from recon6.readers.sdf import read_sdf_data_fields_only
        fh = io.StringIO(_SAMPLE_SDF)
        _read_molecule_block(fh)  # mol_A, skip
        block2 = _read_molecule_block(fh)  # mol_B_broken
        name, data_fields = read_sdf_data_fields_only(block2)
        self.assertEqual(name, 'mol_B_broken')
        self.assertEqual(data_fields, {'RESPONSE': '4.56'})


class TestPositionalAlignment(unittest.TestCase):
    """End-to-end: a molecule that fails to parse/process must be
    dropped from the output entirely (no blank row), while the
    surviving molecules' data fields must remain correctly matched to
    the right molecule."""

    def setUp(self):
        if not HAVE_DATA:
            self.skipTest("Test data not available")
        self.sdf_path = '/tmp/test_alignment_sample.sdf'
        with open(self.sdf_path, 'w') as f:
            f.write(_SAMPLE_SDF)

    def test_failed_molecule_is_dropped_not_blank(self):
        from recon6.recon import ReconConfig, run_recon
        config = ReconConfig(
            data_dir=DATA_DIR, bond_file=BOND_FILE,
            input_files=[self.sdf_path], fmt='sdf', auto_add_h=False,
        )
        results = run_recon(config)
        # 3 input molecules, 1 broken -> only 2 rows in the output
        self.assertEqual(len(results), 2)
        names = [r['Molecule'] for r in results]
        self.assertEqual(names, ['mol_A', 'mol_C'])

    def test_response_data_stays_matched_to_correct_molecule(self):
        from recon6.recon import ReconConfig, run_recon
        config = ReconConfig(
            data_dir=DATA_DIR, bond_file=BOND_FILE,
            input_files=[self.sdf_path], fmt='sdf', auto_add_h=False,
        )
        results = run_recon(config)
        by_name = {r['Molecule']: r for r in results}
        self.assertEqual(by_name['mol_A']['RESPONSE'], '1.23')
        self.assertEqual(by_name['mol_C']['RESPONSE'], '7.89')
        self.assertNotIn('mol_B_broken', by_name)

    def test_csv_output_has_no_blank_descriptor_rows(self):
        from recon6.recon import ReconConfig, run_recon
        out_path = '/tmp/test_alignment_output.csv'
        config = ReconConfig(
            data_dir=DATA_DIR, bond_file=BOND_FILE,
            input_files=[self.sdf_path], fmt='sdf', auto_add_h=False,
            output_csv=out_path,
        )
        run_recon(config)
        csv.field_size_limit(10_000_000)
        with open(out_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        self.assertEqual(len(rows), 2)
        self.assertIn('RESPONSE', reader.fieldnames)
        self.assertNotIn('Skipped', reader.fieldnames)
        for row in rows:
            self.assertNotEqual(row['Energy'], '')

    def test_atom_records_not_leaked_into_csv(self):
        from recon6.recon import ReconConfig, run_recon
        out_path = '/tmp/test_no_leak_output.csv'
        config = ReconConfig(
            data_dir=DATA_DIR, bond_file=BOND_FILE,
            input_files=[self.sdf_path], fmt='sdf', auto_add_h=False,
            output_csv=out_path,
        )
        run_recon(config)
        with open(out_path) as f:
            header = f.readline()
        self.assertNotIn('atom_records', header)


class TestCsvQuotingSafety(unittest.TestCase):
    """Molecule names or data fields containing commas/quotes must be
    safely quoted in the output CSV (and round-trip correctly), since
    these come directly from arbitrary user-authored SDF text."""

    def setUp(self):
        if not HAVE_DATA:
            self.skipTest("Test data not available")

    def test_comma_in_title_and_data_field_round_trips(self):
        from recon6.recon import ReconConfig, run_recon
        sdf = (
            "compound, with comma in name\n"
            "  prog\n\n"
            "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "M  END\n"
            "> <NOTE>\n"
            "value, with, commas\n\n"
            "$$$$\n"
        )
        path = '/tmp/test_comma_quoting.sdf'
        with open(path, 'w') as f:
            f.write(sdf)
        out_path = '/tmp/test_comma_quoting_out.csv'
        config = ReconConfig(
            data_dir=DATA_DIR, bond_file=BOND_FILE,
            input_files=[path], fmt='sdf', auto_add_h=False,
            output_csv=out_path,
        )
        run_recon(config)
        csv.field_size_limit(10_000_000)
        with open(out_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)
        self.assertEqual(row['Molecule'], 'compound, with comma in name')
        self.assertEqual(row['NOTE'], 'value, with, commas')

    def test_double_quote_in_data_field_round_trips(self):
        from recon6.recon import ReconConfig, run_recon
        sdf = (
            'mol_with_quotes\n'
            "  prog\n\n"
            "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "M  END\n"
            "> <NOTE>\n"
            'has "quoted" text\n\n'
            "$$$$\n"
        )
        path = '/tmp/test_quote_quoting.sdf'
        with open(path, 'w') as f:
            f.write(sdf)
        out_path = '/tmp/test_quote_quoting_out.csv'
        config = ReconConfig(
            data_dir=DATA_DIR, bond_file=BOND_FILE,
            input_files=[path], fmt='sdf', auto_add_h=False,
            output_csv=out_path,
        )
        run_recon(config)
        csv.field_size_limit(10_000_000)
        with open(out_path) as f:
            reader = csv.DictReader(f)
            row = next(reader)
        self.assertEqual(row['NOTE'], 'has "quoted" text')

    def test_two_molecules_one_with_comma_one_without(self):
        """Regression guard: a comma in one row's name must not shift
        columns for that row or any other row in the same file."""
        from recon6.recon import ReconConfig, run_recon
        sdf = (
            "plain_name\n"
            "  prog\n\n"
            "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "M  END\n"
            "> <NOTE>\nplain\n\n"
            "$$$$\n"
            "name, with comma\n"
            "  prog\n\n"
            "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "M  END\n"
            "> <NOTE>\ncomma, here\n\n"
            "$$$$\n"
        )
        path = '/tmp/test_two_mol_comma.sdf'
        with open(path, 'w') as f:
            f.write(sdf)
        out_path = '/tmp/test_two_mol_comma_out.csv'
        config = ReconConfig(
            data_dir=DATA_DIR, bond_file=BOND_FILE,
            input_files=[path], fmt='sdf', auto_add_h=False,
            output_csv=out_path,
        )
        run_recon(config)
        csv.field_size_limit(10_000_000)
        with open(out_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]['Molecule'], 'plain_name')
        self.assertEqual(rows[0]['NOTE'], 'plain')
        self.assertEqual(rows[1]['Molecule'], 'name, with comma')
        self.assertEqual(rows[1]['NOTE'], 'comma, here')
        # Column count must be identical and well-formed for both rows
        self.assertEqual(len(rows[0]), len(rows[1]))


class TestImplicitBondFileResolution(unittest.TestCase):
    """bond_file is now optional; when omitted, it should resolve to
    '<data_dir>/bond' automatically."""

    def test_resolved_bond_file_defaults_to_data_dir_bond(self):
        from recon6.recon import ReconConfig
        config = ReconConfig(data_dir='/some/path', input_files=[])
        self.assertEqual(config.resolved_bond_file(), '/some/path/bond')

    def test_explicit_bond_file_overrides_default(self):
        from recon6.recon import ReconConfig
        config = ReconConfig(data_dir='/some/path', input_files=[],
                              bond_file='/other/path/custom_bond')
        self.assertEqual(config.resolved_bond_file(), '/other/path/custom_bond')

    def test_run_recon_works_without_explicit_bond_file(self):
        if not HAVE_DATA:
            self.skipTest("Test data not available")
        from recon6.recon import ReconConfig, run_recon
        config = ReconConfig(
            data_dir=DATA_DIR,
            input_files=[TOXX_SDF],
            fmt='sdf',
        )
        results = run_recon(config)
        self.assertGreater(len(results), 0)


class TestMol2BatchValidation(unittest.TestCase):
    """toxx.mol2 (the SYBYL-converted version of the toxx2.sdf test
    set, 278 molecules) is NOT in the same order as toxx2.sdf - the
    SDF->MOL2 conversion reordered molecules - so this validates by
    comparing the *multiset* of computed Energy values against the
    Fortran reference rather than assuming positional alignment."""

    MOL2_PATH = TOXX_MOL2

    def setUp(self):
        if not HAVE_DATA or not os.path.exists(self.MOL2_PATH):
            self.skipTest("toxx.mol2 or DATA dir not available")

    def test_all_molecules_match_fortran_energy_set(self):
        from recon6.recon import ReconConfig, run_recon
        from collections import Counter

        config = ReconConfig(
            data_dir=DATA_DIR, input_files=[self.MOL2_PATH], fmt='mol2',
        )
        results = run_recon(config)
        self.assertEqual(len(results), 278)

        mol2_energies = Counter(round(r['Energy'], 1) for r in results)

        fortran_energies = Counter()
        with open(TOXX_FF) as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                fortran_energies[round(float(row[1]), 1)] += 1

        matched = sum((mol2_energies & fortran_energies).values())
        self.assertEqual(matched, 278,
                          "every MOL2-derived energy should match some "
                          "Fortran reference energy (as an unordered set)")


if __name__ == '__main__':
    unittest.main()
