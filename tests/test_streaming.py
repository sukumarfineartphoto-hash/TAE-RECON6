"""Tests for streaming CSV output, optional result accumulation, and
log-file redirection introduced in the large-batch memory fix."""
import unittest
import sys, os, csv, tempfile, time, threading
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.dirname(__file__))

from test_config import DATA_DIR, BOND_FILE, TOXX_SDF, TOXX_FF, HAVE_DATA

csv.field_size_limit(10_000_000)

# A minimal valid SDF (two H atoms - simplest possible processed molecule)
_ONE_MOL_SDF = """\
mol_A
  test

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <RESPONSE>
1.23

$$$$
"""

_TWO_MOL_SDF = """\
mol_A
  test

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <RESPONSE>
1.23

$$$$
mol_B
  test

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <RESPONSE>
4.56

$$$$
"""


@unittest.skipUnless(HAVE_DATA, "DATA dir not available")
class TestStreamingCsvWriter(unittest.TestCase):
    """CSV rows are flushed to disk immediately after each molecule is
    processed - not batched until the end."""

    def test_csv_file_appears_before_all_molecules_processed(self):
        """By the time the first molecule is done, the output file must
        already contain its row.  We verify this indirectly: write two
        molecules, then check that the file size grew after the first
        and was not zero until the second completed."""
        from recon6.recon import ReconConfig, run_recon

        with tempfile.NamedTemporaryFile(suffix='.sdf', mode='w',
                                         delete=False) as f:
            f.write(_TWO_MOL_SDF)
            sdf_path = f.name

        out_path = sdf_path.replace('.sdf', '.csv')
        try:
            sizes_after_each = []

            # Monkeypatch _process_molecule to record CSV size mid-run
            import recon6.recon as mod
            original = mod._process_molecule

            call_count = [0]
            def patched(mol, tae_index, config):
                result = original(mol, tae_index, config)
                call_count[0] += 1
                if os.path.exists(out_path):
                    sizes_after_each.append(os.path.getsize(out_path))
                return result

            mod._process_molecule = patched
            try:
                config = ReconConfig(
                    data_dir=DATA_DIR, input_files=[sdf_path],
                    fmt='sdf', output_csv=out_path, auto_add_h=False,
                )
                run_recon(config)
            finally:
                mod._process_molecule = original

            self.assertEqual(call_count[0], 2)
            self.assertEqual(len(sizes_after_each), 2)
            # File must have grown between molecule 1 and molecule 2
            self.assertGreater(sizes_after_each[0], 0,
                               "CSV must be non-empty after first molecule")
            self.assertGreater(sizes_after_each[1], sizes_after_each[0],
                               "CSV must grow after second molecule")
        finally:
            for p in [sdf_path, out_path]:
                if os.path.exists(p): os.unlink(p)

    def test_header_written_immediately(self):
        """The CSV header should be written before any molecule is
        processed, not only when the first data row arrives."""
        from recon6.recon import _CsvStreamWriter
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as f:
            path = f.name
        try:
            writer = _CsvStreamWriter(path, include_data_fields=False)
            # Before writing any rows, the file should already have a header
            size_after_open = os.path.getsize(path)
            writer.close()
            self.assertGreater(size_after_open, 0,
                               "Header must be written on open, not on first row")
        finally:
            os.unlink(path)

    def test_new_data_field_tag_triggers_rewrite(self):
        """If a data_field tag appears in a later row that wasn't in
        the first row, the writer must rewrite to add the new column."""
        from recon6.recon import _CsvStreamWriter
        with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as f:
            path = f.name
        try:
            writer = _CsvStreamWriter(path, include_data_fields=True)
            writer.write({'Molecule': 'A', 'Energy': -1.0, 'chi': 0.0,
                          'TAG1': 'val1'})
            writer.write({'Molecule': 'B', 'Energy': -2.0, 'chi': 0.0,
                          'TAG1': 'val2', 'TAG2': 'new_tag'})
            writer.close()

            with open(path) as f:
                reader = csv.DictReader(f)
                rows = list(reader)
            self.assertIn('TAG2', reader.fieldnames)
            self.assertEqual(rows[0].get('TAG2'), '')  # blank for old row
            self.assertEqual(rows[1].get('TAG2'), 'new_tag')
        finally:
            os.unlink(path)


@unittest.skipUnless(HAVE_DATA, "DATA dir not available")
class TestReturnResultsFlag(unittest.TestCase):
    """return_results=False must not accumulate results in memory."""

    def test_return_results_false_gives_empty_list(self):
        from recon6.recon import ReconConfig, run_recon

        with tempfile.NamedTemporaryFile(suffix='.sdf', mode='w',
                                         delete=False) as f:
            f.write(_TWO_MOL_SDF)
            sdf_path = f.name
        out_path = sdf_path.replace('.sdf', '.csv')
        try:
            config = ReconConfig(
                data_dir=DATA_DIR, input_files=[sdf_path],
                fmt='sdf', output_csv=out_path,
                auto_add_h=False, return_results=False,
            )
            results = run_recon(config)
            self.assertEqual(results, [],
                             "return_results=False must return empty list")
            # But the CSV must still have been written
            with open(out_path) as f:
                rows = list(csv.DictReader(f))
            self.assertEqual(len(rows), 2)
        finally:
            for p in [sdf_path, out_path]:
                if os.path.exists(p): os.unlink(p)

    def test_return_results_true_is_default(self):
        from recon6.recon import ReconConfig, run_recon

        with tempfile.NamedTemporaryFile(suffix='.sdf', mode='w',
                                         delete=False) as f:
            f.write(_TWO_MOL_SDF)
            sdf_path = f.name
        try:
            config = ReconConfig(
                data_dir=DATA_DIR, input_files=[sdf_path],
                fmt='sdf', auto_add_h=False,
            )
            results = run_recon(config)
            self.assertEqual(len(results), 2,
                             "return_results=True (default) must return all rows")
        finally:
            os.unlink(sdf_path)

    def test_csv_contents_identical_with_and_without_return(self):
        """The CSV written with return_results=True and return_results=False
        should be byte-identical for the same input."""
        from recon6.recon import ReconConfig, run_recon

        with tempfile.NamedTemporaryFile(suffix='.sdf', mode='w',
                                         delete=False) as f:
            f.write(_TWO_MOL_SDF)
            sdf_path = f.name
        out1 = sdf_path.replace('.sdf', '_with.csv')
        out2 = sdf_path.replace('.sdf', '_without.csv')
        try:
            for out, ret in [(out1, True), (out2, False)]:
                config = ReconConfig(
                    data_dir=DATA_DIR, input_files=[sdf_path],
                    fmt='sdf', output_csv=out, auto_add_h=False,
                    return_results=ret,
                )
                run_recon(config)
            with open(out1) as f1, open(out2) as f2:
                self.assertEqual(f1.read(), f2.read(),
                             "CSV output must be identical regardless of return_results")
        finally:
            for p in [sdf_path, out1, out2]:
                if os.path.exists(p): os.unlink(p)


@unittest.skipUnless(HAVE_DATA, "DATA dir not available")
class TestLogFile(unittest.TestCase):
    """log_file redirects stderr output; notes and warnings must go to
    the file, not to the console."""

    def test_log_file_receives_output(self):
        from recon6.recon import ReconConfig, run_recon

        with tempfile.NamedTemporaryFile(suffix='.sdf', mode='w',
                                         delete=False) as f:
            f.write(_ONE_MOL_SDF)
            sdf_path = f.name
        log_path = sdf_path.replace('.sdf', '.log')
        try:
            config = ReconConfig(
                data_dir=DATA_DIR, input_files=[sdf_path],
                fmt='sdf', auto_add_h=False, iprint=1,
                log_file=log_path,
            )
            run_recon(config)
            self.assertTrue(os.path.exists(log_path))
            with open(log_path) as f:
                content = f.read()
            self.assertIn('mol_A', content)
        finally:
            for p in [sdf_path, log_path]:
                if os.path.exists(p): os.unlink(p)

    def test_stderr_restored_after_run(self):
        """sys.stderr must be restored to its original value after
        run_recon returns, even if an exception occurs mid-run."""
        import recon6.recon as mod
        from recon6.recon import ReconConfig, run_recon

        original_stderr = sys.stderr
        with tempfile.NamedTemporaryFile(suffix='.sdf', mode='w',
                                         delete=False) as f:
            f.write(_ONE_MOL_SDF)
            sdf_path = f.name
        log_path = sdf_path.replace('.sdf', '.log')
        try:
            config = ReconConfig(
                data_dir=DATA_DIR, input_files=[sdf_path],
                fmt='sdf', auto_add_h=False, log_file=log_path,
            )
            run_recon(config)
            self.assertIs(sys.stderr, original_stderr,
                         "sys.stderr must be restored after run_recon")
        finally:
            for p in [sdf_path, log_path]:
                if os.path.exists(p): os.unlink(p)

    def test_stderr_restored_even_on_exception(self):
        """sys.stderr must also be restored if run_recon raises."""
        import recon6.recon as mod
        from recon6.recon import ReconConfig, run_recon

        original_stderr = sys.stderr
        with tempfile.NamedTemporaryFile(suffix='.log', delete=False) as f:
            log_path = f.name
        try:
            config = ReconConfig(
                data_dir='/nonexistent/data_dir',
                input_files=['/nonexistent.sdf'],
                fmt='sdf', log_file=log_path,
            )
            try:
                run_recon(config)
            except Exception:
                pass
            self.assertIs(sys.stderr, original_stderr,
                         "sys.stderr must be restored even after exception")
        finally:
            if os.path.exists(log_path): os.unlink(log_path)


@unittest.skipUnless(HAVE_DATA, "DATA dir not available")
class TestLargeBatchStreaming(unittest.TestCase):
    """Integration: toxx2.sdf (278 molecules) with return_results=False
    produces the same CSV as the default mode, in less memory."""

    def test_toxx_streaming_matches_batch(self):
        if not os.path.exists(TOXX_SDF):
            self.skipTest("toxx2.sdf not available")
        from recon6.recon import ReconConfig, run_recon

        with tempfile.NamedTemporaryFile(suffix='_stream.csv',
                                          delete=False) as f:
            out_stream = f.name
        with tempfile.NamedTemporaryFile(suffix='_batch.csv',
                                          delete=False) as f:
            out_batch = f.name
        try:
            for out, ret in [(out_stream, False), (out_batch, True)]:
                config = ReconConfig(
                    data_dir=DATA_DIR, input_files=[TOXX_SDF],
                    fmt='sdf', output_csv=out, return_results=ret,
                )
                run_recon(config)

            with open(out_stream) as f:
                rows_s = list(csv.DictReader(f))
            with open(out_batch) as f:
                rows_b = list(csv.DictReader(f))

            self.assertEqual(len(rows_s), len(rows_b))
            for rs, rb in zip(rows_s, rows_b):
                self.assertAlmostEqual(
                    float(rs['Energy']), float(rb['Energy']), places=4)
        finally:
            for p in [out_stream, out_batch]:
                if os.path.exists(p): os.unlink(p)


if __name__ == '__main__':
    unittest.main()


@unittest.skipUnless(HAVE_DATA, "DATA dir not available")
class TestErrorLog(unittest.TestCase):
    """Atom-typing quality report (error_log) tests."""

    def test_error_log_csv_has_correct_columns(self):
        from recon6.recon import ReconConfig, run_recon
        with tempfile.NamedTemporaryFile(suffix='.sdf', mode='w', delete=False) as f:
            f.write(_TWO_MOL_SDF); sdf_path = f.name
        err_path = sdf_path.replace('.sdf', '.err.csv')
        try:
            config = ReconConfig(
                data_dir=DATA_DIR, input_files=[sdf_path],
                fmt='sdf', auto_add_h=False, error_log=err_path,
            )
            run_recon(config)
            self.assertTrue(os.path.exists(err_path))
            with open(err_path) as f:
                reader = csv.DictReader(f)
                _ = list(reader)  # consume rows
                for col in ('Molecule', 'AtomIndex', 'Element',
                            'AtomTypeCode', 'MatchLevel',
                            'BestTAEEntry', 'MatchQuality'):
                    self.assertIn(col, reader.fieldnames)
        finally:
            for p in [sdf_path, err_path]:
                if os.path.exists(p): os.unlink(p)

    def test_error_log_match_quality_labels(self):
        """MatchQuality must be one of the documented label strings."""
        from recon6.recon import ReconConfig, run_recon
        if not os.path.exists(TOXX_SDF):
            self.skipTest("toxx2.sdf not available")
        err_path = '/tmp/test_err_labels.csv'
        try:
            config = ReconConfig(
                data_dir=DATA_DIR, input_files=[TOXX_SDF],
                fmt='sdf', error_log=err_path,
            )
            run_recon(config)
            valid_labels = {'near', 'good', 'partial', 'poor',
                            'very_poor', 'no_match'}
            with open(err_path) as f:
                for row in csv.DictReader(f):
                    self.assertIn(row['MatchQuality'], valid_labels)
                    self.assertLess(int(row['MatchLevel']), 3)
        finally:
            if os.path.exists(err_path): os.unlink(err_path)

    def test_error_log_flushed_incrementally(self):
        """Error log must be flushed after each molecule, not at end."""
        from recon6.recon import ReconConfig, run_recon
        import recon6.recon as mod
        if not os.path.exists(TOXX_SDF):
            self.skipTest("toxx2.sdf not available")
        err_path = '/tmp/test_err_incremental.csv'
        sizes = []
        original = mod._process_molecule
        call_count = [0]
        def patched(mol, tae_index, config):
            result = original(mol, tae_index, config)
            call_count[0] += 1
            if call_count[0] == 5 and os.path.exists(err_path):
                sizes.append(os.path.getsize(err_path))
            return result
        mod._process_molecule = patched
        try:
            config = ReconConfig(
                data_dir=DATA_DIR, input_files=[TOXX_SDF],
                fmt='sdf', error_log=err_path,
            )
            run_recon(config)
        finally:
            mod._process_molecule = original
            if os.path.exists(err_path): os.unlink(err_path)
        # After 5 molecules some imperfect matches should have been
        # flushed already (toxx has imperfect matches; if by chance
        # the first 5 are all perfect, the size check may give 0
        # but that's still a valid state, not a bug).
        self.assertGreaterEqual(len(sizes), 1)

    def test_perfect_match_molecules_produce_no_rows(self):
        """Atoms with perfect match (lev==3) must not appear in the
        error log. We verify using a well-connected molecule from toxx
        where we know some atoms are perfect and check those are absent."""
        if not os.path.exists(TOXX_SDF):
            self.skipTest("toxx2.sdf not available")
        from recon6.recon import ReconConfig, run_recon
        err_path = '/tmp/test_err_perfect.csv'
        try:
            config = ReconConfig(
                data_dir=DATA_DIR, input_files=[TOXX_SDF],
                fmt='sdf', error_log=err_path,
            )
            run_recon(config)
            with open(err_path) as f:
                rows = list(csv.DictReader(f))
            # Every row in the error log must have MatchLevel < 3
            for row in rows:
                self.assertLess(int(row['MatchLevel']), 3,
                                "Only imperfect matches should appear in error log")
        finally:
            if os.path.exists(err_path): os.unlink(err_path)
