"""Tests for recon6/gnn_export.py: schema, JSONL streaming,
per-molecule JSON, and PyG helper (skipped when PyG not installed)."""
import unittest
import sys, os, json, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.dirname(__file__))

from test_config import DATA_DIR, TOXX_SDF, HAVE_DATA

from recon6.gnn_export import (
    QMF_ATOM_COLS, NODE_FEATURE_COLS,
    mol_to_gnn_dict, GnnJsonlWriter, write_mol_json,
)

# Minimal mol + desc dicts sufficient to exercise the export
# without needing the TAE DATA directory.
_MOCK_NUC  = [0, 6, 8, 1, 1]   # 1-based: C, O, H, H
_MOCK_ATOM = [None, 'C ', 'O ', 'H ', 'H ']
_MOCK_COORDS = [[0,0,0], [0.0,0.0,0.0], [1.2,0.0,0.0],
                [-0.5, 0.9,0.0], [-0.5,-0.9,0.0]]
_MOCK_IDCON = [[0]*5 for _ in range(5)]
_MOCK_IVAL  = [0, 3, 1, 1, 1]
# C bonded to O, H, H
_MOCK_IDCON[1][1]=2; _MOCK_IDCON[1][2]=3; _MOCK_IDCON[1][3]=4
_MOCK_IDCON[2][1]=1
_MOCK_IDCON[3][1]=1
_MOCK_IDCON[4][1]=1

from recon6.sparse import SparseMatrix
_MOCK_NBO = SparseMatrix()
_MOCK_NBO[1][2]=2; _MOCK_NBO[2][1]=2  # C=O double bond
_MOCK_NBO[1][3]=1; _MOCK_NBO[3][1]=1
_MOCK_NBO[1][4]=1; _MOCK_NBO[4][1]=1

def _make_mock_mol(name='formaldehyde'):
    return dict(natom=4, atom=_MOCK_ATOM, nuc=_MOCK_NUC,
                coords=_MOCK_COORDS, idcon=_MOCK_IDCON, ival=_MOCK_IVAL,
                nbo=_MOCK_NBO, _name=name)

def _make_mock_desc(natom=4):
    # Minimal desc dict: atom_records with a handful of taedat fields
    atom_records = []
    for _ in range(natom):
        atom_records.append({'energy': -100.0, 'pop': 8.0, 'vol': 50.0,
                              'sa': 20.0})
    return {'atom_records': atom_records, 'Energy': -400.0,
            'Population': 32.0, 'chi': 1.5}


class TestSchemaConstants(unittest.TestCase):
    def test_qmf_atom_cols_count(self):
        self.assertEqual(len(QMF_ATOM_COLS), 172)

    def test_node_feature_cols_count(self):
        self.assertEqual(len(NODE_FEATURE_COLS), 167)

    def test_identity_cols_excluded_from_features(self):
        for col in ('Atom', 'x', 'y', 'z', 'AtNum'):
            self.assertNotIn(col, NODE_FEATURE_COLS)

    def test_canonical_col_names_present(self):
        for col in ('Energy', 'SIDel(Rho)N', 'PIP1', 'LapL', 'piV'):
            self.assertIn(col, QMF_ATOM_COLS)


class TestMolToGnnDict(unittest.TestCase):
    def setUp(self):
        self.mol = _make_mock_mol()
        self.desc = _make_mock_desc()
        self.gnn = mol_to_gnn_dict(self.mol, self.desc)

    def test_molecule_id(self):
        self.assertEqual(self.gnn['molecule_id'], 'formaldehyde')

    def test_atoms_count(self):
        self.assertEqual(len(self.gnn['atoms']), 4)

    def test_atom_identity_fields(self):
        c = self.gnn['atoms'][0]
        self.assertEqual(c['Atom'], 'C')
        self.assertEqual(c['AtNum'], 6)
        self.assertAlmostEqual(c['x'], 0.0)

    def test_atom_feature_field_mapped(self):
        # 'energy' in taedat -> 'Energy' in QMF
        self.assertAlmostEqual(self.gnn['atoms'][0]['Energy'], -100.0)

    def test_bonds_are_directed_pairs(self):
        # 3 bonds in formaldehyde -> 6 directed edges
        self.assertEqual(len(self.gnn['bonds']), 6)

    def test_bond_order_present(self):
        src_dest = {(b['source'], b['dest']): b['bond_order']
                    for b in self.gnn['bonds']}
        # C=O is bond_order 2 (indices 0 and 1, 0-indexed)
        self.assertEqual(src_dest[(0, 1)], 2)
        self.assertEqual(src_dest[(1, 0)], 2)

    def test_bond_indices_zero_based(self):
        for b in self.gnn['bonds']:
            self.assertGreaterEqual(b['source'], 0)
            self.assertGreaterEqual(b['dest'],   0)

    def test_molecule_descriptors_present(self):
        md = self.gnn['molecule_descriptors']
        self.assertAlmostEqual(md['Energy'], -400.0)
        self.assertAlmostEqual(md['chi'], 1.5)

    def test_internal_keys_excluded(self):
        md = self.gnn['molecule_descriptors']
        self.assertNotIn('atom_records', md)
        for k in md:
            self.assertFalse(k.startswith('_'))

    def test_json_serialisable(self):
        s = json.dumps(self.gnn)
        reloaded = json.loads(s)
        self.assertEqual(reloaded['molecule_id'], 'formaldehyde')


class TestGnnJsonlWriter(unittest.TestCase):
    def test_writes_one_line_per_molecule(self):
        with tempfile.NamedTemporaryFile(suffix='.jsonl', mode='w',
                                         delete=False) as f:
            path = f.name
        try:
            mol = _make_mock_mol('mol1')
            desc = _make_mock_desc()
            mol2 = _make_mock_mol('mol2')
            with GnnJsonlWriter(path) as w:
                w.write(mol, desc)
                w.write(mol2, desc)
            lines = open(path).read().strip().split('\n')
            self.assertEqual(len(lines), 2)
            obj0 = json.loads(lines[0])
            obj1 = json.loads(lines[1])
            self.assertEqual(obj0['molecule_id'], 'mol1')
            self.assertEqual(obj1['molecule_id'], 'mol2')
        finally:
            os.unlink(path)

    def test_file_flushed_after_each_write(self):
        with tempfile.NamedTemporaryFile(suffix='.jsonl', mode='w',
                                         delete=False) as f:
            path = f.name
        try:
            mol = _make_mock_mol()
            desc = _make_mock_desc()
            w = GnnJsonlWriter(path)
            w.write(mol, desc)
            # File must have content before close() is called
            size_mid = os.path.getsize(path)
            w.close()
            self.assertGreater(size_mid, 0,
                               "JSONL file must be flushed after each write")
        finally:
            os.unlink(path)

    def test_append_mode_survives_resume(self):
        """Writing to an existing JSONL file appends rather than
        overwriting - safe to resume a partial run."""
        with tempfile.NamedTemporaryFile(suffix='.jsonl', mode='w',
                                         delete=False) as f:
            path = f.name
        try:
            mol = _make_mock_mol()
            desc = _make_mock_desc()
            with GnnJsonlWriter(path) as w:
                w.write(mol, desc)
            with GnnJsonlWriter(path) as w:   # second run
                w.write(_make_mock_mol('mol2'), desc)
            lines = open(path).read().strip().split('\n')
            self.assertEqual(len(lines), 2)
        finally:
            os.unlink(path)


class TestWriteMolJson(unittest.TestCase):
    def test_creates_file_in_directory(self):
        mol = _make_mock_mol('acetone')
        desc = _make_mock_desc()
        with tempfile.TemporaryDirectory() as d:
            path = write_mol_json(mol, desc, d)
            self.assertTrue(os.path.exists(path))
            obj = json.load(open(path))
            self.assertEqual(obj['molecule_id'], 'acetone')

    def test_unsafe_chars_sanitised_in_filename(self):
        mol = _make_mock_mol('mol/with:bad\\chars')
        desc = _make_mock_desc()
        with tempfile.TemporaryDirectory() as d:
            path = write_mol_json(mol, desc, d)
            fname = os.path.basename(path)
            self.assertNotIn('/', fname)
            self.assertNotIn(':', fname)
            self.assertNotIn('\\', fname)


@unittest.skipUnless(HAVE_DATA, "DATA dir not available")
class TestGnnExportIntegration(unittest.TestCase):
    """End-to-end: run_recon with output_gnn set and validate schema
    of a real processed molecule."""

    def test_jsonl_output_from_run_recon(self):
        if not os.path.exists(TOXX_SDF):
            self.skipTest("toxx2.sdf not available")
        from recon6.recon import ReconConfig, run_recon

        with tempfile.NamedTemporaryFile(suffix='.jsonl', delete=False) as f:
            gnn_path = f.name
        try:
            config = ReconConfig(
                data_dir=DATA_DIR,
                input_files=[TOXX_SDF],
                fmt='sdf',
                output_gnn=gnn_path,
                return_results=False,
            )
            run_recon(config)
            lines = open(gnn_path).read().strip().split('\n')
            self.assertEqual(len(lines), 278)
            obj = json.loads(lines[0])
            self.assertIn('molecule_id',          obj)
            self.assertIn('atoms',                obj)
            self.assertIn('bonds',                obj)
            self.assertIn('molecule_descriptors', obj)
            # Every atom must have the expected identity fields
            for atom in obj['atoms']:
                for field in ('Atom', 'x', 'y', 'z', 'AtNum'):
                    self.assertIn(field, atom)
            # Bonds must be 0-indexed and directed
            natom = len(obj['atoms'])
            for bond in obj['bonds']:
                self.assertIn('source',     bond)
                self.assertIn('dest',       bond)
                self.assertIn('bond_order', bond)
                self.assertGreaterEqual(bond['source'], 0)
                self.assertLess(bond['source'], natom)
        finally:
            os.unlink(gnn_path)


@unittest.skipUnless(HAVE_DATA, "DATA dir not available")
class TestPyGHelper(unittest.TestCase):
    def test_to_pyg_data_raises_import_error_without_torch(self):
        """When PyTorch is not installed, to_pyg_data must raise
        ImportError with a helpful message rather than crashing
        opaquely."""
        import sys
        # Mock torch as absent if it isn't installed
        torch_available = 'torch' in sys.modules
        if torch_available:
            self.skipTest("torch is installed - skipping absent-torch test")
        from recon6.gnn_export import to_pyg_data
        gnn = mol_to_gnn_dict(_make_mock_mol(), _make_mock_desc())
        with self.assertRaises(ImportError):
            to_pyg_data(gnn)


if __name__ == '__main__':
    unittest.main()
