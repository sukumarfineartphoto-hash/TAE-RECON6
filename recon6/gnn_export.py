"""
GNN (Graph Neural Network) export for RECON6.

Exports per-atom TAE descriptor data and bond connectivity in a format
ready for graph learning pipelines. The schema mirrors the QMF output
format used by the original Fortran RECON5:

  - Nodes: per-atom identity (element, coords, atomic number) plus the
    full 167-field TAE descriptor vector, using the canonical column
    names from the Fortran QMF format.
  - Edges: source atom, destination atom, bond order (1-indexed atom
    indices converted to 0-indexed for GNN convention). Bonds are
    represented as two directed edges per bond (i->j and j->i), which
    is the standard convention for message-passing GNNs.
  - Graph: whole-molecule descriptors (the same 249-column set already
    written to recon.csv).

Output formats
--------------
JSONL (default, recommended for large batches):
    One JSON object per line, each line self-contained.  Appendable
    and streamable; one `output_gnn` path suffices for the whole run.

JSON (one file per molecule, for small batches or interactive use):
    Pass output_gnn as a directory; one <molecule_id>.json per molecule.

PyTorch Geometric helper
-------------------------
    from recon6.gnn_export import to_pyg_data
    data = to_pyg_data(mol_json)   # returns torch_geometric.data.Data
PyTorch Geometric must be installed separately; it is not a dependency
of RECON6 itself.
"""

import json
import os
import re

# Canonical atom-block column names from the Fortran QMF format,
# in the same order they appear in the !1 header line.  The first 5
# ('Atom','x','y','z','AtNum') are atom identity, not TAE descriptor
# features; the remaining 167 are the per-atom node feature vector.
QMF_ATOM_COLS = [
    'Atom', 'x', 'y', 'z', 'AtNum',
    'Energy', 'Population', 'Volume',
    'DipoleMoment', 'DxMoment', 'DyMoment', 'DzMoment',
    'QAA', 'QBB', 'QCC',
    'QxxMoment', 'QxyMoment', 'QxzMoment',
    'QyyMoment', 'QyzMoment', 'QzzMoment',
    'SurfArea',
    'SIDel(Rho)N', 'Del(Rho)NMin', 'Del(Rho)NMax', 'Del(Rho)NIA',
]
QMF_ATOM_COLS += ['Del(Rho)NA%d' % k for k in range(1, 11)]
QMF_ATOM_COLS += [
    'SIDel(K)N', 'Del(K)Min', 'Del(K)Max', 'Del(K)IA',
]
QMF_ATOM_COLS += ['Del(K)NA%d' % k for k in range(1, 11)]
QMF_ATOM_COLS += ['SIK', 'SIKMin', 'SIKMax', 'SIKIA']
QMF_ATOM_COLS += ['SIKA%d' % k for k in range(1, 11)]
QMF_ATOM_COLS += [
    'SIDel(G)N', 'Del(G)NMin', 'Del(G)NMax', 'Del(G)NIA',
]
QMF_ATOM_COLS += ['Del(G)NA%d' % k for k in range(1, 11)]
QMF_ATOM_COLS += ['SIG', 'SIGMin', 'SIGMax', 'SIGIA']
QMF_ATOM_COLS += ['SIGA%d' % k for k in range(1, 11)]
QMF_ATOM_COLS += ['SIEP', 'SIEPMin', 'SIEPMax', 'SIEPIA']
QMF_ATOM_COLS += ['SIEPA%d' % k for k in range(1, 11)]
QMF_ATOM_COLS += [
    'piV', 'sigmaPV', 'sigmaNV', 'sumsigma', 'sigmanew',
]
QMF_ATOM_COLS += ['EP%d' % k for k in range(1, 11)]
QMF_ATOM_COLS += ['PIPMin', 'PIPMax', 'PIPAvg']
QMF_ATOM_COLS += ['PIP%d' % k for k in range(1, 21)]
QMF_ATOM_COLS += ['Fuk', 'FukMin', 'FukMax', 'FukAvg']
QMF_ATOM_COLS += ['Fuk%d' % k for k in range(1, 11)]
QMF_ATOM_COLS += ['LapL', 'LaplMin', 'LaplMax', 'LaplAvg']
QMF_ATOM_COLS += ['Lapl%d' % k for k in range(1, 11)]

assert len(QMF_ATOM_COLS) == 172, len(QMF_ATOM_COLS)

# The 5 identity fields; everything else is a node feature.
_IDENTITY_COLS = {'Atom', 'x', 'y', 'z', 'AtNum'}
NODE_FEATURE_COLS = [c for c in QMF_ATOM_COLS if c not in _IDENTITY_COLS]
assert len(NODE_FEATURE_COLS) == 167

# Mapping from our internal taedat field names (as stored in
# atom_records by descriptors.py) to the canonical QMF column names.
# taedat.py uses abbreviated lowercase keys; the QMF uses the
# full mixed-case names from the Fortran output.
_TAEDAT_TO_QMF = {
    # scalars
    'energy':   'Energy',
    'pop':      'Population',
    'vol':      'Volume',
    'dip':      'DipoleMoment',
    'dx':       'DxMoment',
    'dy':       'DyMoment',
    'dz':       'DzMoment',
    'qaa':      'QAA',
    'qbb':      'QBB',
    'qcc':      'QCC',
    'qxx':      'QxxMoment',
    'qxy':      'QxyMoment',
    'qxz':      'QxzMoment',
    'qyy':      'QyyMoment',
    'qyz':      'QyzMoment',
    'qzz':      'QzzMoment',
    'sa':       'SurfArea',
    # Del(Rho)N family
    'sidrn':    'SIDel(Rho)N',
    'drnmn':    'Del(Rho)NMin',
    'drnmx':    'Del(Rho)NMax',
    'drnia':    'Del(Rho)NIA',
    # Del(K)N family
    'sidkn':    'SIDel(K)N',
    'dknmn':    'Del(K)Min',
    'dknmx':    'Del(K)Max',
    'dknia':    'Del(K)IA',
    # SIK family
    'sik':      'SIK',
    'sikmn':    'SIKMin',
    'sikmx':    'SIKMax',
    'sikia':    'SIKIA',
    # Del(G)N family
    'sidgn':    'SIDel(G)N',
    'dgnmn':    'Del(G)NMin',
    'dgnmx':    'Del(G)NMax',
    'dgnia':    'Del(G)NIA',
    # SIG family
    'sig':      'SIG',
    'sigmn':    'SIGMin',
    'sigmx':    'SIGMax',
    'sigia':    'SIGIA',
    # SIEP family
    'siep':     'SIEP',
    'siepmn':   'SIEPMin',
    'siepmx':   'SIEPMax',
    'siepia':   'SIEPIA',
    # surface decomposition
    'piv':      'piV',
    'spv':      'sigmaPV',
    'snv':      'sigmaNV',
    'sums':     'sumsigma',
    'sign':     'sigmanew',
    # Fuk
    'fuk':      'Fuk',
    'fukmn':    'FukMin',
    'fukmx':    'FukMax',
    'fukav':    'FukAvg',
    # LapL
    'lapl':     'LapL',
    'lplmn':    'LaplMin',
    'lplmx':    'LaplMax',
    'lplav':    'LaplAvg',
}
# Histogram bin families: taedat key prefix -> QMF name prefix, count
_BIN_FAMILIES = [
    ('drna',  'Del(Rho)NA', 10),
    ('dkna',  'Del(K)NA',   10),
    ('sika',  'SIKA',       10),
    ('dgna',  'Del(G)NA',   10),
    ('siga',  'SIGA',       10),
    ('siepa', 'SIEPA',      10),
    ('ep',    'EP',         10),
    ('pip',   'PIP',        20),
    ('fuka',  'Fuk',        10),
    ('lpla',  'Lapl',       10),
]
for _src, _dst, _n in _BIN_FAMILIES:
    for _k in range(1, _n + 1):
        _TAEDAT_TO_QMF['%s%d' % (_src, _k)] = '%s%d' % (_dst, _k)


def _atom_record_to_qmf_dict(atom_rec, atom_sym, coords, atnum):
    """Convert one atom_records entry (taedat field names) plus its
    identity fields into a dict keyed by canonical QMF column names.
    Missing taedat fields are left as None."""
    d = {
        'Atom':  atom_sym.strip(),
        'x':     coords[0],
        'y':     coords[1],
        'z':     coords[2],
        'AtNum': atnum,
    }
    for tae_key, qmf_key in _TAEDAT_TO_QMF.items():
        d[qmf_key] = atom_rec.get(tae_key)
    return d


def mol_to_gnn_dict(mol, desc, molecule_id=None):
    """
    Build a GNN-ready dict for one molecule.

    Parameters
    ----------
    mol : dict returned by any RECON6 reader
    desc : dict returned by compute_descriptors (or _process_molecule)
    molecule_id : str, optional; uses mol['_name'] if not given

    Returns
    -------
    dict with keys:
        molecule_id : str
        atoms : list of dicts, one per atom, keyed by QMF column names.
                Includes identity fields (Atom, x, y, z, AtNum) and
                all 167 node feature fields.
        bonds : list of dicts with keys source (0-indexed), dest
                (0-indexed), bond_order.  Every bond appears as two
                directed entries (i->j and j->i).
        molecule_descriptors : dict of whole-molecule TAE descriptors
                (the 249-column set from recon.csv), for use as graph-
                level labels or features.
    """
    mol_id = molecule_id or mol.get('_name', '?')
    atom_records = desc.get('atom_records', [])
    natom = mol['natom']
    nuc = mol['nuc']

    # Build per-atom node dicts.
    # atom_records is indexed 0..(natom-1) for atoms 1..natom
    atoms = []
    for i in range(1, natom + 1):
        rec = atom_records[i - 1] if (i - 1) < len(atom_records) else {}
        sym = (mol['atom'][i] or '?').strip()
        coords = mol['coords'][i]
        d = _atom_record_to_qmf_dict(rec, sym, coords, int(nuc[i]))
        atoms.append(d)

    # Build per-bond edge list: two directed edges per bond, 0-indexed.
    bonds = []
    seen = set()
    for i in range(1, natom + 1):
        for j_pos in range(1, mol['ival'][i] + 1):
            j = mol['idcon'][i][j_pos]
            if j == 0:
                continue
            pair = (min(i, j), max(i, j))
            if pair in seen:
                continue
            seen.add(pair)
            bo = mol['nbo'][i][j] if mol['nbo'][i][j] else 1
            bonds.append({'source': i - 1, 'dest': j - 1,
                          'bond_order': int(bo)})
            bonds.append({'source': j - 1, 'dest': i - 1,
                          'bond_order': int(bo)})

    # Whole-molecule descriptors (strip internal _ keys).
    mol_desc = {k: v for k, v in desc.items()
                if not k.startswith('_') and k != 'atom_records'}

    return {
        'molecule_id':           mol_id,
        'atoms':                 atoms,
        'bonds':                 bonds,
        'molecule_descriptors':  mol_desc,
    }


# ---------------------------------------------------------------------------
# JSONL writer (streaming, one line per molecule)
# ---------------------------------------------------------------------------

class GnnJsonlWriter:
    """Appends one JSON-Lines record per molecule to a .jsonl file.

    Each line is a complete, self-contained JSON object produced by
    mol_to_gnn_dict.  The file can be read lazily line-by-line without
    loading the whole batch into memory, making it suitable for very
    large datasets.

    Usage::
        with GnnJsonlWriter('output.jsonl') as w:
            w.write(mol, desc)
    """

    def __init__(self, path):
        self.path = path
        self._fh = open(path, 'a')   # append mode: safe to resume
        self._n = 0

    def write(self, mol, desc, molecule_id=None):
        d = mol_to_gnn_dict(mol, desc, molecule_id=molecule_id)
        self._fh.write(json.dumps(d, allow_nan=False,
                                   separators=(',', ':')) + '\n')
        self._fh.flush()
        self._n += 1

    def close(self):
        if not self._fh.closed:
            self._fh.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


# ---------------------------------------------------------------------------
# Per-molecule JSON writer (one file per molecule)
# ---------------------------------------------------------------------------

def write_mol_json(mol, desc, directory, molecule_id=None,
                   sanitise_filename=True):
    """Write one molecule as a JSON file inside `directory`.

    The filename is <molecule_id>.json, with characters that are
    unsafe in filenames replaced by underscores when
    sanitise_filename=True.
    """
    os.makedirs(directory, exist_ok=True)
    mol_id = molecule_id or mol.get('_name', 'mol')
    fname = re.sub(r'[^\w\-.]', '_', mol_id) + '.json' \
        if sanitise_filename else mol_id + '.json'
    path = os.path.join(directory, fname)
    d = mol_to_gnn_dict(mol, desc, molecule_id=mol_id)
    with open(path, 'w') as fh:
        json.dump(d, fh, allow_nan=False, indent=2)
    return path


# ---------------------------------------------------------------------------
# Optional PyTorch Geometric helper (only importable when PyG is installed)
# ---------------------------------------------------------------------------

def to_pyg_data(gnn_dict):
    """
    Convert a mol_to_gnn_dict result into a PyTorch Geometric Data object.

    Requires torch and torch_geometric to be installed.  RECON6 itself
    does not depend on either.

    Node feature matrix x has shape [num_atoms, 167] using the 167
    NODE_FEATURE_COLS in canonical QMF order.  Any missing values
    (atoms whose TAE atom type was not found in the database) are
    represented as 0.0.

    Edge connectivity edge_index has shape [2, 2*num_bonds] (two
    directed edges per bond).  edge_attr has shape [2*num_bonds, 1]
    holding bond order as a float.

    pos has shape [num_atoms, 3] with x/y/z coordinates.

    y is None unless you add it yourself after the call.
    """
    try:
        import torch
        from torch_geometric.data import Data
    except ImportError as e:
        raise ImportError(
            "to_pyg_data requires PyTorch and PyTorch Geometric. "
            "Install them separately: pip install torch torch_geometric"
        ) from e

    atoms = gnn_dict['atoms']
    bonds = gnn_dict['bonds']

    # Node feature matrix [N, 167]
    x_rows = []
    for atom in atoms:
        row = [float(atom.get(col) or 0.0) for col in NODE_FEATURE_COLS]
        x_rows.append(row)
    x = torch.tensor(x_rows, dtype=torch.float)

    # 3-D positions [N, 3]
    pos = torch.tensor(
        [[float(a['x']), float(a['y']), float(a['z'])] for a in atoms],
        dtype=torch.float)

    # Edge index [2, E] and edge attr [E, 1]
    if bonds:
        edge_index = torch.tensor(
            [[b['source'] for b in bonds],
             [b['dest']   for b in bonds]],
            dtype=torch.long)
        edge_attr = torch.tensor(
            [[float(b['bond_order'])] for b in bonds],
            dtype=torch.float)
    else:
        edge_index = torch.zeros((2, 0), dtype=torch.long)
        edge_attr  = torch.zeros((0, 1), dtype=torch.float)

    return Data(x=x, edge_index=edge_index, edge_attr=edge_attr,
                pos=pos, mol_id=gnn_dict['molecule_id'])
