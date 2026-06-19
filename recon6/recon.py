"""
Top-level RECON6 orchestrator - translates the main program loop from
recon5-5.f. Ties together readers, atom typer, TAE lookup, and
descriptor accumulation.
"""
import os
import csv
import sys
import contextlib
from dataclasses import dataclass, field
from typing import Optional

from .bonds import load_bond_table
from .gettae import TaeIndex, build_modtype, gettae
from .ringid import ringid
from .descriptors import compute_descriptors
from .hydrogenate import needs_hydrogens, add_missing_hydrogens
from .readers.sdf import read_sdf
from .readers.mol2 import read_mol2
from .readers.pdb import read_pdb
from .readers.smiles import parse_smiles
from .readers.gaussian import read_gaussian_com
from .readers.orca import read_orca_xyz


@dataclass
class ReconConfig:
    data_dir: str          # path to DATA/ directory (TAE .dat files)
    input_files: list      # list of file paths (or SMILES strings if fmt='smiles')
    bond_file: Optional[str] = None  # path to bond-length table; defaults to '<data_dir>/bond'
    fmt: str = 'auto'      # 'sdf', 'mol2', 'pdb', 'gaussian', 'orca', 'smiles', or 'auto'
    output_csv: Optional[str] = None   # write recon.csv-style output here
    output_gnn: Optional[str] = None   # write GNN export here; path ending in
                                        # '.jsonl' streams one record per line
                                        # (recommended for large batches); any
                                        # other path is treated as a directory
                                        # and one <molecule_id>.json is written
                                        # per molecule.
    iprint: int = 0        # 0 = quiet, 1 = verbose (per-molecule progress)
    iovr: int = -1         # -1 = auto (CONECT or distance fallback); >0 = force distance
    auto_add_h: bool = True  # add missing H to SDF/MOL2/PDB; never applied to Gaussian/ORCA
    include_data_fields: bool = True  # carry SDF "> <TAG>" fields into output_csv
    return_results: bool = True  # accumulate and return results list; set False for
                                  # large batches to keep memory bounded (returns [])
    log_file: Optional[str] = None  # redirect all notes/warnings/progress to this file
                                     # instead of stderr; useful for unattended large runs

    def resolved_bond_file(self):
        """Return the bond-length table path: the explicit bond_file if
        given, otherwise '<data_dir>/bond'."""
        if self.bond_file:
            return self.bond_file
        return os.path.join(self.data_dir, 'bond')


def _detect_fmt(path):
    ext = os.path.splitext(path)[1].lower()
    return {'.sdf': 'sdf', '.mol2': 'mol2', '.pdb': 'pdb',
            '.smi': 'smiles', '.smiles': 'smiles',
            '.com': 'gaussian', '.gjf': 'gaussian',
            '.inp': 'orca', '.orca': 'orca'}.get(ext, 'sdf')


# ---------------------------------------------------------------------------
# Streaming CSV writer
# ---------------------------------------------------------------------------

# Descriptor columns in canonical order - fully fixed and known in advance,
# so the CSV header can be written immediately without scanning all results.
_BASE_KEYS = ['Molecule', 'Energy', 'Population', 'VOLTAE', 'SurfArea',
              'SIDel_RhoN', 'Del_RhoNMin', 'Del_RhoNMax', 'Del_RhoNIA']
_BASE_KEYS += ['Del_RhoNA%d' % k for k in range(1, 11)]
_BASE_KEYS += ['SIDel_KN', 'Del_KMin', 'Del_KMax', 'Del_KIA']
_BASE_KEYS += ['Del_KNA%d' % k for k in range(1, 11)]
_BASE_KEYS += ['SIK', 'SIKMin', 'SIKMax', 'SIKIA']
_BASE_KEYS += ['SIKA%d' % k for k in range(1, 11)]
_BASE_KEYS += ['SIDel_GN', 'Del_GNMin', 'Del_GNMax', 'Del_GNIA']
_BASE_KEYS += ['Del_GNA%d' % k for k in range(1, 11)]
_BASE_KEYS += ['SIG', 'SIGMin', 'SIGMax', 'SIGIA']
_BASE_KEYS += ['SIGA%d' % k for k in range(1, 11)]
_BASE_KEYS += ['SIEP', 'SIEPMin', 'SIEPMax', 'SIEPIA']
_BASE_KEYS += ['SIEPA%d' % k for k in range(1, 11)]
_BASE_KEYS += ['EP%d' % k for k in range(1, 11)]
_BASE_KEYS += ['PIPMin', 'PIPMax', 'PIPAvg']
_BASE_KEYS += ['PIP%d' % k for k in range(1, 21)]
_BASE_KEYS += ['Fuk', 'FukMin', 'FukMax', 'FukAvg']
_BASE_KEYS += ['Fuk%d' % k for k in range(1, 11)]
_BASE_KEYS += ['Lapl', 'LaplMin', 'LaplMax', 'LaplAvg']
_BASE_KEYS += ['Lapl%d' % k for k in range(1, 11)]
_BASE_KEYS += ['FDRNA%d' % k for k in range(1, 11)]
_BASE_KEYS += ['FDKNA%d' % k for k in range(1, 11)]
_BASE_KEYS += ['FSIKA%d' % k for k in range(1, 11)]
_BASE_KEYS += ['FDGNA%d' % k for k in range(1, 11)]
_BASE_KEYS += ['FSIGA%d' % k for k in range(1, 11)]
_BASE_KEYS += ['FEP%d' % k for k in range(1, 11)]
_BASE_KEYS += ['FPIP%d' % k for k in range(1, 21)]
_BASE_KEYS += ['FFuk%d' % k for k in range(1, 11)]
_BASE_KEYS += ['FLapl%d' % k for k in range(1, 11)]
_BASE_KEYS += ['chi']

_INTERNAL_KEYS = {'atom_records'}


class _CsvStreamWriter:
    """Writes result dicts to a CSV file one row at a time.

    The base descriptor columns are known in advance, so the header is
    written immediately when the writer is opened. Any additional SDF
    data-field columns (e.g. activity/response tags) are discovered on
    the fly as rows arrive. If a new tag appears after the first row has
    already been written, the file is rewritten in-place with the
    expanded header and all previously buffered rows intact - this is
    rare in practice (one SDF file almost always has a consistent tag
    set across all molecules) but handled correctly.
    """

    def __init__(self, path, include_data_fields=True):
        self.path = path
        self.include_data_fields = include_data_fields
        self._fieldnames = list(_BASE_KEYS)
        self._seen = set(_BASE_KEYS)
        self._fh = open(path, 'w', newline='')
        self._writer = csv.DictWriter(
            self._fh, fieldnames=self._fieldnames, extrasaction='ignore')
        self._writer.writeheader()
        self._fh.flush()
        self._n = 0
        # buffer of dicts written so far, kept only when data_fields
        # are enabled and we need to be able to rewrite on new-tag
        # discovery; cleared once it's clear no rewrite is needed.
        self._buffer = [] if include_data_fields else None

    def write(self, row):
        if self.include_data_fields:
            new_tags = [k for k in row
                        if not k.startswith('_')
                        and k != 'Molecule'
                        and k not in _INTERNAL_KEYS
                        and k not in self._seen]
            if new_tags:
                self._fieldnames.extend(new_tags)
                for t in new_tags:
                    self._seen.add(t)
                # New column(s) appeared - rewrite file from scratch.
                self._fh.close()
                rows_so_far = list(self._buffer) if self._buffer else []
                self._fh = open(self.path, 'w', newline='')
                self._writer = csv.DictWriter(
                    self._fh, fieldnames=self._fieldnames,
                    extrasaction='ignore')
                self._writer.writeheader()
                for old_row in rows_so_far:
                    self._writer.writerow(
                        {k: old_row.get(k, '') for k in self._fieldnames})

        out = {k: row.get(k, '') for k in self._fieldnames}
        self._writer.writerow(out)
        self._fh.flush()
        self._n += 1
        if self._buffer is not None:
            self._buffer.append(dict(row))

    def close(self):
        if not self._fh.closed:
            self._fh.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


# ---------------------------------------------------------------------------
# Molecule iteration helpers
# ---------------------------------------------------------------------------

def _iter_molecules(src, fmt, config=None):
    """Yield (mol_or_none, position) pairs one at a time from src."""
    auto_add_h = config.auto_add_h if config is not None else True

    if fmt == 'smiles':
        if os.path.isfile(src):
            with open(src) as fh:
                smiles_lines = [l.strip() for l in fh if l.strip()]
        else:
            smiles_lines = [src]
        for idx, smi in enumerate(smiles_lines):
            try:
                mol = parse_smiles(smi)
                mol['_name'] = str(idx + 1)
                mol['_src'] = src
                yield mol, idx
            except Exception as e:
                print("SMILES parse error: %s - %s" % (smi[:40], e),
                      file=sys.stderr)
                yield None, idx
    elif fmt == 'sdf':
        from .readers.sdf import (_read_molecule_block, _parse_molecule_block,
                                   read_sdf_data_fields_only)
        with open(src) as fh:
            idx = 0
            while True:
                try:
                    block_lines = _read_molecule_block(fh)
                except StopIteration:
                    break
                try:
                    mol = _parse_molecule_block(block_lines)
                except Exception as e:
                    print("SDF read error mol %d: %s" % (idx + 1, e),
                          file=sys.stderr)
                    name, data_fields = read_sdf_data_fields_only(block_lines)
                    yield {'_parse_failed': True, '_name': name or str(idx + 1),
                           '_src': src, 'data_fields': data_fields}, idx
                    idx += 1
                    continue
                mol['_name'] = mol.get('mol_name') or str(idx + 1)
                mol['_src'] = src
                if auto_add_h and needs_hydrogens(mol):
                    n_before = mol['natom']
                    mol = add_missing_hydrogens(mol)
                    print("Note [%s]: added %d hydrogen(s) (input had none; "
                          "%d heavy atoms) using standard valence rules" % (
                              mol['_name'], mol['hydrogens_added'], n_before),
                          file=sys.stderr)
                yield mol, idx
                idx += 1
    elif fmt == 'mol2':
        with open(src) as fh:
            idx = 0
            while True:
                try:
                    mol = read_mol2(fh)
                except StopIteration:
                    break
                except Exception as e:
                    print("MOL2 read error mol %d: %s" % (idx + 1, e),
                          file=sys.stderr)
                    yield None, idx
                    idx += 1
                    continue
                mol['_name'] = str(idx + 1)
                mol['_src'] = src
                if auto_add_h and needs_hydrogens(mol):
                    n_before = mol['natom']
                    mol = add_missing_hydrogens(mol)
                    print("Note [%s]: added %d hydrogen(s) (input had none; "
                          "%d heavy atoms) using standard valence rules" % (
                              mol['_name'], mol['hydrogens_added'], n_before),
                          file=sys.stderr)
                yield mol, idx
                idx += 1
    else:
        print("Unknown format: %s" % fmt, file=sys.stderr)


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

def run_recon(config: ReconConfig):
    """Run RECON6 on the given configuration.

    Returns a list of result dicts (one per successfully processed
    molecule) when config.return_results is True (the default).

    Set config.return_results = False for large batches where holding
    the full list in memory is impractical - run_recon will stream-
    write each result immediately to output_csv (if configured) and
    return an empty list, keeping peak memory bounded to a single
    molecule at a time.

    Set config.log_file to redirect all notes, warnings, and progress
    output to a file instead of stderr - useful for unattended runs.
    """
    bondl = load_bond_table(config.resolved_bond_file())
    tae_index = TaeIndex(config.data_dir)

    results = []
    n_processed = 0

    # Redirect stderr to log_file for the duration of the run if
    # requested; all print(..., file=sys.stderr) calls in this module
    # and in the readers/hydrogenate stack will follow automatically.
    log_fh = None
    original_stderr = sys.stderr
    if config.log_file:
        log_fh = open(config.log_file, 'a', buffering=1)  # line-buffered
        sys.stderr = log_fh

    try:
        csv_writer = None
        if config.output_csv:
            csv_writer = _CsvStreamWriter(
                config.output_csv,
                include_data_fields=config.include_data_fields)

        gnn_writer = None
        if config.output_gnn:
            from .gnn_export import GnnJsonlWriter, write_mol_json
            if config.output_gnn.endswith('.jsonl'):
                gnn_writer = GnnJsonlWriter(config.output_gnn)
            else:
                # directory mode — writer is a sentinel, we'll call
                # write_mol_json per molecule in the loop below
                gnn_writer = config.output_gnn

        try:
            for src in config.input_files:
                fmt = (config.fmt if config.fmt != 'auto'
                       else _detect_fmt(src))

                if fmt == 'pdb':
                    try:
                        mol = read_pdb(src, bondl=bondl, iovr=config.iovr,
                                       auto_add_h=config.auto_add_h)
                        mol['_name'] = os.path.splitext(
                            os.path.basename(src))[0]
                        mol['_src'] = src
                        if mol.get('_pdb_hydrogens_added'):
                            print("Note [%s]: added %d hydrogen(s) (input "
                                  "had none; %d heavy atoms) using standard "
                                  "valence rules" % (
                                      mol['_name'],
                                      mol['_pdb_hydrogens_added'],
                                      mol['_pdb_natom_before_h']),
                                  file=sys.stderr)
                        mol_iter = [(mol, 0)]
                    except Exception as e:
                        print("PDB read error: %s - %s" % (src, e),
                              file=sys.stderr)
                        mol_iter = [(None, 0)]
                elif fmt == 'gaussian':
                    try:
                        mol = read_gaussian_com(src, bondl=bondl)
                        mol['_name'] = os.path.splitext(
                            os.path.basename(src))[0]
                        mol['_src'] = src
                        mol_iter = [(mol, 0)]
                    except Exception as e:
                        print("Gaussian input read error: %s - %s" % (src, e),
                              file=sys.stderr)
                        mol_iter = [(None, 0)]
                elif fmt == 'orca':
                    try:
                        mol = read_orca_xyz(src, bondl=bondl)
                        mol['_name'] = os.path.splitext(
                            os.path.basename(src))[0]
                        mol['_src'] = src
                        mol_iter = [(mol, 0)]
                    except Exception as e:
                        print("ORCA input read error: %s - %s" % (src, e),
                              file=sys.stderr)
                        mol_iter = [(None, 0)]
                else:
                    mol_iter = _iter_molecules(src, fmt, config)

                for mol, position in mol_iter:
                    data_fields = (
                        (mol or {}).get('data_fields', {}) if mol else {})
                    name = (mol or {}).get('_name', str(position + 1))

                    if mol is None or mol.get('_parse_failed'):
                        print("Skipped (excluded from output): %s" % name,
                              file=sys.stderr)
                        continue

                    try:
                        result = _process_molecule(mol, tae_index, config)
                        result.update(data_fields)
                        n_processed += 1
                        if config.iprint > 0:
                            print("Processed: %s  chi=%.4f  Energy=%.6f" % (
                                result['Molecule'], result['chi'],
                                result['Energy']),
                                file=sys.stderr)
                        if csv_writer:
                            csv_writer.write(result)
                        if gnn_writer is not None:
                            if isinstance(gnn_writer, str):
                                # directory mode
                                from .gnn_export import write_mol_json
                                write_mol_json(mol, result, gnn_writer)
                            else:
                                gnn_writer.write(mol, result)
                        if config.return_results:
                            results.append(result)
                    except Exception as e:
                        print("Error processing %s mol %s: %s" % (
                            mol.get('_src', '?'), mol.get('_name', '?'), e),
                            file=sys.stderr)
                        print("Skipped (excluded from output): %s" % name,
                              file=sys.stderr)
                    del mol

        finally:
            if csv_writer:
                csv_writer.close()
                print("Wrote %d rows to %s" % (
                    n_processed, config.output_csv), file=sys.stderr)
            if gnn_writer is not None and not isinstance(gnn_writer, str):
                gnn_writer.close()
                print("Wrote %d GNN records to %s" % (
                    n_processed, config.output_gnn), file=sys.stderr)

    finally:
        sys.stderr = original_stderr
        if log_fh:
            log_fh.close()

    return results


def _process_molecule(mol, tae_index, config):
    natom = mol['natom']
    ival = mol['ival']
    idcon = mol['idcon']
    icon = mol['icon']
    nuc = mol['nuc']

    isize = ringid(natom, ival, idcon)
    modtype = build_modtype(natom, ival, idcon, icon, isize, nuc)
    warnings = []
    atomtype_list, lev = gettae(natom, modtype, tae_index, warnings=warnings)
    for w in warnings:
        print("Warning [%s]: %s" % (mol.get('_name', '?'), w),
              file=sys.stderr)

    lflag = any(lev[i] <= 2 for i in range(1, natom + 1))

    desc = compute_descriptors(mol, atomtype_list, tae_index, lflag=lflag)
    desc['Molecule'] = mol.get('_name', '?')
    desc['_natom'] = natom
    desc['_num_h'] = sum(1 for i in range(1, natom + 1) if nuc[i] == 1)
    desc['_hydrogens_added'] = mol.get('hydrogens_added', 0)
    return desc
