"""
Top-level RECON5 orchestrator - translates the main program loop from
recon5-5.f. Ties together readers, atom typer, TAE lookup, and
descriptor accumulation.
"""
import os
import csv
import sys
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
    iprint: int = 0        # 0 = quiet, 1 = verbose
    iovr: int = -1         # -1 = use file connectivity; >0 = distance criterion
    auto_add_h: bool = True  # add missing hydrogens to SDF/MOL2/PDB input (standard valence rules); never applied to Gaussian/ORCA, which always specify hydrogens explicitly
    include_data_fields: bool = True  # carry SDF "> <TAG>" data-block fields into output_csv

    def resolved_bond_file(self):
        """Return the bond-length table path: the explicit bond_file if
        given, otherwise '<data_dir>/bond' (the standard location the
        TAE DATA directory ships it in)."""
        if self.bond_file:
            return self.bond_file
        return os.path.join(self.data_dir, 'bond')


def _detect_fmt(path):
    ext = os.path.splitext(path)[1].lower()
    return {'.sdf': 'sdf', '.mol2': 'mol2', '.pdb': 'pdb',
            '.smi': 'smiles', '.smiles': 'smiles',
            '.com': 'gaussian', '.gjf': 'gaussian',
            '.inp': 'orca', '.orca': 'orca'}.get(ext, 'sdf')


def _iter_molecules(src, fmt, config=None):
    """Yield (mol_or_none, position) pairs one at a time from src,
    without building a full in-memory list (important for large batch
    files). `mol_or_none` is None when a molecule at that position
    could not be read/parsed at all (so its data_fields, if any could
    still be salvaged, are unavailable too) - the caller is responsible
    for still emitting a row to preserve positional alignment with any
    external response/property data.
    """
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
                print("SMILES parse error: %s - %s" % (smi[:40], e), file=sys.stderr)
                yield None, idx
    elif fmt == 'sdf':
        from .readers.sdf import _read_molecule_block, _parse_molecule_block, read_sdf_data_fields_only
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
                    print("SDF read error mol %d: %s" % (idx + 1, e), file=sys.stderr)
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
                    print("MOL2 read error mol %d: %s" % (idx + 1, e), file=sys.stderr)
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


def run_recon(config: ReconConfig):
    """Run RECON5 on the given configuration. Returns a list of result dicts."""
    bondl = load_bond_table(config.resolved_bond_file())
    tae_index = TaeIndex(config.data_dir)

    results = []

    for src in config.input_files:
        fmt = config.fmt if config.fmt != 'auto' else _detect_fmt(src)

        if fmt == 'pdb':
            # PDB connectivity needs the bond table + iovr setting, so
            # handle it directly rather than through _iter_molecules.
            try:
                mol = read_pdb(src, bondl=bondl, iovr=config.iovr,
                                auto_add_h=config.auto_add_h)
                mol['_name'] = os.path.splitext(os.path.basename(src))[0]
                mol['_src'] = src
                if mol.get('_pdb_hydrogens_added'):
                    print("Note [%s]: added %d hydrogen(s) (input had none; "
                          "%d heavy atoms) using standard valence rules" % (
                              mol['_name'], mol['_pdb_hydrogens_added'],
                              mol['_pdb_natom_before_h']),
                          file=sys.stderr)
                mol_iter = [(mol, 0)]
            except Exception as e:
                print("PDB read error: %s - %s" % (src, e), file=sys.stderr)
                mol_iter = [(None, 0)]
        elif fmt == 'gaussian':
            # Gaussian .com/.gjf connectivity is always distance-based
            # (no bond list in the format at all), so this needs the
            # bond table directly, same as PDB. Note: hydrogens are
            # never added for Gaussian input (see readers/gaussian.py) -
            # this format always specifies all atoms explicitly.
            try:
                mol = read_gaussian_com(src, bondl=bondl)
                mol['_name'] = os.path.splitext(os.path.basename(src))[0]
                mol['_src'] = src
                mol_iter = [(mol, 0)]
            except Exception as e:
                print("Gaussian input read error: %s - %s" % (src, e), file=sys.stderr)
                mol_iter = [(None, 0)]
        elif fmt == 'orca':
            # ORCA "* xyz" input likewise has no bond list, so
            # connectivity is always distance-based. Hydrogens are
            # never added here either, for the same reason as Gaussian.
            try:
                mol = read_orca_xyz(src, bondl=bondl)
                mol['_name'] = os.path.splitext(os.path.basename(src))[0]
                mol['_src'] = src
                mol_iter = [(mol, 0)]
            except Exception as e:
                print("ORCA input read error: %s - %s" % (src, e), file=sys.stderr)
                mol_iter = [(None, 0)]
        else:
            mol_iter = _iter_molecules(src, fmt, config)

        for mol, position in mol_iter:
            data_fields = (mol or {}).get('data_fields', {}) if mol else {}
            name = (mol or {}).get('_name', str(position + 1))

            if mol is None or mol.get('_parse_failed'):
                # Could not read/parse this molecule at all - drop it
                # from the output entirely (no blank row). Any
                # external response/property data tracked by the
                # caller for this input position should likewise be
                # dropped to stay aligned with the surviving rows.
                print("Skipped (excluded from output): %s" % name, file=sys.stderr)
                continue

            try:
                result = _process_molecule(mol, tae_index, config)
                result.update(data_fields)
                results.append(result)
                if config.iprint > 0:
                    print("Processed: %s  chi=%.4f  Energy=%.6f" % (
                        result['Molecule'], result['chi'], result['Energy']))
            except Exception as e:
                print("Error processing %s mol %s: %s" % (
                    mol.get('_src', '?'), mol.get('_name', '?'), e), file=sys.stderr)
                print("Skipped (excluded from output): %s" % name, file=sys.stderr)
            # Drop the molecule dict (and its nbo/idcon arrays) immediately;
            # not strictly necessary now that nbo is sparse, but keeps peak
            # memory low regardless of input size.
            del mol

    if config.output_csv and results:
        _write_csv(results, config.output_csv,
                    include_data_fields=config.include_data_fields)

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
        print("Warning [%s]: %s" % (mol.get('_name', '?'), w), file=sys.stderr)

    lflag = any(lev[i] <= 2 for i in range(1, natom + 1))

    desc = compute_descriptors(mol, atomtype_list, tae_index, lflag=lflag)
    desc['Molecule'] = mol.get('_name', '?')
    desc['_natom'] = natom
    desc['_num_h'] = sum(1 for i in range(1, natom + 1) if nuc[i] == 1)
    desc['_hydrogens_added'] = mol.get('hydrogens_added', 0)
    return desc


def _write_csv(results, path, include_data_fields=True):
    # Collect all keys in deterministic order
    base_keys = ['Molecule', 'Energy', 'Population', 'VOLTAE', 'SurfArea',
                 'SIDel_RhoN', 'Del_RhoNMin', 'Del_RhoNMax', 'Del_RhoNIA']
    base_keys += ['Del_RhoNA%d'%k for k in range(1,11)]
    base_keys += ['SIDel_KN','Del_KMin','Del_KMax','Del_KIA']
    base_keys += ['Del_KNA%d'%k for k in range(1,11)]
    base_keys += ['SIK','SIKMin','SIKMax','SIKIA']
    base_keys += ['SIKA%d'%k for k in range(1,11)]
    base_keys += ['SIDel_GN','Del_GNMin','Del_GNMax','Del_GNIA']
    base_keys += ['Del_GNA%d'%k for k in range(1,11)]
    base_keys += ['SIG','SIGMin','SIGMax','SIGIA']
    base_keys += ['SIGA%d'%k for k in range(1,11)]
    base_keys += ['SIEP','SIEPMin','SIEPMax','SIEPIA']
    base_keys += ['SIEPA%d'%k for k in range(1,11)]
    base_keys += ['EP%d'%k for k in range(1,11)]
    base_keys += ['PIPMin','PIPMax','PIPAvg']
    base_keys += ['PIP%d'%k for k in range(1,21)]
    base_keys += ['Fuk','FukMin','FukMax','FukAvg']
    base_keys += ['Fuk%d'%k for k in range(1,11)]
    base_keys += ['Lapl','LaplMin','LaplMax','LaplAvg']
    base_keys += ['Lapl%d'%k for k in range(1,11)]
    base_keys += ['FDRNA%d'%k for k in range(1,11)]
    base_keys += ['FDKNA%d'%k for k in range(1,11)]
    base_keys += ['FSIKA%d'%k for k in range(1,11)]
    base_keys += ['FDGNA%d'%k for k in range(1,11)]
    base_keys += ['FSIGA%d'%k for k in range(1,11)]
    base_keys += ['FEP%d'%k for k in range(1,11)]
    base_keys += ['FPIP%d'%k for k in range(1,21)]
    base_keys += ['FFuk%d'%k for k in range(1,11)]
    base_keys += ['FLapl%d'%k for k in range(1,11)]
    base_keys += ['chi']

    fieldnames = list(base_keys)
    if include_data_fields:
        # Discover any extra SDF "> <TAG>" fields across all rows, in
        # first-seen order, and append them as trailing columns.
        seen = set(base_keys)
        extra = []
        # Internal bookkeeping keys that are never SDF data fields and
        # should never appear as CSV columns.
        internal_keys = {'atom_records'}
        for row in results:
            for k in row:
                if k.startswith('_') or k == 'Molecule' or k in internal_keys:
                    continue
                if k not in seen:
                    seen.add(k)
                    extra.append(k)
        fieldnames += extra

    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        for row in results:
            writer.writerow({k: row.get(k, '') for k in fieldnames})
    print("Wrote %d rows to %s" % (len(results), path))
