"""CLI entry point: python -m recon6 [options] file1 file2 ..."""
import argparse
import sys
import os
from .recon import ReconConfig, run_recon

def main():
    parser = argparse.ArgumentParser(
        prog='recon6',
        description='RECON6 TAE descriptor calculator (Python port of RECON5)'
    )
    parser.add_argument('inputs', nargs='+',
                        help='Input files (SDF/MOL2/PDB/Gaussian COM/ORCA/SMILES)')
    parser.add_argument('--data-dir', required=True,
                        help="Path to TAE DATA/ directory (bond-length table is read "
                             "automatically from '<data-dir>/bond' unless --bond-file overrides it)")
    parser.add_argument('--bond-file', default=None,
                        help="Override path to the bond-length table file "
                             "(default: '<data-dir>/bond')")
    parser.add_argument('--fmt', default='auto',
                        choices=['auto', 'sdf', 'mol2', 'pdb',
                                 'gaussian', 'orca', 'smiles'],
                        help='Input file format (default: auto-detect from extension)')
    parser.add_argument('--output', '-o', default='recon.csv',
                        help='Output CSV file (default: recon.csv)')
    parser.add_argument('--gnn-output', default=None,
                        help="Write GNN-ready JSON export here. Use a '.jsonl' "
                             "extension for streaming JSONL (one record per line, "
                             "recommended for large batches); any other path is "
                             "treated as a directory and one <molecule_id>.json "
                             "is written per molecule.")
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Print per-molecule progress (Molecule, chi, Energy)')
    parser.add_argument('--log-file', default=None,
                        help='Redirect all notes, warnings and progress output to '
                             'this file instead of stderr (recommended for large batches)')
    parser.add_argument('--error-log', default=None,
                        help='Write atom-typing quality report to this CSV file. '
                             'One row per atom with an imperfect TAE match (MatchLevel < 3), '
                             'with columns: Molecule, AtomIndex, Element, AtomTypeCode, '
                             'MatchLevel, BestTAEEntry, MatchQuality. '
                             'Strongly recommended when building ML models to assess '
                             'descriptor reliability across your dataset.')
    parser.add_argument('--no-return', action='store_true',
                        help='Do not accumulate results in memory - stream-write CSV/GNN '
                             'output only and return nothing (recommended for large batches '
                             'to keep memory bounded)')
    parser.add_argument('--no-add-h', action='store_true',
                        help='Disable automatic hydrogen addition for H-less '
                             'SDF/MOL2/PDB input (Gaussian/ORCA are never modified)')
    parser.add_argument('--no-data-fields', action='store_true',
                        help="Don't carry SDF '> <TAG>' data-block fields into "
                             "the output CSV")
    parser.add_argument('--iovr', type=int, default=-1,
                        help='PDB connectivity: -1 = use CONECT if present, '
                             'falling back to distance automatically (default); '
                             '>0 = force distance-based connectivity')
    args = parser.parse_args()

    config = ReconConfig(
        data_dir=args.data_dir,
        bond_file=args.bond_file,
        input_files=args.inputs,
        fmt=args.fmt,
        output_csv=args.output,
        output_gnn=args.gnn_output,
        iprint=1 if args.verbose else 0,
        auto_add_h=not args.no_add_h,
        include_data_fields=not args.no_data_fields,
        iovr=args.iovr,
        return_results=not args.no_return,
        log_file=args.log_file,
        error_log=args.error_log,
    )

    results = run_recon(config)
    if args.no_return:
        print("Done.")
    else:
        print("Done. %d molecule(s) processed." % len(results))

if __name__ == '__main__':
    main()
