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
    parser.add_argument('inputs', nargs='+', help='Input files (SDF/MOL2/PDB/Gaussian COM/SMILES)')
    parser.add_argument('--data-dir', required=True,
                        help="Path to TAE DATA/ directory (the bond-length table is read "
                             "automatically from '<data-dir>/bond' unless --bond-file overrides it)")
    parser.add_argument('--bond-file', default=None,
                        help="Override path to the bond-length table file "
                             "(default: '<data-dir>/bond')")
    parser.add_argument('--fmt', default='auto',
                        choices=['auto','sdf','mol2','pdb','gaussian','orca','smiles'],
                        help='Input file format (default: auto-detect from extension)')
    parser.add_argument('--output', '-o', default='recon.csv',
                        help='Output CSV file (default: recon.csv)')
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--no-add-h', action='store_true',
                        help='Disable automatic hydrogen addition for H-less SDF/MOL2/PDB input '
                             '(Gaussian/ORCA input is never modified this way)')
    parser.add_argument('--no-data-fields', action='store_true',
                        help="Don't carry SDF '> <TAG>' data-block fields (e.g. activity/response data) into the output CSV")
    parser.add_argument('--iovr', type=int, default=-1,
                        help='PDB connectivity mode: -1 = use CONECT records if present, '
                             'falling back to distance automatically when none exist (default); '
                             '>0 = force distance-based connectivity')
    args = parser.parse_args()

    config = ReconConfig(
        data_dir=args.data_dir,
        bond_file=args.bond_file,
        input_files=args.inputs,
        fmt=args.fmt,
        output_csv=args.output,
        iprint=1 if args.verbose else 0,
        auto_add_h=not args.no_add_h,
        include_data_fields=not args.no_data_fields,
        iovr=args.iovr,
    )

    results = run_recon(config)
    print("Done. %d molecule(s) processed." % len(results))

if __name__ == '__main__':
    main()
