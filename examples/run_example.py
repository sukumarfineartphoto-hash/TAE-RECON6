#!/usr/bin/env python3
"""
Example: compute TAE descriptors for one or more molecules.

Usage:
    python examples/run_example.py --data-dir DATA mol.sdf
    python examples/run_example.py --data-dir DATA ligand.pdb --fmt pdb
    python examples/run_example.py --data-dir DATA opt.com --fmt gaussian
    python examples/run_example.py --data-dir DATA opt.inp --fmt orca
"""
import argparse
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from recon6.recon import ReconConfig, run_recon


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input", help="SDF/MOL2/PDB/Gaussian COM/ORCA file, or a SMILES string/file")
    p.add_argument("--data-dir", required=True)
    p.add_argument("--bond-file", default=None,
                    help="Override path to the bond-length table (default: '<data-dir>/bond')")
    p.add_argument("--fmt", default="auto")
    args = p.parse_args()

    config = ReconConfig(
        data_dir=args.data_dir,
        bond_file=args.bond_file,
        input_files=[args.input],
        fmt=args.fmt,
        iprint=1,
    )
    results = run_recon(config)

    for r in results:
        print("\nMolecule: %s" % r["Molecule"])
        print("  Energy:      %12.4f" % r["Energy"])
        print("  Population:  %12.4f" % r["Population"])
        print("  SurfArea:    %12.4f" % r["SurfArea"])
        print("  Volume:      %12.4f" % r["VOLTAE"])
        print("  chi (Randic):%12.4f" % r["chi"])


if __name__ == "__main__":
    main()
