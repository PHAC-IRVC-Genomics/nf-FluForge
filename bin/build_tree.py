#!/usr/bin/env python3
"""
build_tree.py

This script builds phylogenetic trees from FASTA files (CDS and peptide) using MAFFT and FastTree.
It replicates the logic previously implemented in buildtree.sh, but in Python.

Usage:
  python build_tree.py -s /path/to/fasta_dir -t 8

Options:
  -s, --input     Input directory containing FASTA files (e.g. those produced by build_pseudogenome.py).
  -t, --threads   Number of threads to use (default: 10).
"""

import subprocess
import argparse
from pathlib import Path
import os
import sys
import datetime

def run_command(command):
    """Run a shell command and exit on error."""
    print(f"[Running] {command}")
    result = subprocess.run(command, shell=True)
    if result.returncode != 0:
        print(f"[Error] Command failed: {command}")
        sys.exit(result.returncode)

def main():
    parser = argparse.ArgumentParser(description="Build phylogenetic trees using MAFFT and FastTree.")
    parser.add_argument("-s", "--input", required=True,
                        help="Input directory containing FASTA files (e.g. *.cds.fasta and *.pep.fasta).")
    parser.add_argument("-t", "--threads", default=10, type=int,
                        help="Number of threads (default: 10).")
    args = parser.parse_args()

    input_dir = Path(args.input).resolve()
    if not input_dir.is_dir():
        print(f"[Error] Input directory {input_dir} does not exist.")
        sys.exit(1)

    threads = args.threads

    # Create output directory for tree files (named 'tree')
    tree_dir = Path("tree")
    tree_dir.mkdir(parents=True, exist_ok=True)

    print(f"{datetime.datetime.now()} Begin tree building")

    # Process CDS FASTA files
    for fasta_file in sorted(input_dir.glob("*.cds.fasta")):
        print(f"Processing CDS file: {fasta_file}")
        base_name = fasta_file.stem  # e.g., sample.pseudogenome.cds
        output_file = tree_dir / f"{base_name}.tre"
        cmd = f"mafft --auto --thread {threads} {fasta_file} | fasttree -nt -gtr > {output_file}"
        run_command(cmd)

    # Process peptide FASTA files
    for fasta_file in sorted(input_dir.glob("*.pep.fasta")):
        print(f"Processing peptide file: {fasta_file}")
        base_name = fasta_file.stem  # e.g., sample.pseudogenome.pep
        output_file = tree_dir / f"{base_name}.tre"
        cmd = f"mafft --auto --thread {threads} {fasta_file} | fasttree > {output_file}"
        run_command(cmd)

    print(f"{datetime.datetime.now()} Finished tree building")

if __name__ == "__main__":
    main()
