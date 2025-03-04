#!/usr/bin/env python

import argparse
from pathlib import Path

def sanitize_fasta_file(input_file: Path, output_file: Path):
    """Sanitize a FASTA file by replacing all '|' characters in header lines with '_'."""
    with input_file.open("r") as fin, output_file.open("w") as fout:
        for line in fin:
            if line.startswith(">"):
                sanitized_line = line.replace("|", "_")
                fout.write(sanitized_line)
            else:
                fout.write(line)

def sanitize_directory(input_dir: Path, output_dir: Path):
    """Sanitize all FASTA files (ending with .fasta) in input_dir, writing sanitized copies to output_dir."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for file in input_dir.glob("*.fasta"):
        output_file = output_dir / file.name
        sanitize_fasta_file(file, output_file)

def main():
    parser = argparse.ArgumentParser(
        description="Sanitize FASTA headers by replacing '|' with '_' in all FASTA files in a directory"
    )
    parser.add_argument("-i", "--input-dir", required=True, help="Input directory containing FASTA files")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory for sanitized FASTA files")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    
    sanitize_directory(input_dir, output_dir)

if __name__ == "__main__":
    main()
