#!/usr/bin/env python3
import subprocess
import sys
import argparse
from pathlib import Path

def run(cmd):
    """Run a shell command; exit if it fails."""
    print(f"[Running] {cmd}")
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        print(f"[Error] Command failed with return code {result.returncode}")
        sys.exit(result.returncode)

def main():
    parser = argparse.ArgumentParser(
        description="Run VADR on consensus FASTA files and produce VADR output files."
    )
    parser.add_argument(
        "--consensus_dir",
        required=True,
        help="Directory containing consensus FASTA files (one per sample)."
    )
    parser.add_argument(
        "--vadr_model_dir",
        required=True,
        help="Path to VADR model directory (containing .minfo, .stk, etc.)."
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Output directory where each sample will have its own subfolder (default: current directory)."
    )
    parser.add_argument(
        "--file_ext",
        default=".fasta",
        help="File extension for consensus files (default: .fasta)."
    )
    args = parser.parse_args()

    consensus_dir = Path(args.consensus_dir)
    if not consensus_dir.is_dir():
        print(f"[Error] consensus_dir '{consensus_dir}' is not a valid directory.")
        sys.exit(1)

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Gather all consensus files with the specified extension
    consensus_files = sorted(consensus_dir.glob(f"*{args.file_ext}"))
    if not consensus_files:
        print(f"[Warning] No files ending with '{args.file_ext}' found in '{consensus_dir}'. Exiting.")
        sys.exit(0)

    for fasta_path in consensus_files:
        sample_name = fasta_path.stem  # e.g. "SAMPLE_01" from "SAMPLE_01.fasta"
        sample_outdir = outdir / sample_name
        sample_outdir.mkdir(parents=True, exist_ok=True)

        print(f"\n[INFO] Processing sample '{sample_name}' with FASTA = {fasta_path}\n")

        # Run VADR:
        vadr_cmd = (
            f"v-annotate.pl "
            f"--mdir {args.vadr_model_dir} "
            f"--mkey flu "  # this key is hard-coded; modify if needed
            f"-f {fasta_path} "
            f"--noseqnamemax "
            f"{sample_outdir}"
        )
        run(vadr_cmd)

        # Check for expected output files
        pass_tbl = sample_outdir / f"{sample_name}.vadr.pass.tbl"
        pass_fa  = sample_outdir / f"{sample_name}.vadr.pass.fa"
        if not pass_tbl.is_file():
            print(f"[Warning] {pass_tbl} not found. VADR may have failed to annotate any passing sequences.")
        if not pass_fa.is_file():
            print(f"[Warning] {pass_fa} not found. VADR may have failed to annotate any passing sequences.")
        else:
            print(f"[Done] VADR completed for sample '{sample_name}'.")
    
    print("\nAll samples processed with VADR. Exiting.\n")

if __name__ == "__main__":
    main()