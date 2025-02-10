#!/usr/bin/env python3

import subprocess
import sys
import argparse
from pathlib import Path
import time

def run(cmd):
    """Utility to run a shell command and fail if it returns non-zero."""
    print(f"[Running] {cmd}")
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        print(f"[Error] Command failed with return code {result.returncode}")
        sys.exit(result.returncode)

def main():
    parser = argparse.ArgumentParser(description="Run VADR + Table2asn on all consensus FASTA files in a directory.")
    parser.add_argument(
        "--consensus_dir",
        required=True,
        help="Path to a directory containing consensus FASTA files (one per sample)."
    )
    parser.add_argument(
        "--vadr_model_dir",
        required=True,
        help="Path to VADR model directory (containing .minfo, .stk, etc.)."
    )
    parser.add_argument(
        "--subseqids_script",
        default="sub_seqids_for_table2asn.py",
        help="Path to sub_seqids_for_table2asn.py (PRE_TABLE2ASN)"
    )
    parser.add_argument(
        "--posttable2asn_script",
        default="post_table2asn.py",
        help="Path to post_table2asn.py (POST_TABLE2ASN)"
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Directory for all output. Each sample will have its own subfolder."
    )
    parser.add_argument(
        "--file_ext",
        default=".fasta",
        help="File extension to look for in --consensus_dir (default .fasta)."
    )
    args = parser.parse_args()

    consensus_dir = Path(args.consensus_dir)
    if not consensus_dir.is_dir():
        print(f"[Error] consensus_dir '{consensus_dir}' is not a directory or does not exist.")
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

        # 1) Make a subdirectory for this sample's results
        sample_outdir = outdir / sample_name
        sample_outdir.mkdir(parents=True, exist_ok=True)

        print(f"\n[INFO] Processing sample '{sample_name}' with FASTA = {fasta_path}\n")

        # ---------------------------------------------------------------------
        # 2) Run VADR: creates {sample_name}.vadr.pass.tbl/.fa if successful
        # ---------------------------------------------------------------------
        # We'll store them in sample_outdir
        vadr_cmd = (
            f"v-annotate.pl "
            f"--mdir {args.vadr_model_dir} "
            f"--mkey flu "
            f"-f "
            f"{fasta_path} "
            f"{sample_outdir}"
        )
        run(vadr_cmd)

        pass_tbl = sample_outdir / f"{sample_name}.vadr.pass.tbl"
        pass_fa  = sample_outdir / f"{sample_name}.vadr.pass.fa"
        if not pass_tbl.is_file():
            print(f"[Warning] {pass_tbl} not found. VADR may have failed to annotate any passing sequences.")
            continue
        if not pass_fa.is_file():
            print(f"[Warning] {pass_fa} not found. VADR may have failed to annotate any passing sequences.")
            continue

        # ---------------------------------------------------------------------
        # 3) PRE_TABLE2ASN: runs sub_seqids_for_table2asn.py
        #    Produces {sample_name}.tbl, {sample_name}.fa, {sample_name}.namesub.txt
        # ---------------------------------------------------------------------
        pre_table2asn_cmd = (
            f"python {args.subseqids_script} "
            f"-t {pass_tbl} "
            f"-f {pass_fa} "
            f"-o {sample_outdir} "
            f"-p {sample_name}"
        )
        run(pre_table2asn_cmd)

        new_tbl      = sample_outdir / f"{sample_name}.tbl"
        new_fasta    = sample_outdir / f"{sample_name}.fa"
        namesub_txt  = sample_outdir / f"{sample_name}.namesub.txt"
        if not new_tbl.is_file() or not new_fasta.is_file():
            print("[Error] sub_seqids_for_table2asn.py did not produce expected .tbl or .fa files.")
            continue

        time.sleep(1)  # short wait if needed for FS latency

        # ---------------------------------------------------------------------
        # 4) TABLE2ASN
        #    By default, writes {sample_name}.gbf, ignoring warnings with '|| true'
        # ---------------------------------------------------------------------
        fa_file = sample_outdir / f"{sample_name}.fa"
        fa_file_abs = str(fa_file.resolve())
        table2asn_cmd = (
            f"table2asn "
            f"-indir {sample_outdir} "
            f"-outdir {sample_outdir} "
            f"-V b -c - "
            f"-i {fa_file_abs} "
            f"|| true"
        )
        run(table2asn_cmd)

        gbf_file = sample_outdir / f"{sample_name}.gbf"
        if not gbf_file.is_file():
            print("[Error] table2asn did not produce a .gbf file as expected.")
            continue

        # ---------------------------------------------------------------------
        # 5) POST_TABLE2ASN
        #    Writes final {sample_name}.gbk, .gff, .faa, .ffn
        # ---------------------------------------------------------------------
        post_table2asn_cmd = (
            f"python {args.posttable2asn_script} "
            f"{gbf_file} "
            f"{namesub_txt} "
            f"-p {sample_name} "
            f"-o {sample_outdir}"
        )
        run(post_table2asn_cmd)

        # Final expected outputs in sample_outdir:
        #   {sample_name}.gbk
        #   {sample_name}.gff
        #   {sample_name}.faa
        #   {sample_name}.ffn
        print(f"[Done] Sample '{sample_name}' completed successfully. Results in '{sample_outdir}'.")

    print("\nAll samples processed. Exiting.\n")

if __name__ == "__main__":
    main()
