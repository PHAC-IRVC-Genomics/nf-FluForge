#!/usr/bin/env python3
import subprocess
import sys
import argparse
import time
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
        description="Run Table2asn workflow (PRE, TABLE2ASN, POST) on VADR outputs."
    )
    parser.add_argument(
        "--vadr_outdir",
        required=True,
        help="Directory containing sample subdirectories produced by vadr.py."
    )
    parser.add_argument(
        "--subseqids_script",
        default="sub_seqids_for_table2asn.py",
        help="Path to sub_seqids_for_table2asn.py (PRE_TABLE2ASN script)."
    )
    parser.add_argument(
        "--posttable2asn_script",
        default="post_table2asn.py",
        help="Path to post_table2asn.py (POST_TABLE2ASN script)."
    )
    args = parser.parse_args()

    vadr_outdir = Path(args.vadr_outdir).resolve()
    if not vadr_outdir.is_dir():
        print(f"[Error] vadr_outdir '{vadr_outdir}' is not a directory.")
        sys.exit(1)

    # Iterate over each sample subdirectory (assumed to be created by vadr.py)
    for sample_dir in sorted(vadr_outdir.iterdir()):
        if not sample_dir.is_dir():
            continue
        sample_name = sample_dir.name
        print(f"\n[INFO] Processing sample '{sample_name}' in directory '{sample_dir}'\n")

        # Verify that the VADR output files exist
        pass_tbl = sample_dir / f"{sample_name}.vadr.pass.tbl"
        pass_fa  = sample_dir / f"{sample_name}.vadr.pass.fa"
        if not pass_tbl.is_file() or not pass_fa.is_file():
            print(f"[Warning] Missing VADR output files for sample '{sample_name}'. Skipping.")
            continue

        # ---------------------------------------------------------------------
        # 1) PRE_TABLE2ASN: Run sub_seqids_for_table2asn.py to produce:
        #    {sample_name}.tbl, {sample_name}.fa, {sample_name}.namesub.txt
        # ---------------------------------------------------------------------
        pre_table2asn_cmd = (
            f"python {args.subseqids_script} "
            f"-t {pass_tbl} "
            f"-f {pass_fa} "
            f"-o {sample_dir} "
            f"-p {sample_name}"
        )
        run(pre_table2asn_cmd)

        new_tbl     = sample_dir / f"{sample_name}.tbl"
        new_fasta   = sample_dir / f"{sample_name}.fa"
        namesub_txt = sample_dir / f"{sample_name}.namesub.txt"
        if not new_tbl.is_file() or not new_fasta.is_file():
            print(f"[Error] PRE_TABLE2ASN did not produce expected .tbl or .fa for sample '{sample_name}'. Skipping.")
            continue

        time.sleep(1)  # wait briefly if needed

        # ---------------------------------------------------------------------
        # 2) TABLE2ASN: Run table2asn on the processed files
        # ---------------------------------------------------------------------
        fa_file_abs = str(new_fasta.resolve())
        table2asn_cmd = (
            f"table2asn "
            f"-indir {sample_dir} "
            f"-outdir {sample_dir} "
            f"-V b -c - "  # options; adjust if necessary
            f"-i {fa_file_abs} "
            f"|| true"    # ignore warnings/errors if needed
        )
        run(table2asn_cmd)

        gbf_file = sample_dir / f"{sample_name}.gbf"
        if not gbf_file.is_file():
            print(f"[Error] table2asn did not produce a .gbf file for sample '{sample_name}'. Skipping.")
            continue

        # ---------------------------------------------------------------------
        # 3) POST_TABLE2ASN: Run post_table2asn.py to produce final outputs:
        #    {sample_name}.gbk, .gff, .faa, .ffn
        # ---------------------------------------------------------------------
        post_table2asn_cmd = (
            f"python {args.posttable2asn_script} "
            f"{gbf_file} "
            f"{namesub_txt} "
            f"-p {sample_name} "
            f"-o {sample_dir}"
        )
        run(post_table2asn_cmd)

        print(f"[Done] Sample '{sample_name}' completed successfully. Results in '{sample_dir}'.")

    print("\nAll samples processed with Table2asn workflow. Exiting.\n")

if __name__ == "__main__":
    main()
