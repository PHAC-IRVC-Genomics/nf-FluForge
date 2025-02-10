#!/usr/bin/env python3
"""
build_pseudogenome.py

Revised version:
 - DOES NOT blast to classify A/B.
 - Instead, uses the .mdl file (same logic from tidyone.py) to detect influenza A or B.
 - Uses the same segment synonyms from tidyone.py to classify segments.
 - Concatenates segments in the correct order (A vs. B).
 - Optionally runs buildtree.sh to construct trees.

Usage:
  python build_pseudogenome.py \
    --tidy_output /path/to/tidyone_output \
    --output_dir /path/to/pseudo_output \
    --threads 8 \
    --run_tree
"""

import os
import sys
import argparse
import subprocess
import re
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

###############################################################################
# Segment synonyms (from tidyone.py)
###############################################################################
SEGMENT_SYNONYMS_A = {
    "PB2": ["PB2", "POLYMERASE BASIC 2", "SEGMENT1", "SEGMENT 1"],
    "PB1": ["PB1", "POLYMERASE BASIC 1", "SEGMENT2", "SEGMENT 2"],
    "PA":  ["PA", "POLYMERASE ACIDIC", "SEGMENT3", "SEGMENT 3"],
    "HA":  ["HA", "HEMAGGLUTININ", "H1", "H3", "SEGMENT4", "SEGMENT 4"],
    "NP":  ["NP", "NUCLEOPROTEIN", "SEGMENT5", "SEGMENT 5"],
    "NA":  ["NA", "NEURAMINIDASE", "N1", "N2", "SEGMENT6", "SEGMENT 6"],
    "M":   ["MATRIX", "M1", "M2", "BM2", "AM2", "SEGMENT7", "SEGMENT 7"],
    "NS":  ["NS", "NS1", "NS2", "NEP", "NONSTRUCTURAL", "SEGMENT8", "SEGMENT 8"]
}

SEGMENT_SYNONYMS_B = {
    "PB2": ["PB2", "POLYMERASE BASIC 2", "SEGMENT2", "SEGMENT 2"],
    "PB1": ["PB1", "POLYMERASE BASIC 1", "SEGMENT1", "SEGMENT 1"],
    "PA":  ["PA", "POLYMERASE ACIDIC", "SEGMENT3", "SEGMENT 3"],
    "HA":  ["HA", "HEMAGGLUTININ", "H1", "H3", "SEGMENT4", "SEGMENT 4"],
    "NP":  ["NP", "NUCLEOPROTEIN", "SEGMENT5", "SEGMENT 5"],
    "NA":  ["NA", "NEURAMINIDASE", "N1", "N2", "SEGMENT6", "SEGMENT 6"],
    "M":   ["MATRIX", "M1", "M2", "BM2", "AM2", "SEGMENT7", "SEGMENT 7"],
    "NS":  ["NS", "NS1", "NS2", "NEP", "NONSTRUCTURAL", "SEGMENT 8", "SEGMENT 8"]
}

###############################################################################
# Function to detect A/B from .mdl file
###############################################################################
def detect_influenza_type(mdl_file) -> str:
    """
    Check if the .mdl file indicates fluA or fluB.
    If we see 'fluB-' (or 'flub', etc.) anywhere, return 'B'; else default 'A'.
    """
    if not mdl_file.exists():
        print(f"[WARNING] .mdl file not found: {mdl_file}. Defaulting to A.")
        return "A"

    with open(mdl_file, "r", encoding="utf-8") as f:
        for line in f:
            # Replace potential Unicode dashes with ASCII hyphens:
            line = line.replace("–", "-").replace("—", "-")
            # If we see 'fluB-' or 'flu-b' or 'ibv', classify as B
            if "flub-" in line.lower() or "flu-b" in line.lower():
                return "B"
            if "flub" in line.lower() and "seg" in line.lower():
                return "B"
            if "ibv" in line.lower() or "flub" in line.lower():
                return "B"
    return "A"

###############################################################################
# Classify each segment using robust synonyms matching
###############################################################################
def classify_segment(seqrecord, synonyms_dict):
    """
    Return the canonical segment name (PB2, PB1, etc.) based on synonyms,
    using a robust method similar to tidyone.py.

    This function combines the record id and description (uppercased), normalizes
    the text by replacing underscores and hyphens with spaces, and then checks
    for exact token matches (or, for multiword synonyms, regex word-boundary matches).
    """
    # Combine id and description, and convert to uppercase.
    name_str = f"{seqrecord.id} {seqrecord.description}".upper()
    # Normalize by replacing underscores and hyphens with spaces.
    normalized = re.sub(r'[_\-]+', ' ', name_str)
    tokens = normalized.split()

    matches = []
    for segment, synonyms in synonyms_dict.items():
        for syn in synonyms:
            norm_syn = re.sub(r'[_\-]+', ' ', syn.upper()).strip()
            # If the synonym contains a space, use regex with word boundaries.
            if " " in norm_syn:
                pattern = r'\b' + re.escape(norm_syn) + r'\b'
                if re.search(pattern, normalized):
                    matches.append((segment, norm_syn, len(norm_syn)))
            else:
                # For single-token synonyms, check for an exact token match.
                if norm_syn in tokens:
                    matches.append((segment, norm_syn, len(norm_syn)))
    if matches:
        # Choose the match with the longest synonym length as a heuristic.
        matches.sort(key=lambda x: x[2], reverse=True)
        return matches[0][0]
    return "unknown"

###############################################################################
# Build a pseudogenome
###############################################################################
def build_pseudogenome(segment_records, flu_type):
    """
    For Influenza A: order = [PB2, PB1, PA, HA, NP, NA, M_FULL, NS_FULL] (fallback to M, NS)
    For Influenza B: order = [PB1, PB2, PA, HA, NP, NA, M_FULL, NS_FULL] (fallback to M, NS)
    """
    if flu_type == "B":
        segment_order = ["PB1", "PB2", "PA", "HA", "NP", "NA", "M_FULL", "NS_FULL"]
    else:
        segment_order = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M_FULL", "NS_FULL"]

    # A simpler fallback map if we don't have the merged M_FULL or NS_FULL
    fallback_map = {
        "M_FULL": "M",
        "NS_FULL": "NS"
    }

    glued_seq = Seq("")
    glued_desc = []

    for seg in segment_order:
        if seg in segment_records:
            rec = segment_records[seg]
            if isinstance(rec, list):
                rec = rec[0]  # If multiple, just pick the first
            glued_seq += rec.seq
            glued_desc.append(seg)
        else:
            # fallback? e.g. if seg='M_FULL' not found, try 'M'
            if seg in fallback_map:
                fb = fallback_map[seg]
                if fb in segment_records:
                    rec = segment_records[fb]
                    if isinstance(rec, list):
                        rec = rec[0]
                    glued_seq += rec.seq
                    glued_desc.append(fb)
                else:
                    glued_desc.append(f"missing_{seg}")
            else:
                glued_desc.append(f"missing_{seg}")

    new_rec = SeqRecord(glued_seq, id="pseudogenome", description="|".join(glued_desc))
    return new_rec

###############################################################################
# Main pipeline
###############################################################################
def parse_args():
    parser = argparse.ArgumentParser(
        description="Build pseudogenomes from tidyone.py output, using .mdl file to detect Flu A or B."
    )
    parser.add_argument("--tidy_output", required=True,
                        help="Path to the tidyone.py output dir containing each_sample/*.all.segs.cds.fasta etc.")
    parser.add_argument("--output_dir", required=True,
                        help="Where to put pseudogenomes and optional tree results.")
    parser.add_argument("--threads", type=int, default=8,
                        help="Threads for buildtree.sh if --run_tree is set.")
    parser.add_argument("--run_tree", action="store_true",
                        help="If provided, automatically run buildtree.sh on the pseudogenomes.")
    return parser.parse_args()

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def main():
    args = parse_args()

    tidy_output_dir = Path(args.tidy_output)
    each_sample_dir = tidy_output_dir / "each_sample"

    output_dir = Path(args.output_dir)
    ensure_dir(output_dir)

    # Directory for pseudogenome FASTAs
    pseudo_dir = output_dir / "pseudogenomes"
    ensure_dir(pseudo_dir)

    all_pseudo_cds = []
    all_pseudo_pep = []

    # Iterate each sample's .all.segs.cds.fasta
    for cds_file in sorted(each_sample_dir.glob("*.all.segs.cds.fasta")):
        sample_name = cds_file.stem.replace(".all.segs.cds", "")
        pep_file = each_sample_dir / f"{sample_name}.all.segs.pep.fasta"

        # Detect flu type from .mdl
        mdl_file = each_sample_dir / f"{sample_name}.vadr.mdl"
        flu_type = detect_influenza_type(mdl_file)
        print(f"[INFO] Sample: {sample_name}, classified as Flu {flu_type}")

        # Pick synonyms dict
        if flu_type == "B":
            synonyms_dict = SEGMENT_SYNONYMS_B
        else:
            synonyms_dict = SEGMENT_SYNONYMS_A

        # Read CDS segments
        cds_records = list(SeqIO.parse(cds_file, "fasta"))
        if not cds_records:
            print(f"[WARNING] No CDS records in {cds_file}. Skipping.")
            continue

        # Classify each CDS record
        seg_dict_cds = {}
        for rec in cds_records:
            seg_name = classify_segment(rec, synonyms_dict)
            # Check if this is the merged M or NS from tidyone.py
            if "M_FULL_M1_M2" in rec.id.upper():
                seg_name = "M_FULL"
            elif "NS_FULL_NS1_NEP" in rec.id.upper():
                seg_name = "NS_FULL"

            if seg_name not in seg_dict_cds:
                seg_dict_cds[seg_name] = rec
            else:
                existing = seg_dict_cds[seg_name]
                if isinstance(existing, list):
                    existing.append(rec)
                else:
                    seg_dict_cds[seg_name] = [existing, rec]

        # Build pseudogenome (CDS)
        pseudo_cds = build_pseudogenome(seg_dict_cds, flu_type)
        pseudo_cds.id = sample_name
        pseudo_cds.description = f"Influenza_{flu_type}_pseudogenome"
        SeqIO.write(pseudo_cds, pseudo_dir / f"{sample_name}.pseudogenome.cds.fasta", "fasta")
        all_pseudo_cds.append(pseudo_cds)

        # If peptide file exists, do the same
        if pep_file.exists():
            pep_records = list(SeqIO.parse(pep_file, "fasta"))
            seg_dict_pep = {}
            for rec in pep_records:
                seg_name = classify_segment(rec, synonyms_dict)
                if "M_FULL_M1_M2" in rec.id.upper():
                    seg_name = "M_FULL"
                elif "NS_FULL_NS1_NEP" in rec.id.upper():
                    seg_name = "NS_FULL"

                if seg_name not in seg_dict_pep:
                    seg_dict_pep[seg_name] = rec
                else:
                    existing = seg_dict_pep[seg_name]
                    if isinstance(existing, list):
                        existing.append(rec)
                    else:
                        seg_dict_pep[seg_name] = [existing, rec]

            pseudo_pep = build_pseudogenome(seg_dict_pep, flu_type)
            pseudo_pep.id = sample_name
            pseudo_pep.description = f"Influenza_{flu_type}_pseudoProtein"
            SeqIO.write(pseudo_pep, pseudo_dir / f"{sample_name}.pseudogenome.pep.fasta", "fasta")
            all_pseudo_pep.append(pseudo_pep)

    # Write combined multi-sample FASTAs
    combined_cds_path = output_dir / "all_pseudogenomes.cds.fasta"
    combined_pep_path = output_dir / "all_pseudogenomes.pep.fasta"

    if all_pseudo_cds:
        SeqIO.write(all_pseudo_cds, combined_cds_path, "fasta")
    if all_pseudo_pep:
        SeqIO.write(all_pseudo_pep, combined_pep_path, "fasta")

    # Optionally build trees
    if args.run_tree:
        print(f"[INFO] Running buildtree.sh on {combined_cds_path} ...")
        cmd = [
            "bash", "buildtree.sh",
            "-s", str(combined_cds_path.parent),
            "-t", str(args.threads)
        ]
        subprocess.run(cmd, check=True)

    print(f"[DONE] Pseudogenomes in: {pseudo_dir}")

if __name__ == "__main__":
    main()
