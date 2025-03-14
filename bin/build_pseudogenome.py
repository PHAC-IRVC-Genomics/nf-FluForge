#!/usr/bin/env python3
"""
build_pseudogenome.py

Revised version:
 - Does not blast to classify A/B.
 - Instead, uses the .mdl file (same logic from tidyone.py) to detect influenza A or B.
 - Uses the same segment synonyms from tidyone.py to classify segments.
 - Concatenates segments in the correct order (A vs. B).

Usage:
  python build_pseudogenome.py \
    --tidy_output /path/to/tidyone_output \
    --output_dir /path/to/pseudo_output
"""

import os
import sys
import argparse
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
    "NS":  ["NS", "NS1", "NS2", "NEP", "NONSTRUCTURAL", "SEGMENT 8", "SEGMENT 8"]
}

SEGMENT_SYNONYMS_B = {
    "PB2": ["PB2", "POLYMERASE BASIC 2", "SEGMENT2", "SEGMENT 2"],
    "PB1": ["PB1", "POLYMERASE BASIC 1", "SEGMENT1", "SEGMENT 1"],
    "PA":  ["PA", "POLYMERASE ACIDIC", "SEGMENT3", "SEGMENT 3"],
    "HA":  ["HA", "HEMAGGLUTININ", "H1", "H3", "SEGMENT4", "SEGMENT 4"],
    "NP":  ["NP", "NUCLEOPROTEIN", "SEGMENT5", "SEGMENT 5"],
    # NB is included as a synonym for NA:
    "NA":  ["NA", "NB", "NEURAMINIDASE", "N1", "N2", "SEGMENT6", "SEGMENT 6"],
    "M":   ["MATRIX", "M1", "M2", "BM2", "AM2", "SEGMENT7", "SEGMENT 7"],
    "NS":  ["NS", "NS1", "NS2", "NEP", "NONSTRUCTURAL", "SEGMENT 8", "SEGMENT 8"]
}

###############################################################################
def detect_influenza_type(mdl_file) -> str:
    """
    Check if the .mdl file indicates fluA or fluB.
    If we see 'flub-' or 'flu-b' or 'ibv', return 'B'; else default to 'A'.
    """
    if not mdl_file.exists():
        print(f"[WARNING] .mdl file not found: {mdl_file}. Defaulting to A.")
        return "A"
    with open(mdl_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.replace("–", "-").replace("—", "-")
            if "flub-" in line.lower() or "flu-b" in line.lower():
                return "B"
            if "flub" in line.lower() and "seg" in line.lower():
                return "B"
            if "ibv" in line.lower() or "flub" in line.lower():
                return "B"
    return "A"

def classify_segment(seqrecord, synonyms_dict):
    """
    Return the canonical segment name (PB2, PB1, etc.) based on synonyms,
    using a robust method similar to tidyone.py.
    """
    name_str = f"{seqrecord.id} {seqrecord.description}".upper()
    normalized = re.sub(r'[_\-]+', ' ', name_str)
    tokens = normalized.split()
    matches = []
    for segment, synonyms in synonyms_dict.items():
        for syn in synonyms:
            norm_syn = re.sub(r'[_\-]+', ' ', syn.upper()).strip()
            if " " in norm_syn:
                pattern = r'\b' + re.escape(norm_syn) + r'\b'
                if re.search(pattern, normalized):
                    matches.append((segment, norm_syn, len(norm_syn)))
            else:
                if norm_syn in tokens:
                    matches.append((segment, norm_syn, len(norm_syn)))
    if matches:
        matches.sort(key=lambda x: x[2], reverse=True)
        return matches[0][0]
    return "unknown"

def build_pseudogenome(segment_records, flu_type):
    """
    For Influenza A: order = [PB2, PB1, PA, HA, NP, NA, M_FULL, NS_FULL] (fallback to M, NS)
    For Influenza B: order = [PB1, PB2, PA, HA, NP, NA, M_FULL, NS_FULL] (fallback to M, NS)
    """
    if flu_type == "B":
        segment_order = ["PB1", "PB2", "PA", "HA", "NP", "NA", "M_FULL", "NS_FULL"]
    else:
        segment_order = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M_FULL", "NS_FULL"]

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
                rec = rec[0]  # If multiple, pick the first
            glued_seq += rec.seq
            glued_desc.append(seg)
        else:
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

def parse_args():
    parser = argparse.ArgumentParser(
        description="Build pseudogenomes from tidyone.py output, using .mdl file to detect Flu A or B."
    )
    parser.add_argument("--tidy_output", required=True,
                        help="Path to the tidyone.py output directory containing each_sample/*.all.segs.cds.fasta etc.")
    parser.add_argument("--output_dir", required=True,
                        help="Directory where pseudogenomes will be saved.")
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

    pseudo_dir = output_dir / "pseudogenomes"
    ensure_dir(pseudo_dir)

    all_pseudo_cds = []
    all_pseudo_pep = []

    # Process each sample's CDS FASTA
    for cds_file in sorted(each_sample_dir.glob("*.all.segs.cds.fasta")):
        sample_name = cds_file.stem.replace(".all.segs.cds", "")
        pep_file = each_sample_dir / f"{sample_name}.all.segs.pep.fasta"

        # Detect flu type using .mdl file
        mdl_file = each_sample_dir / f"{sample_name}.vadr.mdl"
        flu_type = detect_influenza_type(mdl_file)
        print(f"[INFO] Sample: {sample_name}, classified as Flu {flu_type}")

        synonyms_dict = SEGMENT_SYNONYMS_B if flu_type == "B" else SEGMENT_SYNONYMS_A

        cds_records = list(SeqIO.parse(cds_file, "fasta"))
        if not cds_records:
            print(f"[WARNING] No CDS records in {cds_file}. Skipping.")
            continue

        seg_dict_cds = {}
        for rec in cds_records:
            seg_name = classify_segment(rec, synonyms_dict)
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

        # Process peptide records if available
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

    combined_cds_path = output_dir / "all_pseudogenomes.cds.fasta"
    combined_pep_path = output_dir / "all_pseudogenomes.pep.fasta"

    if all_pseudo_cds:
        SeqIO.write(all_pseudo_cds, combined_cds_path, "fasta")
    if all_pseudo_pep:
        SeqIO.write(all_pseudo_pep, combined_pep_path, "fasta")

    print(f"[DONE] Pseudogenomes saved in: {pseudo_dir}")

    # Create a CSV file with the lengths of all pseudogenomes produced.
    lengths_csv = output_dir / "pseudogenome_lengths.csv"
    with open(lengths_csv, "w") as out_csv:
        out_csv.write("Sample,CDS_length,PEP_length\n")
        # Build a dictionary for peptide lengths if available
        pep_dict = { rec.id: len(rec.seq) for rec in all_pseudo_pep } if all_pseudo_pep else {}
        for rec in all_pseudo_cds:
            sample = rec.id
            cds_len = len(rec.seq)
            pep_len = pep_dict.get(sample, "NA")
            out_csv.write(f"{sample},{cds_len},{pep_len}\n")
    print(f"[INFO] Pseudogenome lengths saved in: {lengths_csv}")

if __name__ == "__main__":
    main()
