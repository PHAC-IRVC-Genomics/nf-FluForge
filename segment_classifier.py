#!/usr/bin/env python3

"""
tidyone.py

A Python pipeline to:

1) Parse multiple influenza assemblies (FASTA) + annotation (GenBank or GFF).
2) Identify each of the 8 segments (PB2, PB1, PA, HA, NP, NA, M, NS).
3) Extract coding regions (CDS) for each segment, including spliced transcripts (e.g. M2, NEP).
4) Translate them to peptide sequences.
5) Merge M1+M2 into M_FULL_M1_M2, and merge NS1+NEP into NS_FULL_NS1_NEP (ignoring min length).
6) Output for each sample into structured directories, and also produce combined files.

Requires:
  Biopython

Example directory structure (input_dir = "batch_pp/"):

  Samples
  ├── RV04983-24-IAV
  │   ├── RV04983-24-IAV.vadr.pass.fa    <-- FASTA
  │   ├── RV04983-24-IAV.gff            <-- GFF annotation (if --mode gff)
  │   ├── RV04983-24-IAV.gbk            <-- GenBank annotation (if --mode gb)
  │   ├── RV04983-24-IAV.vadr.mdl       <-- VADR mdl file
  │   ...
  ├── RV05560-24-IBV
  ├── sample1
  └── sample2

Usage:
  python tidyone.py \
    --input_dir batch_pp \
    --outdir my_output_dir \
    --mode gff  (or --mode gb)
"""

import os
import sys
import argparse
import shutil
import re
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

##############################################################################
# Dictionaries of segment synonyms for Flu A vs. Flu B
##############################################################################
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
    "NS":  ["NS", "NS1", "NS2", "NEP", "NONSTRUCTURAL", "SEGMENT8", "SEGMENT 8"]
}
##############################################################################

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="A Python pipeline for extracting influenza segments and spliced products from GenBank/GFF."
    )
    parser.add_argument("--input_dir", required=True,
                        help="Directory containing subdirectories for each sample, each with FASTA/annotation files.")
    parser.add_argument("--mode", choices=["gb", "gff"], default="gb",
                        help="Annotation file type: 'gb' for GenBank (.gb/.gbk), 'gff' for GFF (VADR output).")
    parser.add_argument("--outdir", required=True,
                        help="Output directory where results will be saved.")
    parser.add_argument("--min_orf_len", type=int, default=500,
                        help="Minimum CDS length to consider (default=500).")
    parser.add_argument("--max_orf_len", type=int, default=10000,
                        help="Maximum CDS length to consider (default=10000).")
    return parser.parse_args()

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def extract_cds_from_genbank(gb_file, min_len, max_len):
    from Bio.SeqFeature import CompoundLocation

    all_cds_records = []
    for record in SeqIO.parse(gb_file, "genbank"):
        for feat in record.features:
            if feat.type == "CDS":
                if isinstance(feat.location, CompoundLocation):
                    start = min(p.start.position for p in feat.location.parts)
                    end   = max(p.end.position   for p in feat.location.parts)
                else:
                    start = feat.location.start.position
                    end   = feat.location.end.position

                loc_seq = feat.location.extract(record.seq)
                if len(loc_seq) < min_len or len(loc_seq) > max_len:
                    continue

                qualifiers = feat.qualifiers
                gene_name = qualifiers.get("gene", ["unknown_gene"])[0]
                product   = qualifiers.get("product", ["unknown_product"])[0]

                new_id = f"{record.id}_{gene_name}"
                description = product

                new_rec = SeqRecord(
                    loc_seq,
                    id=new_id,
                    description=description
                )
                new_rec.annotations["cds_coords"] = (start, end)
                new_rec.annotations["parent_seq"] = record.id

                all_cds_records.append(new_rec)
    return all_cds_records

def extract_cds_from_gff(fasta_path, gff_path, min_len, max_len):
    from Bio.SeqFeature import CompoundLocation, FeatureLocation

    seq_record = next(SeqIO.parse(fasta_path, "fasta"))
    all_cds_records = []
    features_by_id = {}

    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            ftype = cols[2]
            if ftype.upper() == "CDS":
                start = int(cols[3]) - 1  # GFF is 1-based
                end   = int(cols[4])     # inclusive
                strand = cols[6]
                attributes = cols[8]

                attr_dict = {}
                for kv in attributes.split(";"):
                    if "=" in kv:
                        key, val = kv.split("=", 1)
                        attr_dict[key] = val

                gene_name = attr_dict.get("gene", "unknown_gene")
                product   = attr_dict.get("product", "unknown_product")
                feature_id = attr_dict.get("ID", None)

                if feature_id not in features_by_id:
                    features_by_id[feature_id] = {
                        "gene": gene_name,
                        "product": product,
                        "parts": [],
                        "strand": strand,
                    }
                features_by_id[feature_id]["parts"].append((start, end))

    for feat_id, info in features_by_id.items():
        gene_name = info["gene"]
        product   = info["product"]
        strand    = info["strand"]
        parts     = sorted(info["parts"], key=lambda x: x[0])

        loc_objs = []
        for (s, e) in parts:
            floc = FeatureLocation(s, e, strand=(1 if strand == "+" else -1))
            loc_objs.append(floc)

        if len(loc_objs) == 1:
            full_loc = loc_objs[0]
            start = full_loc.start.position
            end   = full_loc.end.position
        else:
            full_loc = CompoundLocation(loc_objs)
            start = min(p.start for p in loc_objs).position
            end   = max(p.end for p in loc_objs).position

        sub_seq = full_loc.extract(seq_record.seq)
        if len(sub_seq) < min_len or len(sub_seq) > max_len:
            continue

        new_id = f"{seq_record.id}_{gene_name}"
        new_rec = SeqRecord(sub_seq, id=new_id, description=product)
        new_rec.annotations["cds_coords"] = (start, end)
        new_rec.annotations["parent_seq"] = seq_record.id

        all_cds_records.append(new_rec)

    return all_cds_records

def translate_seqrecord(seqrecord):
    prot_seq = seqrecord.seq.translate(to_stop=True)
    return SeqRecord(prot_seq, id=seqrecord.id, description=seqrecord.description)

def detect_influenza_type(mdl_file) -> str:
    """
    Check if the .mdl file indicates fluA or fluB.
    If we see 'fluB-' anywhere in the lines for #idx lines, we return 'B';
    else default to 'A'.
    """
    with open(mdl_file, "r") as f:
        for line in f:
            # skip commented lines if you want, but this example shows
            # lines might begin with '#' 
            if "flub" in line.lower() and "seg" in line.lower():
                return "B"
    # Default to A if not found
    return "A"

def classify_segment(seqrecord, synonyms_dict):
    """
    Classify which influenza segment a given CDS belongs to using a robust method.
    
    This version first combines the record id and description, converts to uppercase,
    and then normalizes by replacing underscores and hyphens with spaces. It then checks
    for exact word matches against the provided synonyms. If a synonym is made up of multiple
    words (or includes a space), a regex search with word boundaries is used. In cases of 
    multiple matches, the synonym with the longest length is chosen.
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

def merge_m1_m2(seg_dict_cds, original_sample_fasta):
    if "M" not in seg_dict_cds:
        return seg_dict_cds

    m_records = seg_dict_cds["M"]
    if not m_records:
        return seg_dict_cds

    m1_recs = [r for r in m_records if "M1" in r.id.upper() or "M1" in r.description.upper()]
    m2_recs = [r for r in m_records if "M2" in r.id.upper() or "M2" in r.description.upper()]

    if not m1_recs or not m2_recs:
        return seg_dict_cds

    starts = []
    ends = []
    for r in (m1_recs + m2_recs):
        s, e = r.annotations.get("cds_coords", (None, None))
        if s is not None and e is not None:
            starts.append(s)
            ends.append(e)

    if not starts or not ends:
        return seg_dict_cds

    merged_start = min(starts)
    merged_end   = max(ends)

    parent_id = m1_recs[0].annotations["parent_seq"]
    genome_record = None
    for rec in SeqIO.parse(original_sample_fasta, "fasta"):
        if rec.id == parent_id:
            genome_record = rec
            break

    if not genome_record:
        print(f"[WARNING] Could not find parent_seq {parent_id} in {original_sample_fasta}. Skipping M_FULL.")
        return seg_dict_cds

    merged_nt = genome_record.seq[merged_start:merged_end]
    merged_rec = SeqRecord(
        merged_nt,
        id=f"{genome_record.id}_M_FULL_M1_M2",
        description="Matrix FULL (M1+M2)"
    )
    merged_rec.annotations["cds_coords"] = (merged_start, merged_end)
    seg_dict_cds["M"].append(merged_rec)

    return seg_dict_cds

def merge_ns1_nep(seg_dict_cds, original_sample_fasta):
    if "NS" not in seg_dict_cds:
        return seg_dict_cds

    ns_records = seg_dict_cds["NS"]
    if not ns_records:
        return seg_dict_cds

    ns1_recs = [r for r in ns_records if "NS1" in r.id.upper() or "NS1" in r.description.upper()]
    nep_recs = [r for r in ns_records if ("NEP" in r.id.upper() or "NS2" in r.id.upper() or
                                          "NEP" in r.description.upper() or "NS2" in r.description.upper())]

    if not ns1_recs or not nep_recs:
        return seg_dict_cds

    starts = []
    ends   = []
    for r in (ns1_recs + nep_recs):
        s, e = r.annotations.get("cds_coords", (None, None))
        if s is not None and e is not None:
            starts.append(s)
            ends.append(e)

    if not starts or not ends:
        return seg_dict_cds

    merged_start = min(starts)
    merged_end   = max(ends)

    parent_id = ns1_recs[0].annotations["parent_seq"]
    genome_record = None
    for rec in SeqIO.parse(original_sample_fasta, "fasta"):
        if rec.id == parent_id:
            genome_record = rec
            break

    if not genome_record:
        print(f"[WARNING] Could not find parent_seq {parent_id} in {original_sample_fasta}. Skipping NS_FULL.")
        return seg_dict_cds

    merged_nt = genome_record.seq[merged_start:merged_end]
    merged_rec = SeqRecord(
        merged_nt,
        id=f"{genome_record.id}_NS_FULL_NS1_NEP",
        description="Nonstructural FULL (NS1+NEP)"
    )
    merged_rec.annotations["cds_coords"] = (merged_start, merged_end)
    seg_dict_cds["NS"].append(merged_rec)

    return seg_dict_cds

def find_annotation_file(sample_dir, sample_name, mode):
    """
    Search for the annotation file in the given 'sample_dir' that matches sample_name
    with either .gb/.gbk if mode='gb', or .gff if mode='gff'.
    Returns the Path or None if not found.
    """
    if mode == "gb":
        candidates = list(sample_dir.glob(f"{sample_name}*.gb")) + \
                     list(sample_dir.glob(f"{sample_name}*.gbk"))
        return candidates[0] if candidates else None
    else:
        candidate = sample_dir / f"{sample_name}.gff"
        if candidate.exists():
            return candidate
        # or .gbf, if you have that extension
        candidate_gbf = sample_dir / f"{sample_name}.gbf"
        if candidate_gbf.exists():
            return candidate_gbf
        return None

def main():
    args = parse_arguments()

    input_dir = Path(args.input_dir)
    outdir = Path(args.outdir)
    ensure_dir(outdir)

    each_sample_dir = outdir / "each_sample"
    segs_dir        = outdir / "segs"
    tree_file_dir   = outdir / "tree_file"

    for d in [each_sample_dir, segs_dir, tree_file_dir]:
        ensure_dir(d)

    combined_cds_records = []
    combined_pep_records = []

    # ---------------------------------------------------------------------
    # MAIN LOOP: For each subdirectory in input_dir, we look for the sample's files
    # ---------------------------------------------------------------------
    for sample_subdir in sorted(input_dir.iterdir()):
        if not sample_subdir.is_dir():
            continue

        sample_name = sample_subdir.name

        # Attempt to detect if it's fluA or fluB via .mdl file
        mdl_file = sample_subdir / f"{sample_name}.vadr.mdl"
        if mdl_file.exists():
            # parse the mdl file to see if 'fluB' is mentioned
            flu_type = detect_influenza_type(mdl_file)
            # also copy the mdl file to each_sample_dir for future reference
            shutil.copy(mdl_file, each_sample_dir / f"{sample_name}.vadr.mdl")
        else:
            # If there's no .mdl, default to "A"
            flu_type = "A"

        # Choose which synonyms dict to use
        if flu_type == "B":
            synonyms_dict = SEGMENT_SYNONYMS_B
        else:
            synonyms_dict = SEGMENT_SYNONYMS_A

        # Find the FASTA file
        fasta_file = sample_subdir / f"{sample_name}.vadr.pass.fa"
        if not fasta_file.exists():
            alt_fasta_file = sample_subdir / f"{sample_name}.fa"
            if alt_fasta_file.exists():
                fasta_file = alt_fasta_file
            else:
                print(f"[WARNING] No FASTA found in {sample_subdir}. Skipping.")
                continue

        # Identify the annotation file (gbk or gff)
        ann_file = find_annotation_file(sample_subdir, sample_name, args.mode)
        if not ann_file:
            print(f"[WARNING] No annotation (.gb/.gbk or .gff) found for sample '{sample_name}' in {sample_subdir}.")
            continue

        # Extract CDS
        if args.mode == "gb":
            cds_list = extract_cds_from_genbank(ann_file, args.min_orf_len, args.max_orf_len)
        else:
            cds_list = extract_cds_from_gff(fasta_file, ann_file, args.min_orf_len, args.max_orf_len)

        # Translate them all
        pep_list = [translate_seqrecord(c) for c in cds_list]

        # Classify segments
        seg_dict_cds = {}
        seg_dict_pep = {}
        for cds_rec, pep_rec in zip(cds_list, pep_list):
            seg_name = classify_segment(cds_rec, synonyms_dict)
            seg_dict_cds.setdefault(seg_name, []).append(cds_rec)
            seg_dict_pep.setdefault(seg_name, []).append(pep_rec)

        # Merge M1+M2
        seg_dict_cds = merge_m1_m2(seg_dict_cds, fasta_file)
        # Merge NS1+NEP
        seg_dict_cds = merge_ns1_nep(seg_dict_cds, fasta_file)

        # Re-translate if M or NS changed
        for seg_name in ["M", "NS"]:
            if seg_name in seg_dict_cds:
                new_peps = [translate_seqrecord(c) for c in seg_dict_cds[seg_name]]
                seg_dict_pep[seg_name] = new_peps

        # Prepare output filenames
        sample_cds_path = each_sample_dir / f"{sample_name}.all.segs.cds.fasta"
        sample_pep_path = each_sample_dir / f"{sample_name}.all.segs.pep.fasta"

        with open(sample_cds_path, "w") as fcds, open(sample_pep_path, "w") as fpep:
            # Write CDS
            for seg_name, recs in seg_dict_cds.items():
                for rec in recs:
                    SeqIO.write(rec, fcds, "fasta")
            # Write peptides
            for seg_name, recs in seg_dict_pep.items():
                for rec in recs:
                    SeqIO.write(rec, fpep, "fasta")

        # Collect for combined output
        all_cds = [r for seg in seg_dict_cds.values() for r in seg]
        all_pep = [r for seg in seg_dict_pep.values() for r in seg]

        combined_cds_records.extend(all_cds)
        combined_pep_records.extend(all_pep)

        # Also append to segment-level FASTAs
        for seg_name, recs in seg_dict_cds.items():
            seg_cds_file = segs_dir / f"{seg_name}.cds.fasta"
            seg_pep_file = segs_dir / f"{seg_name}.pep.fasta"
            with open(seg_cds_file, "a") as fsegcds, open(seg_pep_file, "a") as fsegpep:
                for rec_cds in recs:
                    SeqIO.write(rec_cds, fsegcds, "fasta")
                for rec_pep in seg_dict_pep[seg_name]:
                    SeqIO.write(rec_pep, fsegpep, "fasta")

    # Build a combined tree file of *all* samples' CDS/peptides
    all_cds_path = tree_file_dir / "all.samples.cds.fasta"
    all_pep_path = tree_file_dir / "all.samples.pep.fasta"

    with open(all_cds_path, "w") as f_all_cds:
        SeqIO.write(combined_cds_records, f_all_cds, "fasta")
    with open(all_pep_path, "w") as f_all_pep:
        SeqIO.write(combined_pep_records, f_all_pep, "fasta")

    print("Done. Output in:", outdir)

if __name__ == "__main__":
    main()
