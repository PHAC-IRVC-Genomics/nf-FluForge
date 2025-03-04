#!/usr/bin/env python3

"""
Script to retrieve canonical Influenza A and Influenza B segment sequences 
from NCBI, rename them, and combine them into one FASTA file (fluAB_refs.fasta).

Requirements:
  - python >= 3.x
  - biopython (pip install biopython or conda install biopython)

Usage:
  python fetch_flu_segments.py

Output:
  fluAB_refs.fasta (contains 16 entries: 8 for Flu A, 8 for Flu B)
"""

from Bio import Entrez, SeqIO

# --- STEP 1: Set your NCBI Entrez email here ---
Entrez.email = "abdallahmeknas@gmail.com"

# --- STEP 2: Define accession lists and mapping to "SegmentName" ---
# For Influenza A/Puerto Rico/8/1934 (H1N1):
influenza_a_segments = {
    "InfluenzaA_PB2": "NC_002016.1",
    "InfluenzaA_PB1": "NC_002017.1",
    "InfluenzaA_PA" : "NC_002018.1",
    "InfluenzaA_HA" : "NC_002019.1",
    "InfluenzaA_NP" : "NC_002020.1",
    "InfluenzaA_NA" : "NC_002021.1",
    "InfluenzaA_M"  : "NC_002022.1",
    "InfluenzaA_NS" : "NC_002023.1"
}

# For Influenza B/Brisbane/60/2008:
influenza_b_segments = {
    "InfluenzaB_PB2": "CY121530.1",
    "InfluenzaB_PB1": "CY121531.1",
    "InfluenzaB_PA" : "CY121532.1",
    "InfluenzaB_HA" : "CY121533.1",
    "InfluenzaB_NP" : "CY121534.1",
    "InfluenzaB_NA" : "CY121535.1",
    "InfluenzaB_M"  : "CY121536.1",
    "InfluenzaB_NS" : "CY121537.1"
}

# Combine all into one dictionary if desired
all_influenza_segments = {}
all_influenza_segments.update(influenza_a_segments)
all_influenza_segments.update(influenza_b_segments)

# --- STEP 3: Function to fetch one sequence from NCBI and rename it ---
def fetch_and_rename(acc: str, new_id: str):
    """
    Fetch a sequence by accession from NCBI nucleotide database
    and rename its record ID to 'new_id'.

    Returns a Bio.SeqRecord object.
    """
    handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
    seq_record = SeqIO.read(handle, "fasta")
    handle.close()

    # rename the sequence record
    seq_record.id = new_id
    # clear out the default description so it doesn't clutter the FASTA
    seq_record.description = ""

    return seq_record

def main():
    output_fasta = "fluAB_refs.fasta"
    print(f"Fetching Influenza A and B segments, output will be in '{output_fasta}'")

    with open(output_fasta, "w") as out_handle:
        for new_id, accession in all_influenza_segments.items():
            print(f"  - Fetching {new_id} from accession {accession} ...")
            seq_record = fetch_and_rename(accession, new_id)
            SeqIO.write(seq_record, out_handle, "fasta")

    print(f"Done. All segments written to '{output_fasta}'.")

if __name__ == "__main__":
    main()
