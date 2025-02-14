# nf-FluForge

**nf-FluForge** is a Nextflow-based pipeline for comprehensive influenza virus genomics analysis. The pipeline integrates multiple modules to:

- **Annotate** consensus FASTA files using VADR.
- **Process** VADR outputs with Table2asn to standardize annotations.
- **Extract and classify** influenza segments (e.g., PB2, PB1, PA, HA, NP, NA, M, NS) with a robust segment classifier.
- **Build pseudogenomes** by concatenating the segmented coding sequences in the correct order.
- **Construct phylogenetic trees** from the pseudogenomes using MAFFT and FastTree.

This modular design ensures reproducible and scalable analysis of influenza assemblies.

---

## Prerequisites

- **Conda**  
  nf-FluForge can be run using the Conda profile (`-profile conda`). Other profiles are coming soon.

---

## Quick Start

### Requirements

- **Nextflow**: Version >=22.10.1  
- **Container or Environment**:
  - Conda
  or
  - Singularity (Coming soon)  
  or
  - Docker (Coming soon)

- **Data Inputs**:
  - A directory with consensus FASTA files.

---

## Running nf-FluForge

With your Conda environment activated, run the pipeline using Nextflow. The main workflow is defined in `main.nf`.

### Usage Example

```bash
nextflow run main.nf \
  --consensus_dir /path/to/consensus_fastas \
  --output_dir /path/to/output_directory \
  -profile conda
```

---

## Parameters

- `--consensus_dir`  
  Directory containing consensus FASTA files (one per sample).

- `--output_dir`  
  Directory where the pipeline output will be stored.

- `--max_memory`  
  Maximum memory allocation per process (e.g., `8.GB`).

- `--max_time`  
  Maximum runtime per process (e.g., `2.h`).

- `--max_cpus`  
  Maximum CPUs per process (e.g., `4`).

---

## Pipeline Modules

The pipeline executes the following modules in order:

1. **VADR**: Annotates input FASTA files.
2. **TABLE2ASN**: Processes VADR outputs.
3. **SEGMENT_CLASSIFY**: Extracts and classifies influenza segments.
4. **BUILD_PSEUDOGENOME**: Constructs pseudogenomes.
5. **BUILD_TREE**: Generates phylogenetic trees from the pseudogenomes.

A configuration test process (`TestConfig`) is also included to verify the parameter settings.

---

## Output File Structure and Explanation

After the pipeline completes, your output directory will contain several subdirectories corresponding to each processing stage. An example structure is shown below:

```
/path/to/output_directory/
├── VADR_output/                   
│   ├── Sample1/
│   │   ├── Sample1.vadr.pass.tbl
│   │   ├── Sample1.vadr.pass.fa
│   │   ├── Sample1.gbf
│   │   ├── Sample1.gbk
│   │   ├── Sample1.gff
│   │   ├── Sample1.faa
│   │   └── Sample1.ffn
│   │   └── ...
│   └── Sample2/
│       ├── Sample2.vadr.pass.tbl
│       ├── Sample2.vadr.pass.fa
│   │   ├── Sample2.gbf
│   │   ├── Sample2.gbk
│   │   ├── Sample2.gff
│   │   ├── Sample2.faa
│   │   └── Sample2.ffn
│       └── ...
├── flu_segmented_CDS_peptides/     
│   ├── each_sample/               
│   │   ├── Sample1.all.segs.cds.fasta
│   │   ├── Sample1.all.segs.pep.fasta
│   │   └── ...
│   ├── segs/                      
│   │   ├── PB2.cds.fasta
│   │   ├── PB2.pep.fasta
│   │   └── ...
│   └── all.samples.*.fasta        
├── pseudogenomes/                 
│   ├── Sample1.pseudogenome.cds.fasta
│   ├── Sample1.pseudogenome.pep.fasta
│   └── all_pseudogenomes.*.fasta   
└── tree/                         
    ├── Sample1.pseudogenome.cds.tre
    ├── Sample1.pseudogenome.pep.tre
    └── ...
```

### Directory Explanations

- **`VADR_output`**  
  Contains VADR annotation results for each sample (e.g., `.tbl` and `.fa` files) and final outputs from the Table2asn process (e.g., `.gbf`, `.gbk`, `.gff`, `.faa`, `.ffn`) .

- **`flu_segmented_CDS_peptides`**  
  - `each_sample/`: Per-sample FASTA files containing all extracted CDS and peptide sequences.  
  - `segs/`: Segment-specific FASTA files (e.g., for PB2, PB1, etc.).  
  - `all.samples.*.fasta`: Combined multi-sample FASTA files for downstream analyses.  

- **`pseudogenomes`**  
  Contains pseudogenome sequences (both CDS and peptides) built by concatenating the classified segments in the correct order. Also includes combined multi-sample pseudogenome files.

- **`tree`**  
  Contains phylogenetic tree files (in Newick format, typically with a `.tre` extension) generated from the pseudogenomes using MAFFT and FastTree.

---

## License

Copyright Government of Canada 2025

Written by: National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

---

## Acknowledgments

IRVC 

---

## Contact

For questions or suggestions, please contact **Abdallah Meknas** at **[abdallah.meknas@phac-aspc.gc.ca]**.
