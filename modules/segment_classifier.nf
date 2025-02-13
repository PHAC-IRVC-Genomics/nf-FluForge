process SEGMENT_CLASSIFY {
  tag "CLassifying Segments"

  conda 'conda-forge::python=3.6.15 conda-forge::biopython=1.79'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/biopython:1.79'
  } else {
    container 'quay.io/biocontainers/biopython:1.79'
  }

  input:
  path table2asn_output_directory

  output:
  path ('flu_segmented_CDS_peptides'), emit: segment_classifier_output

  publishDir "${params.output_dir}", mode: 'copy'

  script:
  """
  segment_classifier.py --input_dir ${table2asn_output_directory} --mode gb --outdir flu_segmented_CDS_peptides --min_orf_len 20

  """
}