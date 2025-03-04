process SANITIZE {
  tag "Sanitize fastq headers"

  conda 'conda-forge::python=3.6.15'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/biopython:1.79'
  } else {
    container 'quay.io/biocontainers/biopython:1.79'
  }

  input:
  path consensus_dir

  output:
  path "sanitized_consensus", emit: sanitized_consensus

  script:
  """
  mkdir -p sanitized_consensus
  sanitize.py -i ${consensus_dir} -o sanitized_consensus
  """
}