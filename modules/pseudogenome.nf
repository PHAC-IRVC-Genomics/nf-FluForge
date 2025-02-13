process BUILD_PSEUDOGENOME {
  tag "Building Pseudogenome"

  conda 'conda-forge::python=3.6.15 conda-forge::biopython=1.79'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/biopython:1.79'
  } else {
    container 'quay.io/biocontainers/biopython:1.79'
  }

  input:
  path segment_classifier_output

  output:
  path ('Pseudogenomes'), emit: psuedogenomes_output

  publishDir "${params.output_dir}", mode: 'copy'

  script:
  """
  build_pseudogenome.py --tidy_output ${segment_classifier_output} --output_dir Pseudogenomes
  """
}