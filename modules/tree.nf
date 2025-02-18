process BUILD_TREE {
  tag "Building Tree"

  conda 'conda-forge::mafft=7.526 bioconda::fasttree=2.1.11'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/biopython:1.79'
  } else {
    container 'quay.io/biocontainers/biopython:1.79'
  }

  input:
  path psuedogenomes_output

  output:
  path ('tree'), emit: tree_output

  publishDir "${params.output_dir}", mode: 'copy'

  script:
  """
  build_tree.py -s ${psuedogenomes_output} -t 8
  """
}