process VADR {
  tag "VADR Run"
  
  conda 'bioconda::vadr=1.6.4'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    // use staphb/vadr Docker container due to issues running VADR Bioconda/Biocontainers Singularity container with Nextflow
    container 'staphb/vadr:1.6.3-hav-flu2'
  } else {
    container 'quay.io/biocontainers/vadr:1.6.4--pl5321h031d066_0'
  }

  input:
  path consensus_dir

  output:
  path ('VADR_output'), emit: vadr_output_directory

  publishDir "${params.output_dir}", mode: 'copy'

  script:
  """
  vadr.py --consensus_dir ${consensus_dir} --vadr_model_dir ${params.vadr_model_dir} --outdir VADR_output
  """
}