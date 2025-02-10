process TABLE2ASN {
  tag "Table2ASN Run"

  conda 'bioconda::table2asn=1.28.943 bioconda::gfflu=0.0.2 conda-forge::python=3.10 conda-forge::biopython=1.80 conda-forge::openpyxl=3.1.0 conda-forge::pandas=1.5.3 conda-forge::rich=12.6.0 conda-forge::typer=0.7.0 conda-forge::xlsxwriter=3.0.8 conda-forge::polars=0.17.9 conda-forge::pyarrow=11.0.0'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container 'https://depot.galaxyproject.org/singularity/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0'
  } else {
    container 'quay.io/biocontainers/mulled-v2-cfa20dfeb068db79c8620a11753add64c23d013a:019cd79f70be602ca625a1a0a4eabab462611a3a-0'
  }

  input:
  path vadr_output_directory

  output:
   path ('VADR_output'), emit: table2asn_output_directory

  publishDir "${params.output_dir}", mode: 'copy'

  script:
  """
  table2asn.py --vadr_outdir ${vadr_output_directory} --subseqids_script ${params.pre_table2asn} --posttable2asn_script ${params.post_table2asn} 

  """
}