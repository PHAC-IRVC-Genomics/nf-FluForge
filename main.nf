#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {VADR} from './modules/vadr.nf'
include {TABLE2ASN} from './modules/table2asn.nf'

process TestConfig {
    tag 'Test Config'

    script:
    """
    echo "Testing configuration settings..."
    echo "Input Directory: ${params.consensus_dir}"
    echo "Output Directory: ${params.output_dir}"
    echo "Profile: ${workflow.profile}"
    echo "VADR MODEL DATA: ${params.vadr_model_dir}"
    echo "Max Memory: ${params.max_memory}"
    echo "Max Time: ${params.max_time}"
    echo "Max CPUs: ${params.max_cpus}"
    """
}

workflow {
    // Convert parameters to channels
    ch_versions = Channel.empty()

    // VADR Run
    ch_input = VADR(Channel.fromPath(params.consensus_dir, checkIfExists: true))

    //Table2asn Run
    TABLE2ASN(VADR.out.vadr_output_directory)

    // Test configuration
    TestConfig()
}