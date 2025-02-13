#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {VADR} from './modules/vadr.nf'
include {TABLE2ASN} from './modules/table2asn.nf'
include {SEGMENT_CLASSIFY} from './modules/segment_classifier.nf'
include {BUILD_PSEUDOGENOME} from './modules/pseudogenome.nf'
include {BUILD_TREE} from './modules/tree.nf'

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

    //Segment Classifier Run
    SEGMENT_CLASSIFY(TABLE2ASN.out.table2asn_output_directory)

    //Building Pseudogenomes
    BUILD_PSEUDOGENOME(SEGMENT_CLASSIFY.out.segment_classifier_output)

    //Building Tree
    BUILD_TREE(BUILD_PSEUDOGENOME.out.psuedogenomes_output)

    // Test configuration
    TestConfig()
}