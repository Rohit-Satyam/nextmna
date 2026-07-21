// Nextflow DSL2 process converted from Snakemake rule: danpos_diffNuc_prep
process PREPARE_DANPOS_INPUT {
    tag "${meta.condition}:${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.danpos_input"), emit: staged

    script:
    """
    mkdir -p ${meta.id}.danpos_input/${meta.condition}
    cp "${bam}" ${meta.id}.danpos_input/${meta.condition}/${meta.id}.bam
    cp "${bai}" ${meta.id}.danpos_input/${meta.condition}/${meta.id}.bam.bai
    """
}
