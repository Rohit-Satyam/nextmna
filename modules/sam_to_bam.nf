// Nextflow DSL2 process converted from Snakemake rule: sam_to_bam
process SAM_TO_BAM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${projectDir}/envs/samtools.yml"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam

    script:
    """
    samtools view -b -@ ${task.cpus} -o ${meta.id}.bam "${sam}"
    """
}
