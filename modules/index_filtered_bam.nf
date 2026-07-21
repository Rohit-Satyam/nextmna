// Nextflow DSL2 process converted from Snakemake rule: index_filtered_bam
process INDEX_FILTERED_BAM {
    tag "${meta.id}"
    label 'process_low'

    conda "${projectDir}/envs/samtools.yml"
    publishDir "${params.outdir}/03_alignment_filtered/bam", mode: 'copy', overwrite: true, pattern: '*.bai'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("${bam.name}.bai"), emit: bam_bai

    script:
    """
    samtools index -@ ${task.cpus} "${bam}"
    """
}
