// Nextflow DSL2 process converted from Snakemake rule: index_bam
process INDEX_BAM {
    tag "${meta.id}"
    label 'process_low'

    conda "${projectDir}/envs/samtools.yml"
    publishDir "${params.outdir}/02_alignment/sort", mode: 'copy', overwrite: true, pattern: '*.bai'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("${bam.name}.bai"), emit: bam_bai

    script:
    """
    samtools index -@ ${task.cpus} "${bam}"
    """
}
