// Nextflow DSL2 process converted from Snakemake rule: sort_filtered_bam
process SORT_FILTERED_BAM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${projectDir}/envs/samtools.yml"
    publishDir "${params.outdir}/03_alignment_filtered/bam", mode: 'copy', overwrite: true, pattern: '*.s.bam'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.s.bam"), emit: sorted_bam

    script:
    """
    samtools sort -m ${params.samtools_sort_memory} -@ ${task.cpus} -o ${meta.id}.s.bam "${bam}"
    """
}
