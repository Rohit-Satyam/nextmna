// Nextflow DSL2 process converted from Snakemake rule: qualimap_filtered
process QUALIMAP_FILTERED {
    tag "${meta.id}"
    label 'process_medium'

    conda "${projectDir}/envs/qualimap.yml"
    publishDir "${params.outdir}/00_qc/filtered/qualimap", mode: 'copy', overwrite: true, pattern: '*_qualimap'
    publishDir "${params.outdir}/00_logs/qualimap_filtered", mode: 'copy', overwrite: true, pattern: '*.qualimap_filtered.log'

    input:
    tuple val(meta), path(bam), path(bai)
    path annotation

    output:
    tuple val(meta), path("${meta.id}_qualimap"), emit: qc_dir
    path "${meta.id}.qualimap_filtered.log", emit: log

    script:
    """
    qualimap bamqc \\
      -bam "${bam}" \\
      --java-mem-size=${params.qualimap_java_mem} \\
      -gff "${annotation}" \\
      -outdir ${meta.id}_qualimap \\
      -c > ${meta.id}.qualimap_filtered.log 2>&1
    """
}
