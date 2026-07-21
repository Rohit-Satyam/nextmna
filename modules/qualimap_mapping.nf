// Nextflow DSL2 process converted from Snakemake rule: qualimap
process QUALIMAP_MAPPING {
    tag "${meta.id}"
    label 'process_medium'

    conda "${projectDir}/envs/qualimap.yml"
    publishDir "${params.outdir}/00_qc/mapping/qualimap", mode: 'copy', overwrite: true, pattern: '*_qualimap'
    publishDir "${params.outdir}/00_logs/qualimap", mode: 'copy', overwrite: true, pattern: '*.qualimap.log'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_qualimap"), emit: qc_dir
    path "${meta.id}.qualimap.log", emit: log

    script:
    """
    qualimap bamqc \\
      -bam "${bam}" \\
      --java-mem-size=${params.qualimap_java_mem} \\
      -outdir ${meta.id}_qualimap \\
      -c > ${meta.id}.qualimap.log 2>&1
    """
}
