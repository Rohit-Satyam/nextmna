// Nextflow DSL2 process converted from Snakemake rule: deeptools_alignmentSieve
process ALIGNMENT_SIEVE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${projectDir}/envs/deeptools.yml"
    publishDir "${params.outdir}/00_qc/filtered/deeptools_alignmentSieve", mode: 'copy', overwrite: true, pattern: '*.alignmentSieve.metrics.txt'
    publishDir "${params.outdir}/00_logs/deeptools_alignmentSieve", mode: 'copy', overwrite: true, pattern: '*.alignmentSieve.log'

    input:
    tuple val(meta), path(bam), path(bai)
    val min_fragment_length
    val max_fragment_length
    val min_mapping_quality
    val sam_flag_include

    output:
    tuple val(meta), path("${meta.id}.filtered.unsorted.bam"), emit: filtered_bam
    path "${meta.id}.alignmentSieve.metrics.txt", emit: metrics
    path "${meta.id}.alignmentSieve.log", emit: log

    script:
    """
    alignmentSieve \\
      -b "${bam}" \\
      --minFragmentLength ${min_fragment_length} \\
      --maxFragmentLength ${max_fragment_length} \\
      --samFlagInclude ${sam_flag_include} \\
      --minMappingQuality ${min_mapping_quality} \\
      -p ${task.cpus} \\
      -o ${meta.id}.filtered.unsorted.bam \\
      --filterMetrics ${meta.id}.alignmentSieve.metrics.txt > ${meta.id}.alignmentSieve.log 2>&1
    """
}
