// Nextflow DSL2 process converted from Snakemake rule: alignment_bowtie2
process BOWTIE2_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    conda "${projectDir}/envs/bowtie2.yml"
    publishDir "${params.outdir}/00_logs/bowtie2", mode: 'copy', overwrite: true, pattern: '*.bowtie2.log'

    input:
    tuple val(meta), path(r1), path(r2)
    path index_files
    val rnd_seed
    val extra

    output:
    tuple val(meta), path("${meta.id}.sam"), emit: sam
    path "${meta.id}.bowtie2.log", emit: log

    script:
    """
    bowtie2 \\
      -x reference \\
      -1 "${r1}" \\
      -2 "${r2}" \\
      -S ${meta.id}.sam \\
      -p ${task.cpus} \\
      --seed ${rnd_seed} \\
      ${extra} > ${meta.id}.bowtie2.log 2>&1
    """
}
