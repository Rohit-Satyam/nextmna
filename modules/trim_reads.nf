process TRIM_READS {
    tag "${meta.id}"
    label 'process_low'

    conda "${projectDir}/envs/fastp.yml"

    publishDir "${params.outdir}/01_trimmed_reads",
        mode: 'copy',
        overwrite: true,
        pattern: '*-trimmed-pair*.fastq.gz'

    publishDir "${params.outdir}/00_qc/fastp",
        mode: 'copy',
        overwrite: true,
        pattern: '*.fastp.*'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta),
        path("${meta.id}-trimmed-pair1.fastq.gz"),
        path("${meta.id}-trimmed-pair2.fastq.gz"),
        emit: reads

    path "${meta.id}.fastp.json", emit: json
    path "${meta.id}.fastp.html", emit: html
    path "${meta.id}.fastp.log",  emit: log

    script:
    def extra = params.fastp_extra ?: ''

    """
    fastp \\
      --in1 "${r1}" \\
      --in2 "${r2}" \\
      --out1 "${meta.id}-trimmed-pair1.fastq.gz" \\
      --out2 "${meta.id}-trimmed-pair2.fastq.gz" \\
      --thread ${task.cpus} \\
      --json "${meta.id}.fastp.json" \\
      --html "${meta.id}.fastp.html" \\
      ${extra} \\
      > "${meta.id}.fastp.log" 2>&1

    test -s "${meta.id}-trimmed-pair1.fastq.gz"
    test -s "${meta.id}-trimmed-pair2.fastq.gz"
    test -s "${meta.id}.fastp.json"
    test -s "${meta.id}.fastp.html"
    """
}
