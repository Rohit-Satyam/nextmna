// Nextflow DSL2 process converted from Snakemake rule: raw_fastqc
process FASTQC_RAW {
    tag "${meta.id}"
    label 'process_low'

    conda "${projectDir}/envs/fastqc.yml"
    publishDir "${params.outdir}/00_qc/preprocessing/raw_fastqc", mode: 'copy', overwrite: true, pattern: '*_fastqc.*'
    publishDir "${params.outdir}/00_logs/fastqc", mode: 'copy', overwrite: true, pattern: '*.fastqc.log'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta),
        path("${meta.id}_1_fastqc.html"), path("${meta.id}_1_fastqc.zip"),
        path("${meta.id}_2_fastqc.html"), path("${meta.id}_2_fastqc.zip"), emit: qc
    path "${meta.id}.fastqc.log", emit: log

    script:
    """
    fastqc --threads ${task.cpus} "${r1}" "${r2}" > ${meta.id}.fastqc.log 2>&1

    r1_name=\$(basename "${r1}")
    r1_stem=\${r1_name%.gz}
    r1_stem=\${r1_stem%.fastq}
    r1_stem=\${r1_stem%.fq}

    r2_name=\$(basename "${r2}")
    r2_stem=\${r2_name%.gz}
    r2_stem=\${r2_stem%.fastq}
    r2_stem=\${r2_stem%.fq}

    mv "\${r1_stem}_fastqc.html" ${meta.id}_1_fastqc.html
    mv "\${r1_stem}_fastqc.zip"  ${meta.id}_1_fastqc.zip
    mv "\${r2_stem}_fastqc.html" ${meta.id}_2_fastqc.html
    mv "\${r2_stem}_fastqc.zip"  ${meta.id}_2_fastqc.zip
    """
}
