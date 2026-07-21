// Nextflow DSL2 process converted from Snakemake rule: multiqc_filtered
process MULTIQC_FILTERED {
    tag 'filtered'
    label 'process_low'

    conda "${projectDir}/envs/multiqc.yml"
    publishDir "${params.outdir}/00_qc/filtered/multiqc", mode: 'copy', overwrite: true

    input:
    path qc_dirs

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    script:
    """
    multiqc --force --outdir . .
    """
}
