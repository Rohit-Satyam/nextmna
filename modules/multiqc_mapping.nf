// Nextflow DSL2 process converted from Snakemake rule: multiqc_mapping
process MULTIQC_MAPPING {
    tag 'mapping'
    label 'process_low'

    conda "${projectDir}/envs/multiqc.yml"
    publishDir "${params.outdir}/00_qc/mapping/multiqc", mode: 'copy', overwrite: true

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
