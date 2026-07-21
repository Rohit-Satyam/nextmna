// Nextflow DSL2 process converted from Snakemake rule: adjust_danpos_results
process ADJUST_FDR {
    tag "${treatment}_vs_${control}"
    label 'process_adjust_fdr'

    conda "${projectDir}/envs/adjust_fdr.yml"
    publishDir "${params.outdir}/04_differential_nucleosome_positioning/results", mode: 'copy', overwrite: true, pattern: 'adjustedFDR'
    publishDir "${params.outdir}/00_qc/nucleosome_positioning", mode: 'copy', overwrite: true, pattern: 'adjustFDR.log'

    input:
    path results_dir
    path chrom_sizes
    val treatment
    val control
    val n_simulations
    val sampling_occupancy
    val sampling_fuzziness

    output:
    path "adjustedFDR", emit: adjusted
    path "adjustFDR.log", emit: log

    script:
    def simulations_arg = n_simulations ? "-n ${n_simulations}" : ''
    """
    Rscript ${projectDir}/bin/AdjustFDR.R \\
      -o adjustedFDR \\
      -c "${chrom_sizes}" \\
      -s "${treatment}" \\
      -l "${control}" \\
      --sampling_type_occupancy "${sampling_occupancy}" \\
      --sampling_type_fuzziness "${sampling_fuzziness}" \\
      ${simulations_arg} \\
      "${results_dir}" > adjustFDR.log 2>&1

    test -s "adjustedFDR/${treatment}-${control}.positions.FDR_corrected.xls"
    """
}
