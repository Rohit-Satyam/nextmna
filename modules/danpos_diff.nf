// Nextflow DSL2 process converted from Snakemake rule: danpos_diff
process DANPOS_DIFF {
    tag "${treatment}_vs_${control}"
    label 'process_danpos'

    conda "${projectDir}/envs/danpos.yml"
    publishDir "${params.outdir}/04_differential_nucleosome_positioning",
        mode: 'copy',
        overwrite: true,
        saveAs: { name -> name == 'danpos.log' ? null : name }

    publishDir "${params.outdir}/00_qc/nucleosome_positioning",
        mode: 'copy',
        overwrite: true,
        pattern: 'danpos.log'

    input:
    path staged_dirs
    path danpos_code
    val treatment
    val control
    val extra

    output:
    path "results", emit: results
    path "${treatment}", emit: treatment_dir
    path "${control}", emit: control_dir
    path "danpos.log", emit: log

    script:
    """
    mkdir -p "${treatment}" "${control}"

    for staged in *.danpos_input; do
        if [[ -d "\${staged}/${treatment}" ]]; then
            cp "\${staged}/${treatment}/"* "${treatment}/"
        fi
        if [[ -d "\${staged}/${control}" ]]; then
            cp "\${staged}/${control}/"* "${control}/"
        fi
    done

    test "\$(find "${treatment}" -maxdepth 1 -name '*.bam' | wc -l)" -gt 0
    test "\$(find "${control}" -maxdepth 1 -name '*.bam' | wc -l)" -gt 0

    python3 "${danpos_code}/danpos.py" dpos \\
      "${treatment}/:${control}/" \\
      -o results \\
      ${extra} > danpos.log 2>&1

    test -s "results/${treatment}-${control}.positions.integrative.xls"
    """
}
