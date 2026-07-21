process MERGE_SORT_NUCDYN_BAM {
    tag "${meta.id}"
    label 'process_medium'

    conda 'bioconda::samtools=1.21'

    publishDir "${params.outdir}/NucleosomeDynamics/bam", mode: 'copy'

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("${meta.id}.s.bam"), emit: sorted_bam

    script:
    def bam_args = bams.collect { bam -> "'${bam}'" }.join(' ')
    """
    samtools merge \
        -@ ${task.cpus} \
        -f \
        '${meta.id}.bam' \
        ${bam_args}

    samtools sort \
        -@ ${task.cpus} \
        -m '${params.samtools_sort_memory}' \
        -o '${meta.id}.s.bam' \
        '${meta.id}.bam'

    samtools quickcheck -v '${meta.id}.s.bam'
    """
}


process NUCDYN_READ_BAM {
    tag "${meta.id}"
    label 'process_nucdyn'

    container "${params.nucdyn_container}"
    // The image has its own runNuclDyn ENTRYPOINT. Nextflow needs Bash instead.
    containerOptions '--entrypoint "" -u $(id -u):$(id -g)'

    publishDir "${params.outdir}/NucleosomeDynamics/RData", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.RData"), emit: rdata

    script:
    """
    export HOME=/tmp

    '${params.nucdyn_launcher}' readBAM \
        --input '${bam}' \
        --output '${meta.id}.RData' \
        --type paired \
        2>&1 | tee '${meta.id}.readBAM.log'

    if grep -q 'ERROR while running' '${meta.id}.readBAM.log'; then
        exit 1
    fi

    test -s '${meta.id}.RData'
    """
}


process NUCDYN_NUCLER {
    tag "${meta.id}"
    label 'process_nucdyn'

    container "${params.nucdyn_container}"
    containerOptions '--entrypoint "" -u $(id -u):$(id -g)'

    publishDir "${params.outdir}/NucleosomeDynamics/nucleR", mode: 'copy'

    input:
    tuple val(meta), path(rdata)
    val fragment_length
    val threshold_percentage

    output:
    tuple val(meta), path("${meta.id}.gff"), emit: calls

    script:
    """
    export HOME=/tmp

    '${params.nucdyn_launcher}' nucleR \
        --input '${rdata}' \
        --output '${meta.id}.gff' \
        --type paired \
        --fragmentLen '${fragment_length}' \
        --thresholdPercentage '${threshold_percentage}' \
        2>&1 | tee '${meta.id}.nucleR.log'

    if grep -q 'ERROR while running' '${meta.id}.nucleR.log'; then
        exit 1
    fi

    test -s '${meta.id}.gff'
    """
}


process NUCDYN_NUCLER_STATS {
    tag "${meta.id}"
    label 'process_nucdyn'

    container "${params.nucdyn_container}"
    containerOptions '--entrypoint "" -u $(id -u):$(id -g)'

    publishDir "${params.outdir}/NucleosomeDynamics/statistics/nucleR", mode: 'copy'

    input:
    tuple val(meta), path(calls)
    path annotation

    output:
    tuple val(meta), path("${meta.id}_genes.csv"), path("${meta.id}_genome.csv"), emit: stats

    script:
    """
    export HOME=/tmp

    '${params.nucdyn_launcher}' nucleR_stats \
        --input '${calls}' \
        --genome '${annotation}' \
        --out_genes '${meta.id}_genes.csv' \
        --out_gw '${meta.id}_genome.csv' \
        2>&1 | tee '${meta.id}.nucleR_stats.log'

    if grep -q 'ERROR while running' '${meta.id}.nucleR_stats.log'; then
        exit 1
    fi

    test -s '${meta.id}_genes.csv'
    test -s '${meta.id}_genome.csv'
    """
}


process NUCDYN_TXSTART {
    tag "${meta.id}"
    label 'process_nucdyn'

    container "${params.nucdyn_container}"
    containerOptions '--entrypoint "" -u $(id -u):$(id -g)'

    publishDir "${params.outdir}/NucleosomeDynamics/txstart",
        mode: 'copy'

    input:
    tuple val(meta), path(calls)
    path annotation

    output:
    tuple val(meta),
          path("${meta.id}.txstart.gff"),
          emit: txstart

    script:
    def output_gff = "${meta.id}.txstart.gff"

    """
    export HOME=/tmp

    '${params.nucdyn_launcher}' txstart \
        --calls '${calls}' \
        --genome '${annotation}' \
        --output '${output_gff}' \
        2>&1 | tee '${meta.id}.txstart.log'

    if grep -q 'ERROR while running' '${meta.id}.txstart.log'; then
        exit 1
    fi

    test -s '${output_gff}'
    """
}


process NUCDYN_TXSTART_STATS {
    tag "${meta.id}"
    label 'process_nucdyn'

    container "${params.nucdyn_container}"
    containerOptions '--entrypoint "" -u $(id -u):$(id -g)'

    publishDir "${params.outdir}/NucleosomeDynamics/statistics/txstart", mode: 'copy'

    input:
    tuple val(meta), path(txstart)
    path annotation

    output:
    tuple val(meta), path("${meta.id}_genes.csv"), path("${meta.id}_genome.png"), path("${meta.id}_genome2.png"), emit: stats

    script:
    """
    export HOME=/tmp

    '${params.nucdyn_launcher}' txstart_stats \
        --input '${txstart}' \
        --genome '${annotation}' \
        --out_genes '${meta.id}_genes.csv' \
        --out_gw '${meta.id}_genome.png' \
        --out_gw2 '${meta.id}_genome2.png' \
        2>&1 | tee '${meta.id}.txstart_stats.log'

    if grep -q 'ERROR while running' '${meta.id}.txstart_stats.log'; then
        exit 1
    fi

    test -s '${meta.id}_genes.csv'
    test -s '${meta.id}_genome.png'
    test -s '${meta.id}_genome2.png'
    """
}


process NUCDYN_COMPARE {
    tag "${meta.id}: ${meta.treatment} vs ${meta.control}"
    label 'process_nucdyn_compare'

    container "${params.nucdyn_container}"
    containerOptions '--entrypoint "" -u $(id -u):$(id -g)'

    publishDir "${params.outdir}/NucleosomeDynamics/NucDyn", mode: 'copy'

    input:
    tuple val(meta), path(treatment_rdata), path(control_rdata), path(treatment_calls), path(control_calls)
    path genome_sizes

    output:
    tuple val(meta), path("${meta.id}.gff"), path("${meta.id}.bw"), emit: dynamics

    script:
    """
    export HOME=/tmp

    '${params.nucdyn_launcher}' nucDyn \
        --input2 '${treatment_rdata}' \
        --input1 '${control_rdata}' \
        --calls2 '${treatment_calls}' \
        --calls1 '${control_calls}' \
        --outputGff '${meta.id}.gff' \
        --outputBigWig '${meta.id}.bw' \
        --genome '${genome_sizes}' \
        2>&1 | tee '${meta.id}.nucDyn.log'

    if grep -q 'ERROR while running' '${meta.id}.nucDyn.log'; then
        exit 1
    fi

    test -s '${meta.id}.gff'
    test -s '${meta.id}.bw'
    """
}


process NUCDYN_STATS {
    tag "${meta.id}"
    label 'process_nucdyn'

    container "${params.nucdyn_container}"
    containerOptions '--entrypoint "" -u $(id -u):$(id -g)'

    publishDir "${params.outdir}/NucleosomeDynamics/statistics/NucDyn", mode: 'copy'

    input:
    tuple val(meta), path(dynamics_gff), path(dynamics_bw)
    path annotation

    output:
    tuple val(meta), path("${meta.id}_genes.csv"), path("${meta.id}_genome.png"), emit: stats

    script:
    """
    export HOME=/tmp

    '${params.nucdyn_launcher}' nucDyn_stats \
        --input '${dynamics_gff}' \
        --genome '${annotation}' \
        --out_genes '${meta.id}_genes.csv' \
        --out_gw '${meta.id}_genome.png' \
        2>&1 | tee '${meta.id}.nucDyn_stats.log'

    if grep -q 'ERROR while running' '${meta.id}.nucDyn_stats.log'; then
        exit 1
    fi

    test -s '${meta.id}_genes.csv'
    test -s '${meta.id}_genome.png'
    """
}
