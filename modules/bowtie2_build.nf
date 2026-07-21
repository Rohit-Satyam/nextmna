// Nextflow DSL2 process converted from Snakemake rule: index_bowtie2
process BOWTIE2_BUILD {
    tag "reference"
    label 'process_high'

    conda "${projectDir}/envs/bowtie2.yml"
    publishDir "${params.outdir}/00_index/bowtie2", mode: 'copy', overwrite: true, pattern: 'reference*.bt2'
    publishDir "${params.outdir}/00_logs/bowtie2_index", mode: 'copy', overwrite: true, pattern: '*.log'

    input:
    path reference

    output:
    path "reference*.bt2", emit: index
    path "reference.bowtie2-build.log", emit: log
    path "sizes.genome", emit: sizes

    script:
    """
    bowtie2-build \\
      --seed ${params.rnd_seed} \\
      --threads ${task.cpus} \\
      "${reference}" reference > reference.bowtie2-build.log 2>&1
      samtools faidx ${reference}
      cut -f1,2 ${reference}.fai > sizes.genome

    """
}
