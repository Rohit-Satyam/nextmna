include { FASTQC_RAW }             from './modules/fastqc_raw'
include { TRIM_READS }             from './modules/trim_reads'
include { BOWTIE2_BUILD }          from './modules/bowtie2_build'
include { BOWTIE2_ALIGN }          from './modules/bowtie2_align'
include { SAM_TO_BAM }             from './modules/sam_to_bam'
include { SORT_BAM }               from './modules/sort_bam'
include { INDEX_BAM }              from './modules/index_bam'
include { QUALIMAP_MAPPING }       from './modules/qualimap_mapping'
include { MULTIQC_MAPPING }        from './modules/multiqc_mapping'
include { ALIGNMENT_SIEVE }        from './modules/alignment_sieve'
include { SORT_FILTERED_BAM }      from './modules/sort_filtered_bam'
include { INDEX_FILTERED_BAM }     from './modules/index_filtered_bam'
include { QUALIMAP_FILTERED }      from './modules/qualimap_filtered'
include { MULTIQC_FILTERED }       from './modules/multiqc_filtered'
include { PREPARE_DANPOS_INPUT }   from './modules/prepare_danpos_input'
include { DANPOS_DIFF }            from './modules/danpos_diff'
include { ADJUST_FDR }             from './modules/adjust_fdr'
include { MERGE_SORT_NUCDYN_BAM; NUCDYN_READ_BAM; NUCDYN_NUCLER; NUCDYN_NUCLER_STATS; NUCDYN_TXSTART; NUCDYN_TXSTART_STATS; NUCDYN_COMPARE; NUCDYN_STATS } from './modules/nucleosome_dynamics'

workflow {
    // Validate required inputs before creating any channels or submitting tasks.
    if( !params.input ) {
        error "Missing required parameter: --input <samples.csv|samples.tsv>"
    }
    if( !params.reference ) {
        error "Missing required parameter: --reference <genome.fasta>"
    }
    if( !params.annotation ) {
        error "Missing required parameter: --annotation <annotation.gff>"
    }
    if( !params.nucdyn_annotation ) {
        error "Missing required parameter: --nucdyn_annotation <pf_nucDyn.gff>"
    }
    if( !params.nucdyn_container ) {
        error "Missing required parameter: --nucdyn_container <Docker image>"
    }


    // Parse CSV/TSV metadata and validate the ordered DANPOS contrast.
    def parsed = SamplesheetParser.parse(params.input, params.danpos_comparison)
    def treatment = parsed.treatment
    def control   = parsed.control

    // NucDyn uses input2/calls2 as treatment and input1/calls1 as control.
    def nucdyn_comparison = (params.nucdyn_comparison ?: params.danpos_comparison) as String
    def nucdyn_parts = nucdyn_comparison.tokenize(':')
    if( nucdyn_parts.size() != 2 || nucdyn_parts[0] == nucdyn_parts[1] ) {
        error "--nucdyn_comparison must contain two distinct conditions: treatment:control"
    }
    def nucdyn_treatment = nucdyn_parts[0]
    def nucdyn_control   = nucdyn_parts[1]
    def missing_nucdyn_conditions = [nucdyn_treatment, nucdyn_control].findAll { condition ->
        !parsed.conditions.contains(condition)
    }
    if( missing_nucdyn_conditions ) {
        error "NucDyn conditions not found in the samplesheet: ${missing_nucdyn_conditions.join(', ')}. Available conditions: ${parsed.conditions.join(', ')}"
    }

    log.info "Samplesheet: ${new java.io.File(params.input.toString()).absolutePath}"
    log.info "Samples: ${parsed.records.size()}"
    log.info "Conditions: ${parsed.conditions.join(', ')}"
    log.info "DANPOS contrast: ${treatment} versus ${control}"
    log.info "DANPOS path expression: ${treatment}/:${control}/"
    log.info "NucDyn contrast: ${nucdyn_treatment} versus ${nucdyn_control}"
    log.info "NucDyn result prefix: ${params.nucdyn_contrast_name}"

    // One metadata/read tuple per sample. DSL2 broadcasts this channel to FastQC and trimming.
    def samples_ch = channel
        .fromList(parsed.records)
        .map { rec ->
            def meta = [
                id        : rec.id,
                condition : rec.condition,
                replicate : rec.replicate
            ]
            tuple(
                meta,
                file(rec.r1, checkIfExists: true),
                file(rec.r2, checkIfExists: true)
            )
        }

    def reference_ch  = channel.value(file(params.reference, checkIfExists: true))
    def annotation_ch = channel.value(file(params.annotation, checkIfExists: true))
    def nucdyn_annotation_ch = channel.value(file(params.nucdyn_annotation, checkIfExists: true))
    def danpos_script = file("${projectDir}/bin/DANPOS3/danpos.py")
    if( !danpos_script.exists() ) {
        error "DANPOS3 is not installed at bin/DANPOS3/danpos.py. Run ./scripts/install_danpos3.sh or place DANPOS3 there."
    }
    def danpos_code_ch = channel.value(file("${projectDir}/bin/DANPOS3", checkIfExists: true))

    // Raw-read QC and preprocessing.
    FASTQC_RAW(samples_ch)
    TRIM_READS(samples_ch)

    // Reference indexing, alignment, BAM processing, and mapping QC.
    BOWTIE2_BUILD(reference_ch)
    def bowtie2_index_ch = BOWTIE2_BUILD.out.index

    BOWTIE2_ALIGN(
        TRIM_READS.out.reads,
        bowtie2_index_ch,
        params.rnd_seed as Integer,
        params.bowtie2_extra as String
    )
    SAM_TO_BAM(BOWTIE2_ALIGN.out.sam)
    SORT_BAM(SAM_TO_BAM.out.bam)
    INDEX_BAM(SORT_BAM.out.sorted_bam)

    QUALIMAP_MAPPING(INDEX_BAM.out.bam_bai)
    def mapping_qc_dirs_ch = QUALIMAP_MAPPING.out.qc_dir
        .map { meta, qc_dir -> qc_dir }
        .collect()
    MULTIQC_MAPPING(mapping_qc_dirs_ch)

    // Retain properly paired mononucleosome-sized fragments, then repeat QC.
    ALIGNMENT_SIEVE(
        INDEX_BAM.out.bam_bai,
        params.min_fragment_length as Integer,
        params.max_fragment_length as Integer,
        params.min_mapping_quality as Integer,
        params.sam_flag_include as Integer
    )
    SORT_FILTERED_BAM(ALIGNMENT_SIEVE.out.filtered_bam)
    INDEX_FILTERED_BAM(SORT_FILTERED_BAM.out.sorted_bam)

    QUALIMAP_FILTERED(INDEX_FILTERED_BAM.out.bam_bai, annotation_ch)
    def filtered_qc_dirs_ch = QUALIMAP_FILTERED.out.qc_dir
        .map { meta, qc_dir -> qc_dir }
        .collect()
    MULTIQC_FILTERED(filtered_qc_dirs_ch)

    // Restrict DANPOS to the two ordered conditions requested by the user.
    def selected_bams_ch = INDEX_FILTERED_BAM.out.bam_bai
        .filter { meta, bam, bai -> meta.condition in [treatment, control] }

    PREPARE_DANPOS_INPUT(selected_bams_ch)
    def danpos_input_dirs_ch = PREPARE_DANPOS_INPUT.out.staged
        .map { meta, staged_dir -> staged_dir }
        .collect()

    DANPOS_DIFF(
        danpos_input_dirs_ch,
        danpos_code_ch,
        treatment,
        control,
        params.danpos_extra as String
    )

    // Recalculate occupancy/fuzziness FDR using the authors' R implementation.
    ADJUST_FDR(
        DANPOS_DIFF.out.results,
        BOWTIE2_BUILD.out.sizes,
        treatment,
        control,
        params.adjust_fdr_n_simulations?.toString() ?: '',
        params.adjust_fdr_sampling_occupancy as String,
        params.adjust_fdr_sampling_fuzziness as String
    )


    // Merge biological replicates by condition and run the Nucleosome Dynamics suite.
    // The filtered BAM channel is used directly; no absolute alignment paths are needed.
    def nucdyn_grouped_bams_ch = INDEX_FILTERED_BAM.out.bam_bai
        .filter { meta, bam, bai -> meta.condition in [nucdyn_treatment, nucdyn_control] }
        .map { meta, bam, bai -> tuple(meta.condition as String, bam) }
        .groupTuple()
        .map { condition, bams ->
            def expected_replicates = params.nucdyn_expected_replicates as Integer
            if( expected_replicates > 0 && bams.size() != expected_replicates ) {
                error "NucDyn condition '${condition}' has ${bams.size()} BAM file(s); expected ${expected_replicates}"
            }

            def meta = [
                id        : condition,
                condition : condition
            ]
            tuple(meta, bams)
        }

    MERGE_SORT_NUCDYN_BAM(nucdyn_grouped_bams_ch)
    NUCDYN_READ_BAM(MERGE_SORT_NUCDYN_BAM.out.sorted_bam)
    NUCDYN_NUCLER(
        NUCDYN_READ_BAM.out.rdata,
        params.nucdyn_fragment_length as Integer,
        params.nucdyn_threshold_percentage as Integer
    )

    // Per-condition statistics and TxStart classification can run in parallel.
    NUCDYN_NUCLER_STATS(NUCDYN_NUCLER.out.calls, nucdyn_annotation_ch)
    NUCDYN_TXSTART(NUCDYN_NUCLER.out.calls, nucdyn_annotation_ch)
    NUCDYN_TXSTART_STATS(NUCDYN_TXSTART.out.txstart, nucdyn_annotation_ch)

    // Match each RData object to its nucleR calls by condition, then build one
    // deterministic treatment/control tuple for the final comparison.
    def nucdyn_condition_data_ch = NUCDYN_READ_BAM.out.rdata
        .map { meta, rdata -> tuple(meta.condition, rdata) }
        .join(
            NUCDYN_NUCLER.out.calls.map { meta, calls -> tuple(meta.condition, calls) },
            failOnDuplicate: true,
            failOnMismatch: true
        )

    def nucdyn_comparison_inputs_ch = nucdyn_condition_data_ch
        .toList()
        .map { rows ->
            def by_condition = [:]
            rows.each { row -> by_condition[row[0]] = row }

            if( !by_condition.containsKey(nucdyn_treatment) || !by_condition.containsKey(nucdyn_control) ) {
                error "NucDyn expected conditions '${nucdyn_treatment}' and '${nucdyn_control}', but found: ${by_condition.keySet().join(', ')}"
            }

            def treatment_row = by_condition[nucdyn_treatment]
            def control_row   = by_condition[nucdyn_control]
            def contrast_meta = [
                id        : params.nucdyn_contrast_name as String,
                treatment : nucdyn_treatment,
                control   : nucdyn_control
            ]

            tuple(
                contrast_meta,
                treatment_row[1],
                control_row[1],
                treatment_row[2],
                control_row[2]
            )
        }

    // BOWTIE2_BUILD already produces the genome-sizes file required by nucDyn.
    NUCDYN_COMPARE(nucdyn_comparison_inputs_ch, BOWTIE2_BUILD.out.sizes)
    NUCDYN_STATS(NUCDYN_COMPARE.out.dynamics, nucdyn_annotation_ch)
}
