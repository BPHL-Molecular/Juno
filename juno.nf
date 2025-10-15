#!/usr/bin/env nextflow

/*
  Juno Pipeline (named after Juno Beach in Florida)
  Florida's BPHL Nextflow pipeline for OROV genome assembly
  Author: Arnold Rodriguez Hilario
  Email: arnold.rodriguezhilario@flhealth.gov
*/

nextflow.enable.dsl = 2

// Print pipeline info
log.info """
    Juno - Florida's BPHL Nextflow pipeline for OROV genome assembly
    ==========================================================================
    input dir      : ${params.input_dir}
    output dir     : ${params.output_dir}
    kraken2 db     : ${params.kraken2_db}
    assembly mode  : ${params.assembly_mode}
    skip hrrt      : ${params.skip_hrrt}
    polish contigs : ${params.polish_contigs}
    ==========================================================================
    """

// Input channel
Channel
    .fromFilePairs("${params.input_dir}/*_L001_R{1,2}_*.fastq.gz", checkIfExists: true)
    .map { id, files ->
        def meta = [
            id: id.split('_L001')[0],
            single_end: false
        ]
        [ meta, files ]
    }
    .tap { input_ch }
    .map { meta, files -> meta.id }
    .collect()
    .view { samples ->
        "Found ${samples.size()} samples:\n${samples.sort().join('\n')}"
    }

// Reference channel
ref_files_ch = Channel
    .fromPath("${projectDir}/references/*.{fasta,fa}")
    .ifEmpty { error "No reference files found in ${projectDir}/references" }
    .collect()

Channel
    .fromPath("${projectDir}/references/*.{fasta,fa}")
    .ifEmpty { error "No reference files found in ${projectDir}/references" }
    .map { ref ->
        def ref_id = ref.name.replaceAll(/\.(fasta|fa)$/, '')
        def gff = file("${projectDir}/references/${ref_id}.gff3", checkIfExists: true)
        [ref_id, ref, gff]
    }
    .multiMap { ref_id, ref, gff ->
        basic: tuple(ref_id, ref)
        gff: tuple(ref_id, ref, gff)
        files: ref
        ids: ref_id
    }
    .set { refs }

// Reference segment mapping
segment_ref_mapping = Channel.of(
    ['L', 'PQ064919.1'],
    ['M', 'PQ064920.1'],
    ['S', 'PQ064921.1']
)

// Log references
refs.ids
    .collect()
    .view { ids ->
        "Using ${ids.size()} references:\n${ids.sort().join('\n')}"
    }

// Modules
include { FASTQ_SCAN }               from './modules/fastq_scan'
include { TRIMMOMATIC }              from './modules/trimmomatic'
include { BBDUK_ADAPTERS }           from './modules/bbduk'
include { BBDUK_PHIX }               from './modules/bbduk'
include { HRRT }                     from './modules/hrrt'
include { FASTP_REPORT }             from './modules/fastp'
include { KRAKEN2 }                  from './modules/kraken2'
include { KRAKENTOOLS }              from './modules/krakentools'
include { BWA }                      from './modules/bwa'
include { BWA_VALIDATE }             from './modules/bwa'
include { SAMTOOLS }                 from './modules/samtools'
include { SAMTOOLS_DENOVO }          from './modules/samtools'
include { POLISH_CONTIGS }           from './modules/polish_contigs'
include { TRIM_TERMINALS }           from './modules/trim_terminals'
include { IVAR_VARIANTS }            from './modules/ivar'
include { IVAR_CONSENSUS }           from './modules/ivar'
include { SPADES }                   from './modules/spades'
include { MAKEBLASTDB }              from './modules/blast'
include { BLAST }                    from './modules/blast'
include { CLASSIFY_CONTIGS }         from './modules/classify_contigs'
include { QUAST }                    from './modules/quast'
include { SUMMARY_REPORT_REFERENCE } from './modules/summary_report'
include { SUMMARY_REPORT_DENOVO }    from './modules/summary_report'
include { MULTIQC }                  from './modules/multiqc'

workflow {
    // Raw read statistics with fastq-scan
    ch_fastq_scan = FASTQ_SCAN(input_ch)

    // Read QC Processing
    ch_trimmomatic = input_ch | TRIMMOMATIC
    ch_bbduk_adapters = BBDUK_ADAPTERS(ch_trimmomatic.trimmed_reads)
    ch_bbduk_phix = BBDUK_PHIX(ch_bbduk_adapters.adaptertrim_reads)

    // Conditional human read removal (HRRT)
    if (params.skip_hrrt) {
        // Skip HRRT, use reads from BBDUK
        ch_clean_reads = ch_bbduk_phix.clean_reads
    } else {
        // Apply HRRT to remove human reads
        ch_clean_reads = HRRT(ch_bbduk_phix.clean_reads).dehosted_reads
    }

    // Generate FASTP report on final clean reads
    ch_processed = FASTP_REPORT(ch_clean_reads)

    // Taxonomic classification
    ch_kraken = KRAKEN2(ch_clean_reads)

    // Prepare channels for kraken outputs
    ch_kraken_fixed = ch_kraken.report
        .map { meta, report ->
            def classifications_file = file(report.toString().replace('.kraken2.report', '.kraken2.out'))
            tuple(meta, classifications_file)
        }

    // Update the channel for extraction of OROV reads
    ch_for_extraction = ch_clean_reads
        .join(ch_kraken.report, by: 0)
        .join(ch_kraken_fixed, by: 0)
        .map { meta, clean_reads, kraken_report, kraken_out ->
            tuple(meta, clean_reads, kraken_report, kraken_out)
        }

    // Filtered OROV reads channel
    ch_filtered_reads = KRAKENTOOLS(ch_for_extraction).filtered_reads

    // Conditional assembly workflow based on assembly mode
    if (params.assembly_mode == 'denovo') {
        // DE NOVO ASSEMBLY WORKFLOW

        // BLAST database channel
        ch_blast_db = MAKEBLASTDB(ref_files_ch)

        // De novo assembly channel
        ch_contigs = SPADES(ch_filtered_reads)

        // Validation of assembly through read mapping
        ch_validate = ch_contigs.contigs
            .join(ch_filtered_reads, by: 0)
            .map { meta, contigs, filtered_reads ->
                [meta, contigs, filtered_reads]
            } \
            | BWA_VALIDATE

        // Get validation statistics
        ch_validation_stats = SAMTOOLS_DENOVO(ch_validate.sam)

        // Polish contigs if enabled
        if (params.polish_contigs) {
            ch_polished = ch_contigs.contigs
                .join(ch_validation_stats.sorted_bam, by: 0)
                .join(ch_validation_stats.indexed_bam, by: 0)
                .map { meta, contigs, sorted_bam, indexed_bam, bai ->
                    tuple(meta, contigs, sorted_bam, bai)
                } \
                | POLISH_CONTIGS

            // Trim low-coverage terminal regions
            ch_trimmed = ch_polished.polished_fasta
                .join(ch_polished.vcf, by: 0) \
                | TRIM_TERMINALS

            // Use trimmed contigs for downstream analysis
            ch_contigs_for_analysis = ch_trimmed.trimmed_fasta
        } else {
            // Use original contigs for downstream analysis
            ch_contigs_for_analysis = ch_contigs.contigs
        }

        // BLAST classification channel
        ch_blast = BLAST(
            ch_contigs_for_analysis,
            ch_blast_db.db_files
        )

        // Contig classification channel
        ch_classified = ch_contigs_for_analysis
            .join(ch_blast.blast_results, by: 0)
            .map { meta, contigs, blast_results ->
                tuple(meta, contigs, blast_results)
            } \
            | CLASSIFY_CONTIGS

        // Prepare segments for QUAST
        ch_segments_for_quast = ch_classified.segment_L
            .map { meta, segment_file ->
                tuple(meta, segment_file, 'L')
            }
            .mix(
                ch_classified.segment_M.map { meta, segment_file ->
                    tuple(meta, segment_file, 'M')
                },
                ch_classified.segment_S.map { meta, segment_file ->
                    tuple(meta, segment_file, 'S')
                }
            )
            .filter { meta, segment_file, segment_type ->
                segment_file.size() > 50
            }
            .map { meta, segment_file, segment_type ->
                [segment_type, meta, segment_file]
            }
            .combine(segment_ref_mapping, by: 0)
            .map { segment_type, meta, segment_file, ref_id ->
                [ref_id, meta, segment_file, segment_type]
            }
            .combine(refs.basic, by: 0)
            .map { ref_id, meta, segment_file, segment_type, ref_file ->
                def new_meta = [
                    id: "${meta.id}_${segment_type}",
                    sample_id: meta.id,
                    ref_id: segment_type
                ]
                tuple(new_meta, segment_file, ref_file)
            }

        // Assembly evaluation channel
        ch_quast = ch_segments_for_quast | QUAST

    } else {
        // REFERENCE-BASED ASSEMBLY WORKFLOW

        // Alignment and SAM/BAM channels
        ch_aligned = ch_filtered_reads
            .combine(refs.basic)
            .map { meta, filtered_reads, ref_id, ref ->
                def new_meta = [
                    id: "${meta.id}_${ref_id}",
                    sample_id: meta.id,
                    ref_id: ref_id
                ]
                tuple(new_meta, filtered_reads, ref)
            } \
            | BWA \
            | SAMTOOLS

        // Variant calling channel
        ch_variants = ch_aligned.indexed_bam
            .map { meta, dedup_bam, bai ->
                [meta.ref_id, meta, dedup_bam, bai]
            }
            .combine(refs.gff, by: 0)
            .map { ref_id, meta, dedup_bam, bai, ref, gff ->
                tuple(meta, dedup_bam, bai, ref, gff)
            } \
            | IVAR_VARIANTS

        // Consensus generation and post-assembly QC channels
        ch_consensus = ch_aligned.indexed_bam
            .map { meta, dedup_bam, bai ->
                [meta.ref_id, meta, dedup_bam, bai]
            }
            .combine(refs.basic, by: 0)
            .map { ref_id, meta, dedup_bam, bai, ref ->
                tuple(meta, dedup_bam, bai, ref)
            } \
            | IVAR_CONSENSUS

        ch_quast = ch_consensus.consensus
            .map { meta, cons ->
                [meta.ref_id, meta, cons]
            }
            .combine(refs.basic, by: 0)
            .map { ref_id, meta, cons, ref ->
                tuple(meta, cons, ref)
            } \
            | QUAST
    }

// Summary report channels
    ch_fastq_scan_out = ch_fastq_scan.json.map { meta, json -> json }.collect()
    ch_fastp = ch_processed.json.map { meta, json -> json }.collect()
    ch_kraken2_out = ch_kraken.report.map { meta, report -> report }.collect()

    if (params.assembly_mode == 'denovo') {
        // De novo summary report
        ch_quast_stats = ch_quast.stats.map { meta, stats -> stats }.collect()
        ch_classification_out = ch_classified.classification_summary.map { meta, summary -> summary }.collect()
        ch_validation_coverage = ch_validation_stats.coverage.map { meta, coverage -> coverage }.collect()
        ch_validation_flagstat = ch_validation_stats.flagstat.map { meta, flagstat -> flagstat }.collect()

        SUMMARY_REPORT_DENOVO(
            ch_fastq_scan_out,
            ch_fastp,
            ch_kraken2_out,
            ch_quast_stats,
            ch_classification_out,
            ch_validation_coverage,
            ch_validation_flagstat
        )

        SUMMARY_REPORT_DENOVO.out.summary | MULTIQC

    } else {
        // Reference-based summary report
        ch_coverage = ch_aligned.coverage.map { meta, coverage -> coverage }.collect()
        ch_quast_stats = ch_quast.stats.map { meta, stats -> stats }.collect()
        ch_variants_out = ch_variants.variants.map { meta, variants -> variants }.collect()
        ch_consensus_out = ch_consensus.consensus.map { meta, consensus -> consensus }.collect()
        ch_markdup_stats = ch_aligned.markdup_stats.map { meta, stats -> stats }.collect()

        SUMMARY_REPORT_REFERENCE(
            ch_fastq_scan_out,
            ch_fastp,
            ch_kraken2_out,
            ch_coverage,
            ch_quast_stats,
            ch_variants_out,
            ch_consensus_out,
            ch_markdup_stats
        )

        SUMMARY_REPORT_REFERENCE.out.summary | MULTIQC
    }
}
