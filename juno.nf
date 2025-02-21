#!/usr/bin/env nextflow

/*
  Juno Pipeline
  Florida's BPHL Nextflow pipeline for OROV reference-based assembly from metagenomics reads.
  Author: Arnold Rodriguez Hilario
  Email: arnold.rodriguezhilario@flhealth.gov
*/

nextflow.enable.dsl = 2

// Print pipeline info
log.info """
    Juno - Florida's BPHL Nextflow pipeline for OROV reference-based assembly from metagenomics reads.
    ==========================================================================
    input dir    : ${params.input_dir}
    output dir   : ${params.output_dir}
    references   : ${params.refs_dir}
    kraken2 db   : ${params.kraken2_db}
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
Channel
    .fromPath("${params.refs_dir}/*.{fasta,fa}")
    .ifEmpty { error "No reference files found in ${params.refs_dir}" }
    .map { ref ->
        def ref_id = ref.name.replaceAll(/\.(fasta|fa)$/, '')
        def gff = file("${params.refs_dir}/${ref_id}.gff3", checkIfExists: true)
        [ref_id, ref, gff]
    }
    .multiMap { ref_id, ref, gff ->
        basic: tuple(ref_id, ref)
        gff: tuple(ref_id, ref, gff)
        ids: ref_id
    }
    .set { refs }

// Log found references
refs.ids
    .collect()
    .view { ids ->
        "Found ${ids.size()} references:\n${ids.sort().join('\n')}"
    }

// Modules
include { HRRT }            from './modules/hrrt'
include { FASTP }           from './modules/fastp'
include { KRAKEN2 }         from './modules/kraken2'
include { BWA }             from './modules/bwa'
include { SAMTOOLS }        from './modules/samtools'
include { IVAR_VARIANTS }   from './modules/ivar'
include { IVAR_CONSENSUS }  from './modules/ivar'
include { QUAST }           from './modules/quast'
include { SUMMARY_REPORT }  from './modules/summary_report'
include { MULTIQC }         from './modules/multiqc'

workflow {
    // Initial QC
    ch_fastq = HRRT(input_ch) | FASTP

    // Taxonomic classification
    KRAKEN2(FASTP.out.trimmed_reads)

    // Alignment and SAM/BAM processing
    ch_aligned = FASTP.out.trimmed_reads \
        .combine(refs.basic) \
        .map { meta, trimmed_reads, ref_id, ref ->
            def new_meta = [
                id: "${meta.id}_${ref_id}",
                sample_id: meta.id,
                ref_id: ref_id
            ]
            tuple(new_meta, trimmed_reads, ref)
        } \
        | BWA \
        | SAMTOOLS

    // Variant calling
    ch_variants = ch_aligned.indexed_bam
        .map { meta, sorted_bam, bai ->
            [meta.ref_id, meta, sorted_bam, bai]
        }
        .combine(refs.gff, by: 0)
        .map { ref_id, meta, sorted_bam, bai, ref, gff ->
            tuple(meta, sorted_bam, bai, ref, gff)
        } \
        | IVAR_VARIANTS

    // Consensus generation and post-assembly QC
    ch_consensus = ch_aligned.indexed_bam
        .map { meta, sorted_bam, bai ->
            [meta.ref_id, meta, sorted_bam, bai]
        }
        .combine(refs.basic, by: 0)
        .map { ref_id, meta, sorted_bam, bai, ref ->
            tuple(meta, sorted_bam, bai, ref)
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

    // Collect outputs for reporting
    ch_fastp = FASTP.out.json.map { meta, json -> json }.collect()
    ch_kraken2_out = KRAKEN2.out.report.map { meta, report -> report }.collect()
    ch_coverage = ch_aligned.coverage.map { meta, coverage -> coverage }.collect()
    ch_quast_stats = ch_quast.stats.map { meta, stats -> stats }.collect()
    ch_variants_out = ch_variants.variants.map { meta, variants -> variants }.collect()
    ch_consensus_out = ch_consensus.consensus.map { meta, consensus -> consensus }.collect()

    // Generate summary report
    SUMMARY_REPORT(
        ch_fastp,
        ch_kraken2_out,
        ch_coverage,
        ch_quast_stats,
        ch_variants_out,
        ch_consensus_out
    ) | MULTIQC
}
