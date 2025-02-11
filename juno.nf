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
        tuple(ref_id, ref)
    }
    .tap { refs_ch }
    .map { ref_id, ref -> ref_id }
    .collect()
    .view { refs ->
        "Found ${refs.size()} references:\n${refs.sort().join('\n')}"
    }

// Modules
include { HRRT }            from './modules/hrrt'
include { FASTP }           from './modules/fastp'
include { KRAKEN2 }         from './modules/kraken2'
include { MINIMAP2 }        from './modules/minimap2'
include { SAMTOOLS }        from './modules/samtools'
include { IVAR_VARIANTS }   from './modules/ivar'
include { IVAR_CONSENSUS }  from './modules/ivar'
include { QUAST }           from './modules/quast'
include { SUMMARY_REPORT }  from './modules/summary_report'
include { MULTIQC }         from './modules/multiqc'

workflow {
    // Initial QC and trimming
    HRRT(input_ch)
    FASTP(HRRT.out.dehosted_reads)
    KRAKEN2(FASTP.out.trimmed_reads)

    // Alignment and BAM processing chain
    FASTP.out.trimmed_reads
        .combine(refs_ch)
        .map { meta, trimmed_reads, ref_id, ref ->
            def new_meta = [
                id: "${meta.id}_${ref_id}",
                sample_id: meta.id,
                ref_id: ref_id
            ]
            tuple(new_meta, trimmed_reads, ref)
        }
        .set { reads_with_refs }

    MINIMAP2(reads_with_refs)
    SAMTOOLS(MINIMAP2.out.sam)

    // Prepare input for IVAR processes
    SAMTOOLS.out.indexed_bam
        .map { meta, sorted_bam, bai ->
            [meta.ref_id, meta, sorted_bam, bai]
        }
        .combine(refs_ch, by: 0)
        .map { ref_id, meta, sorted_bam, bai, ref ->
            tuple(meta, sorted_bam, bai, ref)
        }
        .set { ivar_input }

    IVAR_VARIANTS(ivar_input)
    IVAR_CONSENSUS(ivar_input)

    // Prepare input for QUAST process
    IVAR_CONSENSUS.out.consensus
        .map { meta, consensus ->
            [meta.ref_id, meta, consensus]
        }
        .combine(refs_ch, by: 0)
        .map { ref_id, meta, consensus, ref ->
            tuple(meta, consensus, ref)
        }
        .set { consensus_ref_ch }

    QUAST(consensus_ref_ch)

    // Branch channels for summary report
    FASTP.out.json
        .map { meta, json -> json }
        .collect()
        .set { fastp_files }

    SAMTOOLS.out.coverage
        .map { meta, coverage -> coverage }
        .collect()
        .set { coverage_files }

    QUAST.out.stats
        .map { meta, stats -> stats }
        .collect()
        .set { quast_files }

    IVAR_VARIANTS.out.variants
        .map { meta, variants -> variants }
        .collect()
        .set { variants_files }

    IVAR_CONSENSUS.out.consensus
        .map { meta, consensus -> consensus }
        .collect()
        .set { consensus_files }

    // Summary report with key metrics
    SUMMARY_REPORT(fastp_files, coverage_files, quast_files, variants_files, consensus_files)

    // MultiQC report
    MULTIQC(SUMMARY_REPORT.out.summary)
}

