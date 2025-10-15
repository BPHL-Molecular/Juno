process FASTQ_SCAN {
    tag "${meta.id}"
    publishDir "${params.output_dir}/fastq_scan", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.fastq-scan.json"), emit: json

    script:
    prefix = "${meta.id}"
    """
    gzip -dc ${reads[0]} | fastq-scan > ${prefix}.fastq-scan.json
    """
}
