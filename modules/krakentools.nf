process KRAKENTOOLS {
    tag "${meta.id}"
    publishDir "${params.output_dir}/krakentools", mode: 'copy'

    input:
    tuple val(meta), path(reads), path(kraken_report), path(kraken_out)

    output:
    tuple val(meta), path("${prefix}_R{1,2}_filtered.fastq"), emit: filtered_reads

    script:
    prefix = "${meta.id}"
    """
    python3 /KrakenTools/extract_kraken_reads.py \
        -k ${kraken_out} \
        -s1 ${reads[0]} \
        -s2 ${reads[1]} \
        --fastq-output \
        -r ${kraken_report} \
        -t 3052429 2571170 118655 \
        --include-children \
        -o ${prefix}_R1_filtered.fastq \
        -o2 ${prefix}_R2_filtered.fastq
    """
}
