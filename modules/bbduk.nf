process BBDUK_ADAPTERS {
    tag "${meta.id}"
    publishDir "${params.output_dir}/bbduk/bbduk_adapters", mode: 'copy'

    input:
    tuple val(meta), path(trimmed_reads)

    output:
    tuple val(meta), path("${prefix}_R{1,2}_adaptertrim.fastq.gz"), emit: adaptertrim_reads
    tuple val(meta), path("${prefix}_adapter_stats.txt"), emit: adapter_stats

    script:
    prefix = "${meta.id}"
    """
    bbduk.sh \
        in1=${trimmed_reads[0]} \
        in2=${trimmed_reads[1]} \
        out1=${prefix}_R1_adaptertrim.fastq.gz \
        out2=${prefix}_R2_adaptertrim.fastq.gz \
        ref=/bbmap/resources/adapters.fa \
        ktrim=r k=23 mink=11 hdist=1 tpe tbo \
        stats=${prefix}_adapter_stats.txt
    """
}

process BBDUK_PHIX {
    tag "${meta.id}"
    publishDir "${params.output_dir}/bbduk/bbduk_phix", mode: 'copy'

    input:
    tuple val(meta), path(adaptertrim_reads)

    output:
    tuple val(meta), path("${prefix}_R{1,2}_clean.fastq.gz"), emit: clean_reads
    tuple val(meta), path("${prefix}_phix_stats.txt"), emit: phix_stats

    script:
    prefix = "${meta.id}"
    """
    bbduk.sh \
        in1=${adaptertrim_reads[0]} \
        in2=${adaptertrim_reads[1]} \
        out1=${prefix}_R1_clean.fastq.gz \
        out2=${prefix}_R2_clean.fastq.gz \
        ref=/bbmap/resources/phix174_ill.ref.fa.gz \
        k=31 hdist=1 \
        stats=${prefix}_phix_stats.txt
    """
}
