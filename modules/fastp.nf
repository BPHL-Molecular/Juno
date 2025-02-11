process FASTP {
    tag "${meta.id}"
    publishDir "${params.output_dir}/trimmed", mode: 'copy'

    input:
    tuple val(meta), path(dehosted_reads)

    output:
    tuple val(meta), path("${prefix}_R{1,2}_trimmed.fastq.gz"), emit: trimmed_reads
    tuple val(meta), path("${prefix}.fastp.json"), emit: json
    path "${prefix}.fastp.html", emit: html

    script:
    prefix = "${meta.id}"
    """
    fastp \
        --in1 ${dehosted_reads[0]} \
        --in2 ${dehosted_reads[1]} \
        --out1 ${prefix}_R1_trimmed.fastq.gz \
        --out2 ${prefix}_R2_trimmed.fastq.gz \
        --json ${prefix}.fastp.json \
        --html ${prefix}.fastp.html \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        --thread ${task.cpus}
    """
}
