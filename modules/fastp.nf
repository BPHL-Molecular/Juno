process FASTP_REPORT {
    tag "${meta.id}"
    publishDir "${params.output_dir}/fastp", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.fastp.json"), emit: json
    path "${prefix}.fastp.html", emit: html

    script:
    prefix = "${meta.id}"
    """
    fastp \
        --in1 ${reads[0]} \
        --in2 ${reads[1]} \
        --out1 ${prefix}_R1_report_only.fastq.gz \
        --out2 ${prefix}_R2_report_only.fastq.gz \
        --json ${prefix}.fastp.json \
        --html ${prefix}.fastp.html \
        --disable_adapter_trimming \
        --disable_quality_filtering \
        --disable_length_filtering \
        --thread ${task.cpus}
    
    # Remove the temporary output files
    rm ${prefix}_R1_report_only.fastq.gz ${prefix}_R2_report_only.fastq.gz
    """
}
