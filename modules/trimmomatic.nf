process TRIMMOMATIC {
    tag "${meta.id}"
    publishDir "${params.output_dir}/trimmomatic", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}_R{1,2}_trimmed.fastq.gz"), emit: trimmed_reads
    tuple val(meta), path("${prefix}.trimlog.txt"), emit: log

    script:
    prefix = "${meta.id}"
    """
    trimmomatic PE \
        -threads ${task.cpus} \
        -phred33 \
        -trimlog ${prefix}.trimlog.txt \
        ${reads[0]} ${reads[1]} \
        ${prefix}_R1_trimmed.fastq.gz ${prefix}_R1_unpaired.fastq.gz \
        ${prefix}_R2_trimmed.fastq.gz ${prefix}_R2_unpaired.fastq.gz \
        SLIDINGWINDOW:4:20 MINLEN:50 TRAILING:20
    """
}
