process BWA {
    tag "${meta.id}"
    publishDir "${params.output_dir}/alignments", mode: 'copy'

    input:
    tuple val(meta), path(filtered_reads), path(reference)

    output:
    tuple val(meta), path("${prefix}.sam"), emit: sam

    script:
    prefix = "${meta.id}"
    """
    # Index references
    bwa index ${reference}

    # Align reads
    bwa mem \
        -t ${task.cpus} \
        ${reference} \
        ${filtered_reads[0]} ${filtered_reads[1]} \
        > ${prefix}.sam
    """
}

process BWA_VALIDATE {
    tag "${meta.id}"
    publishDir "${params.output_dir}/samtools/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(contigs), path(filtered_reads)

    output:
    tuple val(meta), path("${prefix}_validation.sam"), emit: sam

    when:
    params.assembly_mode == 'denovo'

    script:
    prefix = "${meta.id}"
    """
    # Index contigs
    bwa index ${contigs}

    # Align reads back to contigs
    bwa mem \\
        -t ${task.cpus} \\
        ${contigs} \\
        ${filtered_reads[0]} ${filtered_reads[1]} \\
        > ${prefix}_validation.sam
    """
}
