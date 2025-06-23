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
