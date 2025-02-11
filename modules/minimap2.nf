process MINIMAP2 {
    tag "${meta.id}"
    publishDir "${params.output_dir}/alignments", mode: 'copy'

    input:
    tuple val(meta), path(trimmed_reads), path(reference)

    output:
    tuple val(meta), path("${prefix}.sam"), emit: sam

    script:
    prefix = "${meta.id}"
    """
    minimap2 \
        -ax sr \
        -t ${task.cpus} \
        ${reference} \
        ${trimmed_reads[0]} ${trimmed_reads[1]} \
        > ${prefix}.sam
    """
}
