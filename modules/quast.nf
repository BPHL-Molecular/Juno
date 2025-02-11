process QUAST {
    tag "${meta.id}"
    publishDir "${params.output_dir}/quast", mode: 'copy'

    input:
    tuple val(meta), path(consensus), path(reference)

    output:
    tuple val(meta), path("${prefix}"), emit: report
    tuple val(meta), path("${prefix}.report.tsv"), emit: stats
    tuple val(meta), path("${prefix}/transposed_report.tsv"), emit: transposed_report


    script:
    prefix = "${meta.id}_quast"
    """
    quast.py \
        ${consensus} \
        -r ${reference} \
        -o ${prefix} \
        --threads ${task.cpus} \
        --min-contig 500 \
        --no-plots \
        --no-html \
        --no-icarus

    mv ${prefix}/report.tsv ${prefix}.report.tsv
    """
}
