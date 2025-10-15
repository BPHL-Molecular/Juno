process KRAKEN2 {
    tag "${meta.id}"
    publishDir "${params.output_dir}/kraken2", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.kraken2.report"), emit: report
    tuple val(meta), path("${prefix}.kraken2.out"), emit: classifications

    script:
    prefix = "${meta.id}"
    """
    kraken2 \
        --db ${params.kraken2_db} \
        --threads ${task.cpus} \
        --paired \
        --report ${prefix}.kraken2.report \
        --output ${prefix}.kraken2.out \
        ${reads[0]} ${reads[1]}
    """
}
