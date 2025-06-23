process SPADES {
    tag "$meta.id"
    publishDir "${params.output_dir}/contigs", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(filtered_reads)

    output:
    tuple val(meta), path("${meta.id}_contigs.fasta"), emit: contigs
    tuple val(meta), path("${meta.id}_spades.log"), emit: log

    when:
    params.assembly_mode == 'denovo'

    script:
    prefix = "${meta.id}"
    """
    spades.py \\
        -1 ${filtered_reads[0]} \\
        -2 ${filtered_reads[1]} \\
        --threads ${task.cpus} \\
        --memory ${task.memory.toGiga()} \\
        --isolate \\
        --cov-cutoff auto \\
        -o ${prefix}

    mv ${prefix}/contigs.fasta ${prefix}_contigs.fasta
    mv ${prefix}/spades.log ${prefix}_spades.log
    """
}