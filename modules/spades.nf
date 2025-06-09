process SPADES {
    tag "$meta.id"
    publishDir "${params.output_dir}/${meta.id}/assembly", mode: 'copy'

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
        --only-assembler \\
        --careful \\
        -o spades_output

    # Copy main outputs with proper naming
    cp spades_output/contigs.fasta ${prefix}_contigs.fasta
    cp spades_output/spades.log ${prefix}_spades.log
    """
}