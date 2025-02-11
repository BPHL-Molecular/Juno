process IVAR_VARIANTS {
    tag "${meta.id}"
    publishDir "${params.output_dir}/variants", mode: 'copy'

    input:
    tuple val(meta), path(sorted_bam), path(bai), path(reference)

    output:
    tuple val(meta), path("${prefix}.variants.tsv"), emit: variants

    script:
    prefix = "${meta.id}"
    """
    samtools mpileup \
        -aa -A -d 0 -B -Q 0 \
        -f ${reference} \
        ${sorted_bam} | \
    ivar variants \
        -p ${prefix} \
        -r ${reference} \
        -m 10 \
        -q 20 \
        -t 0.6

    mv ${prefix}.tsv ${prefix}.variants.tsv
    """
}

process IVAR_CONSENSUS {
    tag "${meta.id}"
    publishDir "${params.output_dir}/consensus", mode: 'copy'

    input:
    tuple val(meta), path(sorted_bam), path(bai), path(reference)

    output:
    tuple val(meta), path("${prefix}.consensus.fasta"), emit: consensus
    tuple val(meta), path("${prefix}.consensus.fasta"), val(meta.ref_id), emit: consensus_with_ref

    script:
    prefix = "${meta.id}"
    """
    samtools mpileup \
        -aa -A -d 0 -B -Q 0 \
        -f ${reference} \
        ${sorted_bam} | \
    ivar consensus \
        -p ${prefix} \
        -q 20 \
        -t 0.6 \
        -m 10 \
        -n N

    # Format the consensus FASTA with proper line wrapping
    if [ -f "${prefix}.fa" ]; then
        echo ">${prefix}" > "${prefix}.consensus.fasta"
        grep -v ">" "${prefix}.fa" | fold -w 70 >> "${prefix}.consensus.fasta"
    fi
    """
}
