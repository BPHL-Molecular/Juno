process IVAR_VARIANTS {
    tag "${meta.id}"
    publishDir "${params.output_dir}/variants", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(dedup_bam), path(bai), path(reference), path(gff)

    output:
    tuple val(meta), path("${prefix}.variants.tsv"), optional: true, emit: variants

    script:
    prefix = "${meta.id}"
    """
    samtools mpileup \
        -aa -A -d 100 -B -Q 0 \
        -f ${reference} \
        ${dedup_bam} | \
    ivar variants \
        -p ${prefix} \
        -r ${reference} \
        -g ${gff} \
        -m 10 \
        -q 20 \

    mv ${prefix}.tsv ${prefix}.variants.tsv
    """
}

process IVAR_CONSENSUS {
    tag "${meta.id}"
    publishDir "${params.output_dir}/consensus", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(dedup_bam), path(bai), path(reference)

    output:
    tuple val(meta), path("${prefix}.consensus.fasta"), optional: true, emit: consensus
    tuple val(meta), path("${prefix}.consensus.fasta"), val(meta.ref_id), optional: true, emit: consensus_with_ref

    script:
    prefix = "${meta.id}"
    """
    samtools mpileup \
        -aa -A -d 100 -B -Q 0 \
        -f ${reference} \
        ${dedup_bam} | \
    ivar consensus \
        -p ${prefix} \
        -q 20 \
        -m 10 \
        -n N

    # Format the consensus FASTA with proper line wrapping
    echo ">${prefix}" > "${prefix}.consensus.fasta"
    grep -v ">" "${prefix}.fa" | fold -w 70 >> "${prefix}.consensus.fasta"
    """
}
