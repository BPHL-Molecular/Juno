process SAMTOOLS {
    tag "${meta.id}"
    publishDir "${params.output_dir}/alignments", mode: 'copy', pattern: "*.{bam,bai}"
    publishDir "${params.output_dir}/stats", mode: 'copy', pattern: "*_{stats,flagstat,coverage,depth}.txt"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${prefix}.sorted.bam"), emit: sorted_bam
    tuple val(meta), path("${prefix}.sorted.bam"), path("${prefix}.sorted.bam.bai"), emit: indexed_bam
    tuple val(meta), path("${prefix}_stats.txt"), emit: stats
    tuple val(meta), path("${prefix}_flagstat.txt"), emit: flagstat
    tuple val(meta), path("${prefix}_coverage.txt"), emit: coverage
    tuple val(meta), path("${prefix}_depth.txt"), emit: depth

    script:
    prefix = "${meta.id}"
    """
    # Convert SAM to BAM
    samtools view -@ ${task.cpus} -bS ${sam} | samtools sort -@ ${task.cpus} -o ${prefix}.sorted.bam

    # Index BAM
    samtools index -@ ${task.cpus} ${prefix}.sorted.bam

    # Generate stats
    samtools stats ${prefix}.sorted.bam > ${prefix}_stats.txt
    samtools flagstat ${prefix}.sorted.bam > ${prefix}_flagstat.txt
    samtools coverage ${prefix}.sorted.bam > ${prefix}_coverage.txt
    samtools depth ${prefix}.sorted.bam > ${prefix}_depth.txt
    """
}
