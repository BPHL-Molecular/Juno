process SAMTOOLS {
    tag "${meta.id}"
    publishDir "${params.output_dir}/alignments", mode: 'copy', pattern: "*.{bam,bai}"
    publishDir "${params.output_dir}/stats", mode: 'copy', pattern: "*_{stats,flagstat,coverage,depth,markdup_stats}.txt"

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${prefix}.dedup.bam"), emit: dedup_bam
    tuple val(meta), path("${prefix}.dedup.bam"), path("${prefix}.dedup.bam.bai"), emit: indexed_bam
    tuple val(meta), path("${prefix}_stats.txt"), emit: stats
    tuple val(meta), path("${prefix}_flagstat.txt"), emit: flagstat
    tuple val(meta), path("${prefix}_coverage.txt"), emit: coverage
    tuple val(meta), path("${prefix}_depth.txt"), emit: depth
    tuple val(meta), path("${prefix}_markdup_stats.txt"), emit: markdup_stats

    script:
    prefix = "${meta.id}"
    """
    # Convert SAM to BAM and sort by read name
    samtools view -@ ${task.cpus} -bS ${sam} | samtools sort -@ ${task.cpus} -n -o ${prefix}.namesorted.bam

    # Fill in information about paired end reads
    samtools fixmate -@ ${task.cpus} -m ${prefix}.namesorted.bam ${prefix}.fixmate.bam

    # Sort by coordinate
    samtools sort -@ ${task.cpus} -o ${prefix}.sorted.bam ${prefix}.fixmate.bam

    # Mark and remove duplicates
    samtools markdup -@ ${task.cpus} -r -s ${prefix}.sorted.bam ${prefix}.dedup.bam 2> ${prefix}_markdup_stats.txt

    # Index deduplicated BAM
    samtools index -@ ${task.cpus} ${prefix}.dedup.bam

    # Generate stats on deduplicated BAM
    samtools stats ${prefix}.dedup.bam > ${prefix}_stats.txt
    samtools flagstat ${prefix}.dedup.bam > ${prefix}_flagstat.txt
    samtools coverage ${prefix}.dedup.bam > ${prefix}_coverage.txt
    samtools depth ${prefix}.dedup.bam > ${prefix}_depth.txt

    # Remove intermediate files to save space
    rm ${prefix}.namesorted.bam ${prefix}.fixmate.bam ${prefix}.sorted.bam
    """
}

process SAMTOOLS_DENOVO {
    tag "${meta.id}"
    publishDir "${params.output_dir}/samtools/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${prefix}_validation.sorted.bam"), path("${prefix}_validation.sorted.bam.bai"), emit: indexed_bam
    tuple val(meta), path("${prefix}_validation.sorted.bam"), emit: sorted_bam
    tuple val(meta), path("${prefix}_validation_coverage.txt"), emit: coverage
    tuple val(meta), path("${prefix}_validation_flagstat.txt"), emit: flagstat

    when:
    params.assembly_mode == 'denovo'

    script:
    prefix = "${meta.id}"
    """
    # Convert SAM to BAM and sort by read name
    samtools view -@ ${task.cpus} -bS ${sam} | samtools sort -@ ${task.cpus} -n -o ${prefix}_validation.namesorted.bam

    # Fill in information about paired end reads
    samtools fixmate -@ ${task.cpus} -m ${prefix}_validation.namesorted.bam ${prefix}_validation.fixmate.bam

    # Sort by coordinate
    samtools sort -@ ${task.cpus} -o ${prefix}_validation.sorted.bam ${prefix}_validation.fixmate.bam

    # Index BAM
    samtools index -@ ${task.cpus} ${prefix}_validation.sorted.bam

    # Generate essential statistics
    samtools coverage ${prefix}_validation.sorted.bam > ${prefix}_validation_coverage.txt
    samtools flagstat ${prefix}_validation.sorted.bam > ${prefix}_validation_flagstat.txt
    
    # Remove intermediate files
    rm ${prefix}_validation.namesorted.bam ${prefix}_validation.fixmate.bam
    """
}
