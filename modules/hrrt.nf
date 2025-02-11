process HRRT {
    tag "${meta.id}"
    publishDir "${params.output_dir}/dehosted", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}_R{1,2}_cleaned.fastq.gz"), emit: dehosted_reads

    script:
    prefix = "${meta.id}"
    """
    # Decompress files
    gzip -d -c ${reads[0]} > ${prefix}_R1.fastq
    gzip -d -c ${reads[1]} > ${prefix}_R2.fastq

    # Run human read removal
    /opt/scrubber/scripts/scrub.sh \
        -r \
        -p ${task.cpus} \
        -i ${prefix}_R1.fastq \
        -o ${prefix}_R1_cleaned.fastq

    /opt/scrubber/scripts/scrub.sh \
        -r \
        -p ${task.cpus} \
        -i ${prefix}_R2.fastq \
        -o ${prefix}_R2_cleaned.fastq

    # Compress cleaned files
    gzip ${prefix}_R1_cleaned.fastq
    gzip ${prefix}_R2_cleaned.fastq

    # Cleanup intermediate files
    rm ${prefix}_R1.fastq ${prefix}_R2.fastq
    """
}
