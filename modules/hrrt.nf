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
    (
        gzip -dc ${reads[0]} | /opt/scrubber/scripts/scrub.sh -r -p ${task.cpus.intdiv(2)} -i - -o ${prefix}_R1_cleaned.fastq
        gzip ${prefix}_R1_cleaned.fastq
    ) &

    (
        gzip -dc ${reads[1]} | /opt/scrubber/scripts/scrub.sh -r -p ${task.cpus.intdiv(2)} -i - -o ${prefix}_R2_cleaned.fastq
        gzip ${prefix}_R2_cleaned.fastq
    ) &

    wait
    """
}
