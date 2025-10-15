process HRRT {
    tag "${meta.id}"
    publishDir "${params.output_dir}/dehosted", mode: 'copy'

    input:
    tuple val(meta), path(clean_reads)

    output:
    tuple val(meta), path("${prefix}_R{1,2}_dehosted.fastq.gz"), emit: dehosted_reads

    script:
    prefix = "${meta.id}"
    """
    (
        gzip -dc ${clean_reads[0]} | /opt/scrubber/scripts/scrub.sh -r -p ${task.cpus.intdiv(2)} -i - -o ${prefix}_R1_dehosted.fastq
        gzip ${prefix}_R1_dehosted.fastq
    ) &

    (
        gzip -dc ${clean_reads[1]} | /opt/scrubber/scripts/scrub.sh -r -p ${task.cpus.intdiv(2)} -i - -o ${prefix}_R2_dehosted.fastq
        gzip ${prefix}_R2_dehosted.fastq
    ) &

    wait
    """
}
