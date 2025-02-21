process HRRT {
    tag "${meta.id}"
    publishDir "${params.output_dir}/dehosted", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}_R{1,2}_cleaned.fastq.gz"), emit: dehosted_reads

    script:
    prefix = "${meta.id}"
    if (params.parallel_hrrt)
        """
        # Process R1 and R2 in parallel (HPC environment)
        (
            gzip -d -c ${reads[0]} > ${prefix}_R1.fastq
            /opt/scrubber/scripts/scrub.sh -r -p ${task.cpus} -i ${prefix}_R1.fastq -o ${prefix}_R1_cleaned.fastq
            gzip ${prefix}_R1_cleaned.fastq
            rm ${prefix}_R1.fastq
        ) &

        (
            gzip -d -c ${reads[1]} > ${prefix}_R2.fastq
            /opt/scrubber/scripts/scrub.sh -r -p ${task.cpus} -i ${prefix}_R2.fastq -o ${prefix}_R2_cleaned.fastq
            gzip ${prefix}_R2_cleaned.fastq
            rm ${prefix}_R2.fastq
        ) &

        wait
        """
    else
        """
        # Process R1 and R2 sequentially
        gzip -d -c ${reads[0]} > ${prefix}_R1.fastq
        /opt/scrubber/scripts/scrub.sh -r -p ${task.cpus} -i ${prefix}_R1.fastq -o ${prefix}_R1_cleaned.fastq
        gzip ${prefix}_R1_cleaned.fastq
        rm ${prefix}_R1.fastq

        gzip -d -c ${reads[1]} > ${prefix}_R2.fastq
        /opt/scrubber/scripts/scrub.sh -r -p ${task.cpus} -i ${prefix}_R2.fastq -o ${prefix}_R2_cleaned.fastq
        gzip ${prefix}_R2_cleaned.fastq
        rm ${prefix}_R2.fastq
        """
}
