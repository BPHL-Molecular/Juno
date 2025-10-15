process POLISH_CONTIGS {
    tag "${meta.id}"
    publishDir "${params.output_dir}/pilon/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(contigs), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}_polished.fasta"), emit: polished_fasta
    tuple val(meta), path("${prefix}_pilon.changes"), emit: changes
    tuple val(meta), path("${prefix}_pilon.vcf"), emit: vcf

    when:
    params.assembly_mode == 'denovo' && params.polish_contigs

    script:
    prefix = "${meta.id}"
    memory = task.memory.toGiga()
    """
    # Run Pilon to polish the assembly
    java -jar -Xmx${memory}G /pilon/pilon.jar \\
         --genome ${contigs} \\
         --bam ${bam} \\
         --output ${prefix}_pilon \\
         --changes \\
         --vcf \\
         --threads ${task.cpus} \\
         --fix all \\
         --mindepth 0.1
    
    # Rename the output fasta to a consistent naming scheme
    mv ${prefix}_pilon.fasta ${prefix}_polished.fasta
    """
}
