process MAKEBLASTDB {
    tag "OROV_References"
    publishDir "${params.output_dir}/blast/blast_db", mode: 'copy'

    input:
    path(references)

    output:
    path("orov_references*"), emit: db_files

    script:
    """
    # Combine all reference files
    cat ${references.join(' ')} > orov_references.fasta

    # Create BLAST database
    makeblastdb \\
        -in orov_references.fasta \\
        -dbtype nucl \\
        -out orov_references \\
        -title "OROV_References"
    """
}

process BLAST {
    tag "${meta.id}"
    publishDir "${params.output_dir}/blast/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(contigs)
    path(db_files)

    output:
    tuple val(meta), path("${meta.id}_blast_results.tsv"), emit: blast_results

    when:
    params.assembly_mode == 'denovo'

    script:
    prefix = "${meta.id}"
    """
    # Run BLAST against pre-built database
    blastn \\
        -query ${contigs} \\
        -db orov_references \\
        -out ${prefix}_blast_results.tsv \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \\
        -max_target_seqs 10 \\
        -evalue 1e-10 \\
        -task megablast \\
        -num_threads ${task.cpus}
    """
}
