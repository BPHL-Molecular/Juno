process BLAST {
    tag "${meta.id}"
    publishDir "${params.output_dir}/${meta.id}/blast", mode: 'copy'

    input:
    tuple val(meta), path(contigs), path(reference)

    output:
    tuple val(meta), path("${meta.id}_blast_results.tsv"), emit: blast_results

    when:
    params.assembly_mode == 'denovo'

    script:
    prefix = "${meta.id}"
    """
    # Create blast database directory
    mkdir -p blast_db

    # Combine reference files into single FASTA in blast_db directory
    cat ${reference} > blast_db/orov_references.fasta

    # Create BLAST database
    makeblastdb \\
        -in blast_db/orov_references.fasta \\
        -dbtype nucl \\
        -out blast_db/orov_references \\
        -title "OROV_References"

    # Run BLAST against reference database
    blastn \\
        -query ${contigs} \\
        -db blast_db/orov_references \\
        -out ${prefix}_blast_results.tsv \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs" \\
        -max_target_seqs 1 \\
        -evalue 1e-10 \\
        -task megablast \\
        -num_threads ${task.cpus}
    """
}