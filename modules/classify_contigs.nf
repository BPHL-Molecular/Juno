process CLASSIFY_CONTIGS {
    tag "${meta.id}"
    publishDir "${params.output_dir}/assemblies/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(contigs), path(blast_results)

    output:
    tuple val(meta), path("${meta.id}_L.fasta"), emit: segment_L
    tuple val(meta), path("${meta.id}_M.fasta"), emit: segment_M
    tuple val(meta), path("${meta.id}_S.fasta"), emit: segment_S
    tuple val(meta), path("${meta.id}_unassigned.fasta"), emit: unassigned
    tuple val(meta), path("${meta.id}_classification_summary.txt"), emit: classification_summary

    when:
    params.assembly_mode == 'denovo'

    script:
    prefix = "${meta.id}"
    def identity_threshold = 85
    def coverage_threshold = 70
    """
#!/usr/bin/env python3
import os

def format_fasta_with_wrapping(sequences, output_file):
    \"\"\"Format FASTA sequences with 70-character line wrapping.\"\"\"
    with open(output_file, 'w') as f:
        if sequences:
            for header, sequence in sequences.items():
                f.write(f">{header}\\n")
                for i in range(0, len(sequence), 70):
                    f.write(sequence[i:i+70] + "\\n")

def parse_blast_results(blast_file):
    \"\"\"Parse BLAST results and return dictionary of best hits per contig.\"\"\"
    blast_hits = {}
    try:
        with open(blast_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                fields = line.split('\\t')
                if len(fields) < 15:
                    continue

                qseqid = fields[0]
                sseqid = fields[1]
                pident = float(fields[2])
                qcovs = float(fields[14])
                bitscore = float(fields[11])

                if qseqid not in blast_hits or bitscore > blast_hits[qseqid]['bitscore']:
                    blast_hits[qseqid] = {
                        'sseqid': sseqid,
                        'pident': pident,
                        'qcovs': qcovs,
                        'bitscore': bitscore
                    }

        return blast_hits
    except Exception as e:
        print(f"Warning: Could not parse BLAST results: {e}")
        return {}

def classify_contig(contig_id, blast_hit, identity_threshold, coverage_threshold):
    \"\"\"Classify a contig based on BLAST hit.\"\"\"
    if not blast_hit:
        return 'unassigned', 'No BLAST hit'

    identity = blast_hit['pident']
    coverage = blast_hit['qcovs']
    subject = blast_hit['sseqid']

    # Check thresholds
    if identity < identity_threshold:
        return 'unassigned', f'Low identity ({identity:.1f}%)'

    if coverage < coverage_threshold:
        return 'unassigned', f'Low coverage ({coverage:.1f}%)'

    if subject == 'PQ064919.1':
        segment = 'L'
    elif subject == 'PQ064920.1':
        segment = 'M'
    elif subject == 'PQ064921.1':
        segment = 'S'
    else:
        return 'unassigned', f'Unknown reference: {subject}'

    return segment, f'Identity: {identity:.1f}%, Coverage: {coverage:.1f}%, Reference: {subject}'

# Parse BLAST results
blast_hits = parse_blast_results('${blast_results}')

# Initialize output data structures
segment_sequences = {
    'L': {},
    'M': {},
    'S': {},
    'unassigned': {}
}

classification_stats = {'L': 0, 'M': 0, 'S': 0, 'unassigned': 0}
classification_details = []

# Read contigs and classify
with open('${contigs}', 'r') as fasta_file:
    current_seq = ""
    current_id = ""

    for line in fasta_file:
        line = line.strip()
        if line.startswith('>'):
            if current_id and current_seq:
                blast_hit = blast_hits.get(current_id)
                segment, reason = classify_contig(current_id, blast_hit, ${identity_threshold}, ${coverage_threshold})

                segment_sequences[segment][current_id] = current_seq
                classification_stats[segment] += 1
                classification_details.append([current_id, segment, reason, len(current_seq)])

            current_id = line[1:].split()[0]
            current_seq = ""
        else:
            current_seq += line

    if current_id and current_seq:
        blast_hit = blast_hits.get(current_id)
        segment, reason = classify_contig(current_id, blast_hit, ${identity_threshold}, ${coverage_threshold})
        segment_sequences[segment][current_id] = current_seq
        classification_stats[segment] += 1
        classification_details.append([current_id, segment, reason, len(current_seq)])

# Write formatted FASTA files for each segment
for segment in ['L', 'M', 'S', 'unassigned']:
    output_file = f'${prefix}_' + segment + '.fasta'
    format_fasta_with_wrapping(segment_sequences[segment], output_file)

# Write classification summary
with open('${prefix}_classification_summary.txt', 'w') as summary:
    summary.write("OROV Contig Classification Summary\\n")
    summary.write("=" * 40 + "\\n\\n")

    summary.write("Classification Statistics:\\n")
    for segment in ['L', 'M', 'S', 'unassigned']:
        count = classification_stats[segment]
        summary.write(f"  {segment}: {count} contigs\\n")

    total_contigs = sum(classification_stats.values())
    summary.write(f"\\nTotal contigs processed: {total_contigs}\\n\\n")

    summary.write("Detailed Classification:\\n")
    summary.write("Contig_ID\\tSegment\\tReason\\tLength\\n")
    for detail in classification_details:
        summary.write("\\t".join(map(str, detail)) + "\\n")
"""
}