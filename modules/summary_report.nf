process SUMMARY_REPORT_REFERENCE {
    tag "SummaryReport"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        file fastp_files
        file kraken2_files
        file coverage_files
        file quast_dirs
        file variants_files
        file consensus_files

    output:
        path "summary_report.tsv", emit: summary

    when:
    params.assembly_mode == 'reference'

    script:
    '''
#!/usr/bin/env python3
import os, glob, json

def get_orov_reads(kraken_report):
    \"\"\"Get OROV reads count from kraken2 report - using parent taxon count.\"\"\"
    try:
        with open(kraken_report) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 6 and "Orthobunyavirus oropoucheense" in fields[5]:
                    return int(fields[1])
    except Exception as e:
        print(f"Error processing kraken report {kraken_report}: {e}")
    return 0

# Index fastp files
fastp_dict = {}
for f in glob.glob("*.fastp.json"):
    sample_id = f.replace(".fastp.json", "")
    fastp_dict[sample_id] = f

# Index kraken2 files
kraken2_dict = {}
for f in glob.glob("*.kraken2.report"):
    sample_id = f.replace(".kraken2.report", "")
    kraken2_dict[sample_id] = f

# Index coverage files
coverage_dict = {}
for f in glob.glob("*_coverage.txt"):
    key = f.replace("_coverage.txt", "")
    coverage_dict[key] = f

# Index variants files
variants_dict = {}
for f in glob.glob("*.variants.tsv"):
    key = f.replace(".variants.tsv", "")
    variants_dict[key] = f

# Index consensus files
consensus_dict = {}
for f in glob.glob("*.consensus.fasta"):
    key = f.replace(".consensus.fasta", "")
    consensus_dict[key] = f

# Index quast files
quast_dict = {}
for f in glob.glob("*_quast.report.tsv"):
    key = f.replace("_quast.report.tsv", "")
    quast_dict[key] = f

# QC thresholds
min_coverage = 90.0
min_depth = 15.0
high_n_threshold = 5.0

# Prepare output
header = ["sampleID", "reference", "num_raw_reads", "num_clean_reads", "kraken2_reads",
          "num_mapped_reads", "percent_mapped_clean_reads", "mean_base_qual", "mean_map_qual",
          "percent_reference_covered", "mean_ref_depth", "reference_length", "assembly_length",
          "num_variants", "total_n_bases", "qc_pass_fail"]
out_lines = ["\t".join(header)]

# Process each sample-reference combination
for composite_key in sorted(coverage_dict.keys()):
    parts = composite_key.rsplit("_", 1)
    if len(parts) != 2:
        continue

    sample_id, ref_id = parts

    if sample_id not in fastp_dict:
        continue

    with open(fastp_dict[sample_id]) as f:
        fastp_data = json.load(f)
        num_raw_reads = fastp_data["summary"]["before_filtering"]["total_reads"]
        num_clean_reads = fastp_data["summary"]["after_filtering"]["total_reads"]

    kraken2_reads = 0
    if sample_id in kraken2_dict:
        kraken2_reads = get_orov_reads(kraken2_dict[sample_id])

    with open(coverage_dict[composite_key]) as f:
        coverage_data = None
        for line in f:
            if not line.startswith("#"):
                coverage_data = line.strip().split("\t")
                break

    if not coverage_data:
        continue

    num_mapped_reads = float(coverage_data[3])
    percent_reference_covered = float(coverage_data[5])
    mean_ref_depth = float(coverage_data[6])
    mean_base_qual = float(coverage_data[7])
    mean_map_qual = float(coverage_data[8])
    percent_mapped_clean_reads = (num_mapped_reads / num_clean_reads * 100) if num_clean_reads > 0 else 0

    reference_length = assembly_length = 0
    if composite_key in quast_dict:
        with open(quast_dict[composite_key]) as f:
            for line in f:
                if "Reference length" in line:
                    reference_length = int(line.strip().split("\t")[1])
                elif line.startswith("Total length"):
                    assembly_length = int(line.strip().split("\t")[1])

    num_variants = 0
    if composite_key in variants_dict:
        with open(variants_dict[composite_key]) as f:
            num_variants = sum(1 for line in f if not line.startswith("#"))

    total_n_bases = 0
    if composite_key in consensus_dict:
        with open(consensus_dict[composite_key]) as f:
            for line in f:
                if not line.startswith(">"):
                    total_n_bases += line.upper().count("N")

    # Determine QC status
    n_base_percentage = (total_n_bases / assembly_length * 100) if assembly_length > 0 else 0
    basic_pass = (percent_reference_covered >= min_coverage and mean_ref_depth >= min_depth)
    
    if not basic_pass:
        qc_pass_fail = "FAIL"
    elif n_base_percentage > high_n_threshold:
        qc_pass_fail = "PASS_W_HIGH_N_BASES"
    else:
        qc_pass_fail = "PASS"

    # Create output row
    row = [
        composite_key,
        ref_id,
        str(num_raw_reads),
        str(num_clean_reads),
        str(kraken2_reads),
        str(int(num_mapped_reads)),
        f"{percent_mapped_clean_reads:.2f}",
        f"{mean_base_qual:.1f}",
        f"{mean_map_qual:.1f}",
        f"{percent_reference_covered:.1f}",
        f"{mean_ref_depth:.2f}",
        str(reference_length),
        str(assembly_length),
        str(num_variants),
        str(total_n_bases),
        qc_pass_fail
    ]

    out_lines.append("\t".join(row))

with open("summary_report.tsv", "w") as f:
    f.write(os.linesep.join(out_lines))
'''
}

process SUMMARY_REPORT_DENOVO {
    tag "SummaryReport_Denovo"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        file fastp_files
        file kraken2_files
        file quast_dirs
        file classification_files

    output:
        path "summary_report.tsv", emit: summary

    when:
    params.assembly_mode == 'denovo'

    script:
    '''
#!/usr/bin/env python3
import os, glob, json

def get_orov_reads(kraken_report):
    \"\"\"Get OROV reads count from kraken2 report - using parent taxon count.\"\"\"
    try:
        with open(kraken_report) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 6 and "Orthobunyavirus oropoucheense" in fields[5]:
                    return int(fields[1])
    except Exception as e:
        print(f"Error processing kraken report {kraken_report}: {e}")
    return 0

def get_classification_stats(classification_file):
    \"\"\"Extract classification statistics from de novo classification summary.\"\"\"
    stats = {'L': 0, 'M': 0, 'S': 0, 'unassigned': 0}
    try:
        with open(classification_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Look for lines containing segment info and contigs
                if 'L:' in line and 'contigs' in line:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part.isdigit():
                            stats['L'] = int(part)
                            break
                            
                elif 'M:' in line and 'contigs' in line:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part.isdigit():
                            stats['M'] = int(part)
                            break
                            
                elif 'S:' in line and 'contigs' in line:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part.isdigit():
                            stats['S'] = int(part)
                            break
                            
                elif 'unassigned:' in line and 'contigs' in line:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part.isdigit():
                            stats['unassigned'] = int(part)
                            break
                            
    except Exception as e:
        print(f"Error processing classification file {classification_file}: {e}")
    return stats

def get_blast_quality_from_classification(classification_file, segment):
    \"\"\"Extract BLAST quality metrics for a specific segment from classification summary.\"\"\"
    try:
        with open(classification_file, 'r') as f:
            in_details = False
            for line in f:
                line = line.strip()
                if line.startswith('Detailed Classification:'):
                    in_details = True
                    continue
                elif in_details and line and not line.startswith('Contig_ID'):
                    parts = line.split('\t')
                    if len(parts) >= 4:
                        contig_segment = parts[1]
                        reason = parts[2]
                        
                        if contig_segment == segment and 'Identity:' in reason:
                            identity = float(reason.split('Identity: ')[1].split('%')[0])
                            coverage = float(reason.split('Coverage: ')[1].split('%')[0])
                            return identity, coverage
    except Exception as e:
        print(f"Error extracting BLAST quality for {segment}: {e}")
    return 0, 0

def get_assembly_status(segment, classification_stats, quast_data):
    \"\"\"Determine assembly status for de novo assembly.\"\"\"
    
    # Check if segment was assembled
    if classification_stats[segment] == 0:
        return "NO_ASSEMBLY"
    
    # Check if we have assembly data
    assembly_length = quast_data.get('assembly_length', 0)
    if assembly_length == 0:
        return "NO_ASSEMBLY"
    
    return "ASSEMBLED"

# Index files
fastp_dict = {}
for f in glob.glob("*.fastp.json"):
    sample_id = f.replace(".fastp.json", "")
    fastp_dict[sample_id] = f

kraken2_dict = {}
for f in glob.glob("*.kraken2.report"):
    sample_id = f.replace(".kraken2.report", "")
    kraken2_dict[sample_id] = f

classification_dict = {}
for f in glob.glob("*_classification_summary.txt"):
    sample_id = f.replace("_classification_summary.txt", "")
    classification_dict[sample_id] = f

quast_dict = {}
for f in glob.glob("*_quast.report.tsv"):
    key = f.replace("_quast.report.tsv", "")
    quast_dict[key] = f

# Prepare output
header = ["sampleID", "segment", "num_raw_reads", "num_clean_reads", "kraken2_reads",
          "num_contigs_L", "num_contigs_M", "num_contigs_S",
          "assembly_length", "reference_length", "na50", "blast_identity", "blast_coverage", "assembly_status"]
out_lines = ["\t".join(header)]

# Process each sample
for sample_id in sorted(fastp_dict.keys()):
    with open(fastp_dict[sample_id]) as f:
        fastp_data = json.load(f)
        num_raw_reads = fastp_data["summary"]["before_filtering"]["total_reads"]
        num_clean_reads = fastp_data["summary"]["after_filtering"]["total_reads"]

    kraken2_reads = 0
    if sample_id in kraken2_dict:
        kraken2_reads = get_orov_reads(kraken2_dict[sample_id])

    classification_stats = {'L': 0, 'M': 0, 'S': 0, 'unassigned': 0}
    if sample_id in classification_dict:
        classification_stats = get_classification_stats(classification_dict[sample_id])

    for segment in ['L', 'M', 'S']:
        segment_key = f"{sample_id}_{segment}"
        
        quast_data = {}
        if segment_key in quast_dict:
            with open(quast_dict[segment_key]) as f:
                for line in f:
                    if "Reference length" in line:
                        quast_data['reference_length'] = int(line.strip().split("\t")[1])
                    elif line.startswith("Total length"):
                        quast_data['assembly_length'] = int(line.strip().split("\t")[1])
                    elif "NA50" in line and not "NGA50" in line:
                        quast_data['na50'] = int(line.strip().split("\t")[1])

        blast_identity, blast_coverage = 0, 0
        if sample_id in classification_dict:
            blast_identity, blast_coverage = get_blast_quality_from_classification(
                classification_dict[sample_id], segment
            )

        assembly_status = get_assembly_status(segment, classification_stats, quast_data)

        num_contigs_L = classification_stats['L'] if segment == 'L' else 0
        num_contigs_M = classification_stats['M'] if segment == 'M' else 0
        num_contigs_S = classification_stats['S'] if segment == 'S' else 0

        row = [
            sample_id,
            segment,
            str(num_raw_reads),
            str(num_clean_reads),
            str(kraken2_reads),
            str(num_contigs_L),
            str(num_contigs_M),
            str(num_contigs_S),
            str(quast_data.get('assembly_length', 0)),
            str(quast_data.get('reference_length', 0)),
            str(quast_data.get('na50', 0)),
            f"{blast_identity:.1f}",
            f"{blast_coverage:.1f}",
            assembly_status
        ]

        out_lines.append("\t".join(row))

with open("summary_report.tsv", "w") as f:
    f.write(os.linesep.join(out_lines))
'''
}
