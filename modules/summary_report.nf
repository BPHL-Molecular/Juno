process SUMMARY_REPORT_REFERENCE {
    tag "SummaryReport"
    publishDir "${params.output_dir}", mode: 'copy'

    input:
        file fastq_scan_files
        file fastp_files
        file kraken2_files
        file coverage_files
        file quast_dirs
        file variants_files
        file consensus_files
        file markdup_stats_files

    output:
        path "summary_report.tsv", emit: summary

    when:
    params.assembly_mode == 'reference'

    script:
    '''
#!/usr/bin/env python3
import os, glob, json

def get_orov_reads(kraken_report):
    """Get OROV reads count from kraken2 report - using parent taxon count."""
    try:
        with open(kraken_report) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 6 and "Orthobunyavirus oropoucheense" in fields[5]:
                    return int(fields[1])
    except Exception as e:
        print(f"Error processing kraken report {kraken_report}: {e}")
    return 0

def get_duplicate_percent(markdup_stats_file):
    """Extract duplicate percentage from samtools markdup stats file."""
    try:
        duplicate_total = 0
        written_total = 0
        
        with open(markdup_stats_file, 'r') as f:
            for line in f:
                if "DUPLICATE TOTAL:" in line:
                    parts = line.strip().split(":")
                    if len(parts) >= 2:
                        duplicate_total = int(parts[1].strip())
                elif "WRITTEN:" in line:
                    parts = line.strip().split(":")
                    if len(parts) >= 2:
                        written_total = int(parts[1].strip())
        
        if written_total > 0:
            return (duplicate_total / written_total) * 100
    except Exception as e:
        print(f"Error processing markdup stats file {markdup_stats_file}: {e}")
    return 0.0

# Index fastq-scan files (raw read counts)
fastq_scan_dict = {}
for f in glob.glob("*.fastq-scan.json"):
    sample_id = f.replace(".fastq-scan.json", "")
    fastq_scan_dict[sample_id] = f

# Index fastp files (clean read counts)
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

# Index markdup stats files
markdup_stats_dict = {}
for f in glob.glob("*_markdup_stats.txt"):
    key = f.replace("_markdup_stats.txt", "")
    markdup_stats_dict[key] = f

# QC thresholds
min_coverage = 90.0
min_depth = 15.0
high_n_threshold = 5.0

# Prepare output
header = ["sampleID", "reference", "num_raw_reads", "num_clean_reads", "kraken2_reads",
          "num_mapped_reads", "percent_mapped_clean_reads", "duplicate_percent", "mean_base_qual", "mean_map_qual",
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

    # Get raw reads from fastq-scan
    num_raw_reads = 0
    if sample_id in fastq_scan_dict:
        with open(fastq_scan_dict[sample_id]) as f:
            fastq_scan_data = json.load(f)
            num_raw_reads = fastq_scan_data["qc_stats"]["read_total"]

    # Get clean reads from fastp
    with open(fastp_dict[sample_id]) as f:
        fastp_data = json.load(f)
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

    # Get duplicate percentage
    duplicate_percent = 0.0
    if composite_key in markdup_stats_dict:
        duplicate_percent = get_duplicate_percent(markdup_stats_dict[composite_key])

    # Create output row
    row = [
        composite_key,
        ref_id,
        str(num_raw_reads),
        str(num_clean_reads),
        str(kraken2_reads),
        str(int(num_mapped_reads)),
        f"{percent_mapped_clean_reads:.2f}",
        f"{duplicate_percent:.2f}",
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
        file fastq_scan_files
        file fastp_files
        file kraken2_files
        file quast_dirs
        file classification_files
        file validation_coverage_files
        file validation_flagstat_files

    output:
        path "summary_report.tsv", emit: summary

    when:
    params.assembly_mode == 'denovo'

    script:
    '''
#!/usr/bin/env python3
import os, glob, json

def get_orov_reads(kraken_report):
    """Get OROV reads count from kraken2 report - using parent taxon count."""
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
    """Extract classification statistics from de novo classification summary."""
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
    """Extract BLAST quality metrics for a specific segment from classification summary."""
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
    """Determine assembly status for de novo assembly with quality assessment."""
    
    # Check 1: Are there contigs classified to this segment?
    if classification_stats[segment] == 0:
        return "NO_ASSEMBLY"
    
    # Check 2: Does QUAST have valid largest_contig?
    largest_contig = quast_data.get('largest_contig', 0)
    if largest_contig == 0:
        return "NO_ASSEMBLY"
    
    # Check 3: Quality assessment - compare largest contig to reference length
    reference_length = quast_data.get('reference_length', 0)
    length_ratio = (largest_contig / reference_length) * 100
    
    # Accept assemblies where largest contig is 90-150% of reference length
    if 90.0 <= length_ratio <= 150.0:
        return "ASSEMBLED"
    else:
        return "FRAGMENTED"

def get_validation_stats(coverage_file, flagstat_file):
    """Extract overall validation statistics from coverage and flagstat files."""
    stats = {'percent_mapped': 0.0, 'mean_depth': 0.0, 'mapped_reads': 0, 'total_reads': 0}
    
    try:
        # Process flagstat file for overall mapping stats
        with open(flagstat_file, 'r') as f:
            for line in f:
                if "mapped (" in line and "primary" not in line:
                    fields = line.strip().split()
                    stats['mapped_reads'] = int(fields[0])
                    stats['percent_mapped'] = float(fields[4].strip('()%'))
                elif "in total" in line:
                    fields = line.strip().split()
                    stats['total_reads'] = int(fields[0])
    except Exception as e:
        print(f"Error processing validation files: {e}")
        
    return stats

def get_contigs_by_segment(classification_file):
    """Extract mapping of contigs to segments from classification summary file."""
    segment_contigs = {'L': [], 'M': [], 'S': [], 'unassigned': []}
    
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
                        contig_id = parts[0].strip()
                        segment = parts[1].strip()
                        if segment in segment_contigs:
                            segment_contigs[segment].append(contig_id)
    except Exception as e:
        print(f"Error parsing classification file {classification_file}: {e}")
    
    return segment_contigs

def parse_coverage_by_contig(coverage_file):
    """Parse coverage data by contig from samtools coverage output."""
    contig_coverage = {}
    
    try:
        with open(coverage_file, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    fields = line.strip().split('\t')
                    if len(fields) >= 7:
                        contig_id = fields[0].strip()
                        coverage_pct = float(fields[5])
                        mean_depth = float(fields[6])
                        num_reads = int(fields[3])
                        contig_coverage[contig_id] = {
                            'coverage_pct': coverage_pct,
                            'mean_depth': mean_depth,
                            'num_reads': num_reads
                        }
    except Exception as e:
        print(f"Error parsing coverage file {coverage_file}: {e}")
    
    return contig_coverage

def calculate_segment_validation_stats(segment_contigs, contig_coverage, flagstat_data):
    """Calculate validation statistics per segment based on contig assignments."""
    segment_stats = {}
    total_mapped_reads = flagstat_data.get('mapped_reads', 0)
    
    for segment, contigs in segment_contigs.items():
        if segment != 'unassigned' and contigs:
            # Strip _pilon suffix from contig names to match validation coverage file
            # Classification uses polished names, validation uses original SPAdes names
            contigs_no_suffix = [contig.replace('_pilon', '') for contig in contigs]
            
            # Sum metrics for all contigs in this segment
            segment_reads = sum(contig_coverage.get(contig, {}).get('num_reads', 0) for contig in contigs_no_suffix)
            
            # Extract lengths from contig IDs (use original contig names without suffix)
            contig_lengths = {}
            for contig in contigs_no_suffix:
                parts = contig.split('_')
                if len(parts) > 3 and 'length' in parts[2]:
                    try:
                        contig_lengths[contig] = int(parts[3])
                    except (ValueError, IndexError):
                        contig_lengths[contig] = 0
            
            total_length = sum(contig_lengths.values())
            
            # Calculate weighted average depth
            weighted_depth = 0
            if total_length > 0:
                for contig in contigs_no_suffix:
                    if contig in contig_coverage and contig in contig_lengths:
                        weighted_depth += contig_coverage[contig]['mean_depth'] * (contig_lengths[contig] / total_length)
            
            # Calculate percentage of reads mapped to this segment
            segment_pct = (segment_reads / total_mapped_reads * 100) if total_mapped_reads > 0 else 0
            
            segment_stats[segment] = {
                'mapped_reads': segment_reads,
                'percent_mapped': segment_pct,
                'mean_depth': weighted_depth
            }
        else:
            # Default values for segments with no contigs
            segment_stats[segment] = {
                'mapped_reads': 0,
                'percent_mapped': 0.0,
                'mean_depth': 0.0
            }
    
    return segment_stats

# Index fastq-scan files (raw read counts)
fastq_scan_dict = {}
for f in glob.glob("*.fastq-scan.json"):
    sample_id = f.replace(".fastq-scan.json", "")
    fastq_scan_dict[sample_id] = f

# Index fastp files (clean read counts)
fastp_dict = {}
for f in glob.glob("*.fastp.json"):
    sample_id = f.replace(".fastp.json", "")
    fastp_dict[sample_id] = f

# Index kraken2 files
kraken2_dict = {}
for f in glob.glob("*.kraken2.report"):
    sample_id = f.replace(".kraken2.report", "")
    kraken2_dict[sample_id] = f

# Index classification files
classification_dict = {}
for f in glob.glob("*_classification_summary.txt"):
    sample_id = f.replace("_classification_summary.txt", "")
    classification_dict[sample_id] = f

# Index quast files
quast_dict = {}
for f in glob.glob("*_quast.report.tsv"):
    key = f.replace("_quast.report.tsv", "")
    quast_dict[key] = f

# Index validation files
validation_coverage_dict = {}
for f in glob.glob("*_validation_coverage.txt"):
    sample_id = f.replace("_validation_coverage.txt", "")
    validation_coverage_dict[sample_id] = f

validation_flagstat_dict = {}
for f in glob.glob("*_validation_flagstat.txt"):
    sample_id = f.replace("_validation_flagstat.txt", "")
    validation_flagstat_dict[sample_id] = f

# Prepare output
header = ["sampleID", "segment", "num_raw_reads", "num_clean_reads", "kraken2_reads",
          "num_contigs_L", "num_contigs_M", "num_contigs_S",
          "largest_contig", "reference_length", "na50", "blast_identity", "blast_coverage", 
          "validation_mapped_reads", "validation_percent_mapped", "validation_mean_depth",
          "assembly_status"]
out_lines = ["\t".join(header)]

# Process each sample
for sample_id in sorted(fastp_dict.keys()):
    # Get raw reads from fastq-scan
    num_raw_reads = 0
    if sample_id in fastq_scan_dict:
        with open(fastq_scan_dict[sample_id]) as f:
            fastq_scan_data = json.load(f)
            num_raw_reads = fastq_scan_data["qc_stats"]["read_total"]

    # Get clean reads from fastp
    with open(fastp_dict[sample_id]) as f:
        fastp_data = json.load(f)
        num_clean_reads = fastp_data["summary"]["after_filtering"]["total_reads"]

    kraken2_reads = 0
    if sample_id in kraken2_dict:
        kraken2_reads = get_orov_reads(kraken2_dict[sample_id])

    classification_stats = {'L': 0, 'M': 0, 'S': 0, 'unassigned': 0}
    if sample_id in classification_dict:
        classification_stats = get_classification_stats(classification_dict[sample_id])

    # Get contig-to-segment mappings
    segment_contigs = {'L': [], 'M': [], 'S': [], 'unassigned': []}
    if sample_id in classification_dict:
        segment_contigs = get_contigs_by_segment(classification_dict[sample_id])
    
    # Get contig-specific coverage data
    contig_coverage = {}
    if sample_id in validation_coverage_dict:
        contig_coverage = parse_coverage_by_contig(validation_coverage_dict[sample_id])
    
    # Get overall validation stats from flagstat
    overall_validation_stats = {'mapped_reads': 0, 'total_reads': 0, 'percent_mapped': 0.0}
    if sample_id in validation_flagstat_dict and sample_id in validation_coverage_dict:
        overall_validation_stats = get_validation_stats(
            validation_coverage_dict[sample_id],
            validation_flagstat_dict[sample_id]
        )
    
    # Calculate segment-specific validation statistics
    segment_validation = calculate_segment_validation_stats(
        segment_contigs,
        contig_coverage,
        overall_validation_stats
    )
    
    for segment in ['L', 'M', 'S']:
        segment_key = f"{sample_id}_{segment}"
        
        quast_data = {}
        if segment_key in quast_dict:
            with open(quast_dict[segment_key]) as f:
                for line in f:
                    if "Reference length" in line:
                        quast_data['reference_length'] = int(line.strip().split("\t")[1])
                    elif line.startswith("Largest contig"):
                        quast_data['largest_contig'] = int(line.strip().split("\t")[1])
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
        
        # Use segment-specific validation metrics
        segment_validation_data = segment_validation.get(segment, {'mapped_reads': 0, 'percent_mapped': 0.0, 'mean_depth': 0.0})
        
        row = [
            sample_id,
            segment,
            str(num_raw_reads),
            str(num_clean_reads),
            str(kraken2_reads),
            str(num_contigs_L),
            str(num_contigs_M),
            str(num_contigs_S),
            str(quast_data.get('largest_contig', 0)),
            str(quast_data.get('reference_length', 0)),
            str(quast_data.get('na50', 0)),
            f"{blast_identity:.1f}",
            f"{blast_coverage:.1f}",
            f"{int(segment_validation_data['mapped_reads'])}",
            f"{segment_validation_data['percent_mapped']:.2f}",
            f"{segment_validation_data['mean_depth']:.2f}",
            assembly_status
        ]

        out_lines.append("\t".join(row))

with open("summary_report.tsv", "w") as f:
    f.write(os.linesep.join(out_lines))
'''
}
