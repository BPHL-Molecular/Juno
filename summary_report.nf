process SUMMARY_REPORT {
    tag "SummaryReport"

    publishDir "${params.output_dir}", mode: 'copy'

    input:
        file fastp_files
        file kraken2_files
        file coverage_files
        file quast_dirs
        file variants_files
        file consensus_files
        file classification_files

    output:
        path "summary_report.tsv", emit: summary

    script:
    def assembly_mode = params.assembly_mode
    def min_coverage = params.qc_thresholds.min_coverage
    def min_depth = params.qc_thresholds.min_depth
    '''
#!/usr/bin/env python3
import os, glob, json

def debug_print(msg):
    print(f"DEBUG: {msg}", flush=True)

def get_orov_reads(kraken_report):
    try:
        with open(kraken_report) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 6 and "Orthobunyavirus oropoucheense" in fields[5]:
                    return int(fields[2])
    except Exception as e:
        debug_print(f"Error processing kraken report {kraken_report}: {e}")
    return 0

def get_classification_stats(classification_file):
    """Extract classification statistics from de novo classification summary."""
    stats = {'L': 0, 'M': 0, 'S': 0, 'unassigned': 0}
    try:
        with open(classification_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('L:'):
                    stats['L'] = int(line.split(':')[1].split()[0])
                elif line.startswith('M:'):
                    stats['M'] = int(line.split(':')[1].split()[0])
                elif line.startswith('S:'):
                    stats['S'] = int(line.split(':')[1].split()[0])
                elif line.startswith('unassigned:'):
                    stats['unassigned'] = int(line.split(':')[1].split()[0])
    except Exception as e:
        debug_print(f"Error processing classification file {classification_file}: {e}")
    return stats

# Assembly mode and QC thresholds
assembly_mode = "${assembly_mode}"
min_coverage = float("${min_coverage}")
min_depth = float("${min_depth}")
debug_print(f"Assembly mode: {assembly_mode}")
debug_print(f"QC thresholds - Coverage: {min_coverage}%, Depth: {min_depth}x")

# Index fastp files
fastp_dict = {}
for f in glob.glob("*.fastp.json"):
    sample_id = f.replace(".fastp.json", "")
    fastp_dict[sample_id] = f
    debug_print(f"Found FASTP file: {f} for sample {sample_id}")

# Index kraken2 files
kraken2_dict = {}
for f in glob.glob("*.kraken2.report"):
    sample_id = f.replace(".kraken2.report", "")
    kraken2_dict[sample_id] = f
    debug_print(f"Found Kraken2 report: {f} for sample {sample_id}")

# Index classification files (de novo mode only)
classification_dict = {}
if assembly_mode == "denovo":
    for f in glob.glob("*_classification_summary.txt"):
        sample_id = f.replace("_classification_summary.txt", "")
        classification_dict[sample_id] = f
        debug_print(f"Found classification file: {f} for sample {sample_id}")

if assembly_mode == "denovo":
    # DE NOVO ASSEMBLY REPORTING
    
    # Index quast files for de novo (per segment)
    quast_dict = {}
    for f in glob.glob("*_quast.report.tsv"):
        key = f.replace("_quast.report.tsv", "")
        quast_dict[key] = f
        debug_print(f"Found QUAST file: {f} for key {key}")

    # Prepare headers for de novo mode
    header = ["sampleID", "segment", "num_raw_reads", "num_clean_reads", "kraken2_reads",
              "num_contigs_L", "num_contigs_M", "num_contigs_S", "num_contigs_unassigned",
              "assembly_length", "reference_length", "qc_pass_fail"]
    out_lines = ["\t".join(header)]

    # Process each sample
    for sample_id in sorted(fastp_dict.keys()):
        debug_print(f"Processing sample: {sample_id}")

        # Get fastp data
        with open(fastp_dict[sample_id]) as f:
            fastp_data = json.load(f)
            num_raw_reads = fastp_data["summary"]["before_filtering"]["total_reads"]
            num_clean_reads = fastp_data["summary"]["after_filtering"]["total_reads"]

        # Get kraken2 data
        kraken2_reads = 0
        if sample_id in kraken2_dict:
            kraken2_reads = get_orov_reads(kraken2_dict[sample_id])

        # Get classification statistics
        classification_stats = {'L': 0, 'M': 0, 'S': 0, 'unassigned': 0}
        if sample_id in classification_dict:
            classification_stats = get_classification_stats(classification_dict[sample_id])

        # Process each segment
        for segment in ['L', 'M', 'S']:
            segment_key = f"{sample_id}_{segment}"
            
            # Get quast data for this segment
            assembly_length = reference_length = 0
            if segment_key in quast_dict:
                with open(quast_dict[segment_key]) as f:
                    for line in f:
                        if "Reference length" in line:
                            reference_length = int(line.strip().split("\t")[1])
                        elif line.startswith("Total length"):
                            assembly_length = int(line.strip().split("\t")[1])

            # Determine QC status (simplified for de novo - based on whether segment was assembled)
            qc_pass_fail = "PASS" if classification_stats[segment] > 0 and assembly_length > 0 else "FAIL"

            # Create output row
            row = [
                sample_id,
                segment,
                str(num_raw_reads),
                str(num_clean_reads),
                str(kraken2_reads),
                str(classification_stats['L']),
                str(classification_stats['M']),
                str(classification_stats['S']),
                str(classification_stats['unassigned']),
                str(assembly_length),
                str(reference_length),
                qc_pass_fail
            ]

            out_lines.append("\t".join(row))
            debug_print(f"Added row for {sample_id}_{segment}")

else:
    # REFERENCE-BASED ASSEMBLY REPORTING (existing logic)
    
    # Index coverage files
    coverage_dict = {}
    for f in glob.glob("*_coverage.txt"):
        key = f.replace("_coverage.txt", "")
        coverage_dict[key] = f
        debug_print(f"Found coverage file: {f} for key {key}")

    # Index variants files
    variants_dict = {}
    for f in glob.glob("*.variants.tsv"):
        key = f.replace(".variants.tsv", "")
        variants_dict[key] = f
        debug_print(f"Found variants file: {f} for key {key}")

    # Index consensus files
    consensus_dict = {}
    for f in glob.glob("*.consensus.fasta"):
        key = f.replace(".consensus.fasta", "")
        consensus_dict[key] = f
        debug_print(f"Found consensus file: {f} for key {key}")

    # Index quast files
    quast_dict = {}
    for f in glob.glob("*_quast.report.tsv"):
        key = f.replace("_quast.report.tsv", "")
        quast_dict[key] = f
        debug_print(f"Found QUAST file: {f} for key {key}")

    # Prepare headers and output lines
    header = ["sampleID", "reference", "num_raw_reads", "num_clean_reads", "kraken2_reads",
              "num_mapped_reads", "percent_mapped_clean_reads", "mean_base_qual", "mean_map_qual",
              "percent_reference_covered", "mean_ref_depth", "reference_length", "assembly_length",
              "num_variants", "total_n_bases", "qc_pass_fail"]
    out_lines = ["\t".join(header)]

    # Process each sample-reference combination
    for composite_key in sorted(coverage_dict.keys()):
        debug_print(f"Processing composite key: {composite_key}")

        parts = composite_key.rsplit("_", 1)
        if len(parts) != 2:
            debug_print(f"Skipping malformed key: {composite_key}")
            continue

        sample_id, ref_id = parts
        debug_print(f"Sample ID: {sample_id}, Reference ID: {ref_id}")

        # Get fastp data
        if sample_id not in fastp_dict:
            debug_print(f"No FASTP data for {sample_id}")
            continue

        with open(fastp_dict[sample_id]) as f:
            fastp_data = json.load(f)
            num_raw_reads = fastp_data["summary"]["before_filtering"]["total_reads"]
            num_clean_reads = fastp_data["summary"]["after_filtering"]["total_reads"]

        # Get kraken2 data
        kraken2_reads = 0
        if sample_id in kraken2_dict:
            kraken2_reads = get_orov_reads(kraken2_dict[sample_id])

        # Get coverage data
        with open(coverage_dict[composite_key]) as f:
            coverage_data = None
            for line in f:
                if not line.startswith("#"):
                    coverage_data = line.strip().split("\t")
                    break

        if not coverage_data:
            debug_print(f"No coverage data for {composite_key}")
            continue

        num_mapped_reads = float(coverage_data[3])
        percent_reference_covered = float(coverage_data[5])
        mean_ref_depth = float(coverage_data[6])
        mean_base_qual = float(coverage_data[7])
        mean_map_qual = float(coverage_data[8])
        percent_mapped_clean_reads = (num_mapped_reads / num_clean_reads * 100) if num_clean_reads > 0 else 0

        # Get quast data
        reference_length = assembly_length = 0
        if composite_key in quast_dict:
            with open(quast_dict[composite_key]) as f:
                for line in f:
                    if "Reference length" in line:
                        reference_length = int(line.strip().split("\t")[1])
                    elif line.startswith("Total length"):
                        assembly_length = int(line.strip().split("\t")[1])

        # Count variants
        num_variants = 0
        if composite_key in variants_dict:
            with open(variants_dict[composite_key]) as f:
                num_variants = sum(1 for line in f if not line.startswith("#"))

        # Count N bases in consensus
        total_n_bases = 0
        if composite_key in consensus_dict:
            with open(consensus_dict[composite_key]) as f:
                for line in f:
                    if not line.startswith(">"):
                        total_n_bases += line.upper().count("N")

        # Determine QC status
        qc_pass_fail = "PASS" if (percent_reference_covered >= min_coverage and mean_ref_depth >= min_depth) else "FAIL"

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
        debug_print(f"Added row for {composite_key}")

# Write output file
with open("summary_report.tsv", "w") as f:
    f.write(os.linesep.join(out_lines))
debug_print("Finished writing summary report")
'''
}
