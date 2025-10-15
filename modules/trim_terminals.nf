process TRIM_TERMINALS {
    tag "${meta.id}"
    publishDir "${params.output_dir}/trimmed_contigs", mode: 'copy'

    input:
    tuple val(meta), path(polished_fasta), path(pilon_vcf)

    output:
    tuple val(meta), path("${prefix}_trimmed.fasta"), emit: trimmed_fasta
    tuple val(meta), path("${prefix}_trim_report.txt"), emit: trim_report

    when:
    params.polish_contigs

    script:
    prefix = "${meta.id}"
    """
#!/usr/bin/env python3
import sys
from collections import defaultdict

def parse_pilon_vcf(vcf_file):
    \"\"\"Parse Pilon VCF file and extract LowCov positions per contig.\"\"\"
    lowcov_positions = defaultdict(list)
    
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\\t')
            if len(fields) < 7:
                continue
            
            chrom = fields[0]
            pos = int(fields[1])
            filter_field = fields[6]
            
            if filter_field == "LowCov":
                lowcov_positions[chrom].append(pos)
    
    return lowcov_positions

def find_terminal_lowcov_regions(positions, contig_length):
    \"\"\"Identify only terminal (5' and 3') consecutive LowCov runs.\"\"\"
    if not positions:
        return 0, contig_length + 1
    
    sorted_pos = sorted(positions)
    
    # Find 5' terminal LowCov run
    trim_5prime = 0
    for i, pos in enumerate(sorted_pos):
        if pos == i + 1:
            trim_5prime = pos
        else:
            break
    
    # Find 3' terminal LowCov run
    trim_3prime = contig_length + 1
    for i, pos in enumerate(reversed(sorted_pos)):
        expected_pos = contig_length - i
        if pos == expected_pos:
            trim_3prime = pos
        else:
            break
    
    return trim_5prime, trim_3prime

def format_fasta_with_wrapping(sequences, output_file):
    \"\"\"Format FASTA sequences with 70-character line wrapping.\"\"\"
    with open(output_file, 'w') as f:
        if sequences:
            for header, sequence in sequences.items():
                f.write(f">{header}\\n")
                for i in range(0, len(sequence), 70):
                    f.write(sequence[i:i+70] + "\\n")

# Parse VCF for LowCov positions
lowcov_data = parse_pilon_vcf('${pilon_vcf}')

# Process each contig
trimmed_sequences = {}
report_lines = ["Contig\\tOriginal_Length\\tTrim_5prime\\tTrim_3prime\\tFinal_Length\\tStatus"]

with open('${polished_fasta}', 'r') as fasta_file:
    current_seq = ""
    current_id = ""
    
    for line in fasta_file:
        line = line.strip()
        if line.startswith('>'):
            if current_id and current_seq:
                # Process previous contig
                contig_name = current_id.split()[0]
                original_length = len(current_seq)
                
                # Strip _pilon suffix to match VCF contig names
                vcf_contig_name = contig_name.replace('_pilon', '')
                lowcov_pos = lowcov_data.get(vcf_contig_name, [])
                
                if not lowcov_pos:
                    # No LowCov positions, keep as-is
                    trimmed_sequences[current_id] = current_seq
                    report_lines.append(
                        f"{contig_name}\\t{original_length}\\t0\\t0\\t{original_length}\\tNo_Trimming_Needed"
                    )
                else:
                    # Find terminal LowCov regions
                    trim_5prime, trim_3prime = find_terminal_lowcov_regions(lowcov_pos, original_length)
                    
                    trimmed_5prime_bases = trim_5prime
                    trimmed_3prime_bases = original_length - (trim_3prime - 1) if trim_3prime <= original_length else 0
                    
                    if trimmed_5prime_bases == 0 and trimmed_3prime_bases == 0:
                        # LowCov positions are internal only, keep as-is
                        trimmed_sequences[current_id] = current_seq
                        report_lines.append(
                            f"{contig_name}\\t{original_length}\\t0\\t0\\t{original_length}\\tInternal_LowCov_Only"
                        )
                    else:
                        # Perform trimming
                        start = trim_5prime
                        end = trim_3prime - 1
                        trimmed_seq = current_seq[start:end]
                        final_length = len(trimmed_seq)
                        
                        # Update header with trim info
                        new_header = f"{current_id} trimmed_5prime:{trimmed_5prime_bases}bp trimmed_3prime:{trimmed_3prime_bases}bp final_length:{final_length}bp"
                        trimmed_sequences[new_header] = trimmed_seq
                        
                        report_lines.append(
                            f"{contig_name}\\t{original_length}\\t{trimmed_5prime_bases}\\t{trimmed_3prime_bases}\\t{final_length}\\tTrimmed"
                        )
                        
                        print(f"Trimmed {contig_name}: {original_length} -> {final_length} bp (5':{trimmed_5prime_bases} bp, 3':{trimmed_3prime_bases} bp)", file=sys.stderr)
            
            current_id = line[1:]
            current_seq = ""
        else:
            current_seq += line
    
    # Process last contig
    if current_id and current_seq:
        contig_name = current_id.split()[0]
        original_length = len(current_seq)
        
        # Strip _pilon suffix to match VCF contig names
        vcf_contig_name = contig_name.replace('_pilon', '')
        lowcov_pos = lowcov_data.get(vcf_contig_name, [])
        
        if not lowcov_pos:
            trimmed_sequences[current_id] = current_seq
            report_lines.append(
                f"{contig_name}\\t{original_length}\\t0\\t0\\t{original_length}\\tNo_Trimming_Needed"
            )
        else:
            trim_5prime, trim_3prime = find_terminal_lowcov_regions(lowcov_pos, original_length)
            
            trimmed_5prime_bases = trim_5prime
            trimmed_3prime_bases = original_length - (trim_3prime - 1) if trim_3prime <= original_length else 0
            
            if trimmed_5prime_bases == 0 and trimmed_3prime_bases == 0:
                trimmed_sequences[current_id] = current_seq
                report_lines.append(
                    f"{contig_name}\\t{original_length}\\t0\\t0\\t{original_length}\\tInternal_LowCov_Only"
                )
            else:
                start = trim_5prime
                end = trim_3prime - 1
                trimmed_seq = current_seq[start:end]
                final_length = len(trimmed_seq)
                
                new_header = f"{current_id} trimmed_5prime:{trimmed_5prime_bases}bp trimmed_3prime:{trimmed_3prime_bases}bp final_length:{final_length}bp"
                trimmed_sequences[new_header] = trimmed_seq
                
                report_lines.append(
                    f"{contig_name}\\t{original_length}\\t{trimmed_5prime_bases}\\t{trimmed_3prime_bases}\\t{final_length}\\tTrimmed"
                )
                
                print(f"Trimmed {contig_name}: {original_length} -> {final_length} bp (5':{trimmed_5prime_bases} bp, 3':{trimmed_3prime_bases} bp)", file=sys.stderr)

# Write output files
format_fasta_with_wrapping(trimmed_sequences, '${prefix}_trimmed.fasta')

with open('${prefix}_trim_report.txt', 'w') as report:
    report.write('\\n'.join(report_lines) + '\\n')

print(f"Done! Processed {len(trimmed_sequences)} contigs", file=sys.stderr)
"""
}
