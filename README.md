# Juno 🦟🦠🧬📊
**A Nextflow Pipeline for Reference-Based Assembly of Oropouche Virus (OROV)**

A containerized Nextflow pipeline designed for processing Illumina paired-end metagenomics sequencing data against OROV reference genomes, performing QC, taxonomic classification, variant calling, and consensus generation.

## ⚡ Usage
```bash
$ nextflow run juno.nf -profile singularity -params-file params.yaml
```

## 🐊 HiPerGator Usage
```bash
$ sbatch ./juno.sh
```

## 📦 Dependencies
- [Nextflow 23.04.0+](https://www.nextflow.io/docs/latest/install.html)
- [Singularity](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html#quick-installation-steps) or [Docker](https://docs.docker.com/engine/install/)
- [Python 3.6+](https://docs.python.org/3/using/unix.html)
- [Slurm](https://slurm.schedmd.com/documentation.html) (only if HiPerGator will be used)

## ⚙️ Configuration

**Important:** All pipeline parameters **must be set in the `params.yaml` file**. Make sure you edit this file to provide the correct paths and values before running the pipeline. You will also need to download the kraken2 viral database from the [BenLangmead Index zone](https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20241228.tar.gz) link.

### Required Parameters (in `params.yaml`)

Below is a sample configuration. **Modify these values as needed:**

```yaml
# Directory containing paired-end FASTQ files.
# Files must follow the naming convention including `_R1_` and `_R2_` to denote forward and reverse reads.
input_dir: "/path/to/fastq" 

# Directory where all pipeline results will be stored.
output_dir: "/path/to/output_dir"

# Path to the Kraken2 viral database.
kraken2_db: "/path/to/kraken2_db"

# Number of CPU threads to use.
threads: 8

# Quality control thresholds
qc_thresholds:
    min_coverage: 90
    min_depth: 15
```

## 🛠️ Pipeline Steps
1. **Quality Control**
   - Human Read Removal (`sra-human-scrubber`)
   - Read QC and trimming (`fastp`)
2. **Taxonomic Classification**
   - Read classification (`kraken2`)
3. **Assembly**
   - Reference alignment (`minimap2`)
   - SAM/BAM processing (`samtools`)
   - Variant calling & consensus (`ivar`)
4. **Quality Assessment**
   - Assembly evaluation (`quast`)
   - Report generation (`multiqc`)

## 📂 Output Structure
```
output_dir/
├── dehosted/         # Cleaned reads
├── trimmed/          # Trimmed reads
├── kraken2/          # Classification results
├── alignments/       # BAM files & indices
├── stats/            # Alignment statistics
├── variants/         # Variant calls
├── consensus/        # Consensus sequences
├── quast/            # Assembly metrics
├── multiqc/          # Combined QC report
└── summary_report.tsv
```

## 📋 Summary Report Metrics
- Sample and reference identifiers
- Cleaned read counts
- Mapping statistics
- Coverage metrics
- Variant counts
- Assembly quality metrics
- Overall QC status

## 🐛 Troubleshooting
**Pipeline Errors:**  
   Check Nextflow execution logs in .nextflow.log 
   
**Low Coverage Regions:**  
   Regions with coverage below thresholds will be filled with 'N' in consensus sequences.
   
**Quality Thresholds:**  
   Default quality thresholds can be modified in params.yaml as needed.

## 🤝 Contributing
We welcome contributions to make Juno better! Feel free to open issues or submit pull requests to suggest any additional features or enhancements!

## 📧 Contact
**Email**: bphl-sebioinformatics@flhealth.gov

## ⚖️ License
Juno is licensed under the [MIT License](LICENSE).

