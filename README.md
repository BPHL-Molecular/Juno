# Juno 🦟🦠🧬📊
**A Nextflow Pipeline for Reference-Based Assembly of Oropouche Virus (OROV)**

A containerized Nextflow pipeline designed for processing Illumina paired-end metagenomics sequencing data against OROV reference genomes, performing QC, taxonomic classification, variant calling, and consensus generation.

## ⚡ Usage
```bash
$ nextflow run juno.nf -profile singularity -params-file params.yaml
```

## 📦 Dependencies and Required Parameters
- [Nextflow 23.04.0+](https://www.nextflow.io/docs/latest/install.html)
- [Docker](https://docs.docker.com/engine/install/]/(Singularity)[https://docs.sylabs.io/guides/latest/user-guide/quick_start.html#quick-installation-steps)
- [Python 3.6+](https://docs.python.org/3/using/unix.html)
- [Kraken2 viral database](https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20241228.tar.gz)


--input: Directory containing paired-end FASTQ files (must follow the _R1_ and _R2_ naming convention).

--output: Directory where all pipeline results will be stored.

--k2-db: Path to the Kraken2 database.

--threads: Number of CPU threads to use.


***It is highly recommended to remove human reads before running this pipeline. We include a python script in this repository for this purpose (host_removal.py). Although, this script will only work on HiPerGator.***


## 🛠️ Pipeline Steps
1. **Quality Control**
   - Raw read QC (`fastqc`)
   - Read trimming (`fastp`)
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
├── fastqc/           # Raw read QC reports
├── trimmed/          # Cleaned reads
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
- Raw and cleaned read counts
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

