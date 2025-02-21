# Juno 🦟🦠🧬📊 - A Nextflow Pipeline for Reference-Based Assembly of Oropouche Virus (OROV) Genomes

Juno is designed for processing Illumina paired-end metagenomics sequencing data against OROV reference genomes, performing QC, taxonomic classification, alignment, variant calling, and consensus generation.

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

**Important:** All pipeline parameters **must be set in the `params.yaml` file**. Make sure you edit this file to provide the correct paths and values before running the pipeline. 

You will also need to download the kraken2 viral database from the BenLangmead Index zone [link](https://benlangmead.github.io/aws-indexes/k2).

#### Required Parameters (in `params.yaml`):

```yaml
# Input/Output paths
input_dir: "/path/to/fastq"
output_dir: "/path/to/output_dir"

# References path, default reference directory, DO NOT change.
refs_dir: "${projectDir}/references"

# Database path
kraken2_db: "/path/to/kraken2_db"

# Resource configuration, default number of threads per process
threads: 32

# Human scrubber processing option, set to true for HPC environments
parallel_hrrt: false

# Quality control thresholds
qc_thresholds:
    min_coverage: 90
    min_depth: 15
```
###### Please see the [notes](https://github.com/BPHL-Molecular/Juno/tree/main/references) on the references sequences used in this pipeline.

## 🛠️ Pipeline Steps
1. **Quality Control**
   - Human Read Removal (`sra-human-scrubber`)
   - Read QC and trimming (`fastp`)
2. **Taxonomic Classification**
   - Read classification (`kraken2`)
3. **Assembly**
   - Reference alignment (`bwa`)
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
├── alignments/       # SAM/BAM files & indices
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
- Classification read counts
- Mapping statistics
- Coverage metrics
- Variant counts
- Assembly quality metrics
- Overall QC status

## 🐛 Troubleshooting
**Pipeline Errors:**  
   Check Nextflow execution logs in .nextflow.log 
   
**Low Coverage Regions:**  
   Regions with low coverage (<10x) will be filled with 'N' in consensus sequences.
   
**Quality Thresholds:**  
   Default quality thresholds can be modified in params.yaml as needed.

## 🤝 Contributing
We welcome contributions to make Juno better! Feel free to open issues or submit pull requests to suggest any additional features or enhancements!

## 📧 Contact
**Email**: bphl-sebioinformatics@flhealth.gov

## ⚖️ License
Juno is licensed under the [MIT License](LICENSE).

