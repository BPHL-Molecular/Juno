# Juno 🦟🦠🧬📊 - A Nextflow Pipeline for the Assembly of Oropouche Virus (OROV) Genomes

Juno is a flexible pipeline for Illumina data that supports both reference-based and de novo assembly of Oropouche virus genomes.

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
- [Slurm](https://slurm.schedmd.com/documentation.html) (This applies only if HiPerGator is used)

## ⚙️ Configuration

#### 1. Clone this repository

```bash
git clone https://github.com/BPHL-Molecular/Juno.git
cd Juno
```

#### 2. Create a directory for Input FASTQ Files

```bash
mkdir fastq
# move or copy your FASTQ files into this directory
```

#### 3. Set required parameters:
**Important:** All pipeline parameters **must be set in the `params.yaml` file**. Make sure you edit this file to provide the correct paths and values before running the pipeline. 

```yaml
# Input/Output paths
input_dir: "/path/to/fastq"
output_dir: "/path/to/juno_output"

# Kraken2 database path
kraken2_db: "/path/to/kraken2_db"

# Human read removal using NCBI's SRA human read removal tool (HRRT)
skip_hrrt: false

# Quality control thresholds
qc_thresholds:
    min_coverage: 90 # Minimum coverage percentage
    min_depth: 15 # Minimum depth of coverage (x)
```

You will need to download the kraken2/bracken viral database from the BenLangmead Index zone [link](https://benlangmead.github.io/aws-indexes/k2) for read classification.

The pipeline includes a step for removing human reads using [NCBI’s SRA Human Read Removal Tool (HRRT)](https://ncbiinsights.ncbi.nlm.nih.gov/2023/02/02/scrubbing-human-sequences-sra-submissions/). This step is enabled by default, but please note that it significantly increases runtime due to the large container size and the intensive I/O involved in decompressing input files and recompressing cleaned outputs.

To skip this step and use raw reads directly, set: ```--skip_hrrt: true```

Note: Skipping HRRT may be appropriate for:
- Non-human samples
- Pre-cleaned datasets
- Testing/development workflows

###### Please see the [notes](https://github.com/BPHL-Molecular/Juno/tree/main/references) on the reference genomes used in this pipeline.

## 🛠️ Pipeline Steps
1. **Quality Control**
   - Human Read Removal - [`sra-human-scrubber`](https://github.com/ncbi/sra-human-scrubber)
   - Read QC and trimming - [`fastp`](https://github.com/OpenGene/fastp)
2. **Taxonomic Classification**
   - Read classification - [`kraken2`](https://github.com/DerrickWood/kraken2)
   - Filter classified OROV reads - [`krakentools`](https://github.com/jenniferlu717/KrakenTools)
3. **Assembly**
   - Reference alignment - [`bwa`](https://github.com/lh3/bwa)
   - SAM/BAM processing - [`samtools`](https://github.com/samtools/samtools)
   - Variant calling & consensus - [`ivar`](https://github.com/andersen-lab/ivar)
4. **Quality Assessment**
   - Assembly evaluation - [`quast`](https://github.com/ablab/quast)
5. **Aggregate and Summarize Results**
   - Aggregate results from bioinformatics analyses - [`multiqc`](https://github.com/MultiQC/MultiQC)
   - Summary report generation - python-based script to summarize key assembly and QC metrics

## 📂 Output Structure
```
output_dir/
├── dehosted/         # Cleaned reads
├── trimmed/          # Trimmed reads
├── kraken2/          # Classification results
├── filtered_reads/   # Classified OROV reads
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
- Raw read counts
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
