#!/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=juno_original_reference
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arnold.rodriguezhilario@flhealth.gov
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100gb
#SBATCH --time=5:00:00
#SBATCH --output=juno_original_reference.%j.out
#SBATCH --error=juno_original_reference.%j.err

# Load required modules
module load apptainer nextflow python3

# Set up singularity/apptainer cache
export NXF_SINGULARITY_CACHEDIR=/blue/bphl-florida/arnold.rodriguez/singularity
export PATH

# Run nextflow command
nextflow run juno.nf -profile singularity -params-file params.yaml

