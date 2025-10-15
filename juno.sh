#!/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=juno_orov_assembly
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=user@email.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100gb
#SBATCH --time=2:00:00
#SBATCH --output=juno_orov_assembly.%j.out
#SBATCH --error=juno_orov_assembly.%j.err

# Load required modules
module load apptainer nextflow python3

# Set up singularity/apptainer cache
export NXF_SINGULARITY_CACHEDIR=/path/to/singularity/cache
export PATH

# Run nextflow command
nextflow run juno.nf -profile singularity -params-file params.yaml
