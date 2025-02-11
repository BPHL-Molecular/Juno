#!/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=juno_orov_assembl
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=user@email.gov
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200gb
#SBATCH --time=5:00:00
#SBATCH --output=juno_orov_assembly.%j.out
#SBATCH --error=juno_orov_assembly.%j.err

# Load required modules
module load apptainer nextflow python3

# Run nextflow command
nextflow run juno.nf -profile singularity -params-file params.yaml

