#!/bin/bash
#SBATCH --job-name=vms55_occu_200km
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vms55@georgetown.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=15gb
#SBATCH --time=120:00:00

echo "Working directory is 'pwd'"
echo "Starting R at 'date'."

module load jags
module load R
export OMP_NUM_THREADS=1

Rscript scripts/004_RunModel_100.R

echo "Finished R at 'date'."


