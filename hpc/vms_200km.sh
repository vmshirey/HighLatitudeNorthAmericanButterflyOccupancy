#!/bin/bash
#SBATCH --job-name=vms55_occu_200km
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vms55@georgetown.edu
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --cpus-per-task=12
#SBATCH --time=120:00:00

echo "Start  Date              = $(date)"

echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"

module load jags
module load R
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

Rscript scripts/004_RunModel_200.R

echo "Finish Date              = $(date)"
