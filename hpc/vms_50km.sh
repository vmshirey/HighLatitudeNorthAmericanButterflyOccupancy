#!/bin/bash
#SBATCH --job-name=vms55_occu_50km
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vms55@georgetown.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=30
#SBATCH --time=140:00:00

echo "Start  Date              = $(date)"

echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"

module load jags
module load R
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

Rscript scripts/004_RunModel_50.R

echo "Finish Date              = $(date)"


