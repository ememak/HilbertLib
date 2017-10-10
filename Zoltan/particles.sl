#!/bin/sh
#SBATCH --job-name=ZoltanLib
#SBATCH --nodes=4
#SBATCH --tasks-per-node=48
#SBATCH --time=30:00
#SBATCH --account=icm-praktyki
#SBATCH --output=./results/MainOutput.out
#SBATCH --error=./results/MainOutput.err