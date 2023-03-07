#!/bin/bash
#SBATCH -n 1
#SBATCH --ntasks-per-node=24
#SBATCH -t 10:00:00
#SBATCH --mem 100G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/join_sex_AllFreq_positions.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/join_sex_AllFreq_positions.err
#SBATCH -p standard
#SBATCH -A coxlab






module load goolf/7.1.0_3.1.4 R

Rscript /scratch/dtn2an/workspace/refseq_aligned/commands/join_sex_AllFreq_positions.R
