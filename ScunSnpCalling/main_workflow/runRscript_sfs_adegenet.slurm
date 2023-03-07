#!/bin/bash
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/sfs.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/sfs.err
#SBATCH -p standard
#SBATCH -A coxlab






module load goolf/7.1.0_3.1.4 R

Rscript /scratch/dtn2an/workspace/refseq_aligned/commands/scun_site_frequency_spectrum.R
