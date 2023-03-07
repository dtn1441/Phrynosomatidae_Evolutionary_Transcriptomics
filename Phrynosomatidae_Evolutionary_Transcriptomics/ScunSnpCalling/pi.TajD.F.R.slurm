#!/bin/bash
#
#
#SBATCH -J pi.TajD.F.R
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 0-20:00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/pi.TajD.F.R.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/pi.TajD.F.R.err
#SBATCH --account=coxlab
#SBATCH --partition=standard

module load goolf/7.1.0_3.1.4 R

Rscript /scratch/dtn2an/workspace/refseq_aligned/commands/pi.TajD.F.R
