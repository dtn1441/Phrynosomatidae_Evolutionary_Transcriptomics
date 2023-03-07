#!/bin/bash
#
#
#SBATCH -J bam2vcfTest
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 0-100:00:00
#SBATCH --mem 100G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/bam2vcfTest.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/bam2vcfTest.err
#SBATCH --account=coxlab
#SBATCH --partition=standard
#SBATCH --array=1-9000%100

cd /scratch/dtn2an/workspace/refseq_aligned/commands/

./bam2vcf.sh
