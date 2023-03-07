#!/bin/bash
#
#
#SBATCH -J PoolSNP
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 0-40:00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/PoolSNP_%A_%a.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/PoolSNP_%A_%a.err
#SBATCH --account=coxlab
#SBATCH --partition=standard
#SBATCH --array=1-380

module load parallel/20200322

cd /scratch/dtn2an/workspace/refseq_aligned/merge_runs/mpileups/finished/chrom/

mkdir PoolSNP

input_files=/scratch/dtn2an/workspace/refseq_aligned/input_files.txt

name_chrom=$(cat ${input_files} | awk '{print $1}' | sed "${SLURM_ARRAY_TASK_ID}q;d")
name_only=$(cat ${input_files} | awk '{print $2}' | sed "${SLURM_ARRAY_TASK_ID}q;d")
echo $name_chrom
echo $name_only
echo $output_name

bash /scratch/dtn2an/tools/PoolSNP/PoolSNP.sh \
mpileup=/scratch/dtn2an/workspace/refseq_aligned/merge_runs/mpileups/finished/chrom/${name_chrom}.mpileup \
reference=/scratch/dtn2an/genomes/refseq_genome/ncbi-genomes-2021-09-15/GCF_019175285.1_SceUnd_v1.1_genomic.fna \
names=${name_only} \
max-cov=0.98 \
min-cov=10 \
min-count=10 \
min-freq=0.01 \
miss-frac=0.2 \
badsites=1 \
allsites=1 \
output=/scratch/dtn2an/workspace/refseq_aligned/merge_runs/mpileups/finished/chrom/PoolSNP/${name_chrom}_PoolSNP.vcf

echo ${output_name}
