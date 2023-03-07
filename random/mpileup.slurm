#!/bin/bash
#
#
#SBATCH -J mpileup
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 100:00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/mpileup.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/mpileup.err
#SBATCH --account=coxlab
#SBATCH --partition=standard

cd /scratch/dtn2an/workspace/refseq_aligned/merge_runs/mpileups/

module load samtools/1.12
for i in *_prehaplotypecaller_sort.bam
do 
	name=$(echo $i | awk -F _ '{print $1}')
	samtools mpileup -f /scratch/dtn2an/genomes/refseq_genome/ncbi-genomes-2021-09-15/GCF_019175285.1_SceUnd_v1.1_genomic.fna $i > ${name}.mpileup 
	echo ${name}
done



