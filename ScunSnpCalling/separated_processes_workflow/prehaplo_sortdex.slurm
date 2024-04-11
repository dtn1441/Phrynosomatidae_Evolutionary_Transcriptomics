#!/bin/bash
#
#
#SBATCH -J initial_sortdex
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 0-20:00:00
#SBATCH --mem 25G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/sort_index.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/sort_index.err
#SBATCH --account=coxlab
#SBATCH --partition=standard

pre_processing_dir=""


cd ${pre_processing_dir}

module load samtools

for i in *_splitNcigars.bam
do
	name=$(echo $i | awk -F "_split" '{print $1}')
	samtools sort $i > ${name}_prehaplo_sort.bam
	samtools index ${name}_prehaplo_sort.bam
done
