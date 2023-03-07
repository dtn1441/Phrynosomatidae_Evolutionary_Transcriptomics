#!/bin/bash
#
#
#SBATCH -J FinalSortAndIndex_scun_cat
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 0-20:00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/sortdex__scun_cat.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/sortdex_scun_cat.err
#SBATCH --account=coxlab
#SBATCH --partition=standard

pre_processing_dir=""



cd ${pre_processing_dir}

module load samtools/1.12

for i in *_splitNcigars.bam
do
	name=$(echo $i | awk -F "_split" '{print $1}')
	samtools sort $i > ${name}_prehaplotypecaller_sort.bam
	samtools index ${name}_prehaplotypecaller_sort.bam
done
