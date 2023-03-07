#!/bin/bash
#
#
#SBATCH -J scun_cat_dup_removal
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 0-40:00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/scun_cat_dup_removal.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/scun_cat_dup_removal.err
#SBATCH --account=coxlab
#SBATCH --partition=standard

cd /scratch/dtn2an/workspace/refseq_aligned/scun_cat/

module load samtools/1.12

for i in *.bam
do
	name=$(echo $i | awk -F . '{print $1}')
	samtools sort $i > ${name}_sorted.bam
done 


module load gatk/4.2.3.0

for i in *_sorted.bam
do
	name=$(echo $i | awk -F "_sorted" '{print $1}')
	gatk MarkDuplicates \
		-I $i \
		-O ${name}_removed_dups.bam \
		-M marked_dups_metrics.txt \
		-REMOVE_DUPLICATES true \
		-READ_NAME_REGEX null
done

echo finished
