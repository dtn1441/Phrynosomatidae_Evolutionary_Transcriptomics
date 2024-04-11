#!/bin/bash
#
#
#SBATCH -J MarkRemoveDuplicates
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 0-15:00:00
#SBATCH --mem 20G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/MarkRemoveDuplicates.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/MarkRemoveDuplicates.err
#SBATCH --account=coxlab
#SBATCH --partition=standard

pre_processing_dir=""

##following modules may change according to sys admin

module load gcc/11.4.0 openmpi/4.1.4
module load python/3.11.4
module load gatk/4.2.3.0

for i in *_initial_sort.bam
do
	name=$(echo $i | awk -F "_initial" '{print $1}')
	gatk MarkDuplicates \
		-I $i \
		-O ${name}_removed_dups.bam \
		-M marked_dups_metrics.txt \
		-REMOVE_DUPLICATES true \
		-READ_NAME_REGEX null
done

