#!/bin/bash
#
#
#SBATCH -J PreRecalHaplotypecaller
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH --cpus-per-task 20
#SBATCH -t 72:00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/PreRecalHaplotypecaller.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/PreRecalHaplotypecaller.err
#SBATCH --account=coxlab
#SBATCH --partition=standard


pre_processing_dir=""


cd ${pre_processing_dir}

## Modules may chance according to sys admin

module load gcc/11.4.0 openmpi/4.1.4
module load python/3.11.4
module load gatk/4.2.3.0

refgen=""

for i in *_prehaplo_sort.bam
do
	name=$(echo $i | awk -F "_prehaplo_sort" '{print $1}')
	gatk --java-options "-Xmx50g" HaplotypeCaller \
		-R ${refgen} \
		-I $i \
		-O ${name}.vcf \
		--standard-min-confidence-threshold-for-calling 30 \
		--min-base-quality-score 25 \
		-ERC GVCF
done
			
