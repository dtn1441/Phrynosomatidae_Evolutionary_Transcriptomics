#!/bin/bash
#
#
#SBATCH -J PostRecalHaplotypecaller
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH --cpus-per-task 20
#SBATCH -t 72:00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/PostRecalHaplotypecaller.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/PostRecalHaplotypecaller.err
#SBATCH --account=coxlab
#SBATCH --partition=standard


bqsr_dir=""
ref_fasta=""

cd ${bqsr_dir}

## Modules may chance according to sys admin

module load gcc/11.4.0 openmpi/4.1.4
module load python/3.11.4
module load gatk/4.2.3.0



for i in *_postRecalPreHaploSort.bam
do
	name=$(echo $i | awk -F "_post" '{print $1}')
	gatk --java-options "-Xmx50g" HaplotypeCaller \
		-R ${ref_fasta} \
		-I $i \
		-O ${name}.vcf \
		--standard-min-confidence-threshold-for-calling 30 \
		--min-base-quality-score 25 \
		-ERC GVCF
done


