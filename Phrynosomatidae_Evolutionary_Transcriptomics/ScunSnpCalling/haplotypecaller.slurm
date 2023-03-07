#!/bin/bash
#
#
#SBATCH -J Haplotypecaller_scun_cat
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH --cpus-per-task 20
#SBATCH -t 72:00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/haplotypecaller_scun_cat_gr1.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/haplotypecaller_scun_cat_gr1.err
#SBATCH --account=coxlab
#SBATCH --partition=standard


pre_processing_dir=""


cd ${pre_processing_dir}

module load gatk/4.2.3.0
refgen="/scratch/dtn2an/genomes/refseq_genome/ncbi-genomes-2021-09-15/GCF_019175285.1_SceUnd_v1.1_genomic.fna"

for i in *_prehaplotypecaller_RG.bam
do
	name=$(echo $i | awk -F "_prehaplotypecaller" '{print $1}')
	gatk --java-options "-Xmx50g" HaplotypeCaller \
		-R ${refgen} \
		-I $i \
		-O ${name}.vcf \
		--standard-min-confidence-threshold-for-calling 30 \
		--min-base-quality-score 25 \
		-ERC GVCF
done
			
