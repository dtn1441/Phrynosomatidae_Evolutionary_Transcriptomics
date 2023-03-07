#!/bin/bash
#
#
#SBATCH -J Haplotypecaller
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH --cpus-per-task 20
#SBATCH -t 72:00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/haplotypecaller.%A_%a.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/haplotypecaller.%A_%a.err
#SBATCH --account=coxlab
#SBATCH --partition=standard
#SBATCH --array=1-11



cd /scratch/dtn2an/workspace/refseq_aligned/scja_3/

interval_file="/scratch/dtn2an/workspace/refseq_aligned/interval_list.txt"
myint=$(cat ${interval_file} | sed "${SLURM_ARRAY_TASK_ID}q;d")

module load gatk/4.2.3.0
refgen="/scratch/dtn2an/genomes/refseq_genome/ncbi-genomes-2021-09-15/GCF_019175285.1_SceUnd_v1.1_genomic.fna"

for i in *_prehaplotypecaller_rdgrps.bam
do
	name=$(echo $i | awk -F "_prehaplotypecaller" '{print $1}')
	gatk --java-options "-Xmx50g" HaplotypeCaller \
		-R ${refgen} \
		-I $i \
		-O ./${name}_${myint}.vcf \
		--standard-min-confidence-threshold-for-calling 30 \
		--min-base-quality-score 25 \
		-ERC GVCF
		-L ${interval_file}
done
			
