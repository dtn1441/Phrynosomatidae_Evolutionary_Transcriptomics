#!/bin/bash
#
#
#SBATCH -J Haplotypecaller
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH --cpus-per-task 20
#SBATCH -t 72:00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/haplotypecaller.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/haplotypecaller.err
#SBATCH --account=coxlab
#SBATCH --partition=standard

cd /scratch/dtn2an/workspace/refseq_aligned/merge_runs/new_tissue_merged/

module load gatk/4.2.3.0
refgen="/scratch/dtn2an/genomes/refseq_genome/ncbi-genomes-2021-09-15/GCF_019175285.1_SceUnd_v1.1_genomic.fna"

for i in *_prehaplotypecaller_sort.bam
do
	name=$(echo $i | awk -F "_prehaplotypecaller" '{print $1}')
	/scratch/dtn2an/tools/samtools-1.16.1/samtools samples $i | column | awk -F '\t' '{for (i = 1; i <= NF; i++) {if ($i ~ /Scun/) print $i} }' > temp.txt
	samples=$(cat temp.txt)
	echo ${samples}
	for sample in $samples
	do
		echo ${sample} 
		gatk --java-options "-Xmx50g" HaplotypeCaller \
			-R ${refgen} \
			-I $i \
			-O ./variant_calls/${name}_${sample}.vcf \
			--sample-name ${sample}
			--standard-min-confidence-threshold-for-calling 30 \
			--min-base-quality-score 25 \
			-ERC GVCF
	done
done
			
