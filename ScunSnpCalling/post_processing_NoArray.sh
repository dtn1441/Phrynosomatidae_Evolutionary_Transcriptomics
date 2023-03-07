#!/bin/bash
#
#
#SBATCH -J scja_post_processing
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=120gb
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/scja_post_processing.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/scja_post_processing.err
#SBATCH --account=coxlab
#SBATCH --partition=standard


refgen="/scratch/dtn2an/genomes/refseq_genome/ncbi-genomes-2021-09-15/GCF_019175285.1_SceUnd_v1.1_genomic.fna"


##Separate snps and indels

cd /scratch/dtn2an/workspace/refseq_aligned/scja_snps/

module load gatk

for i in *.vcf
do
	name=$(echo $i | awk -F "_jc" '{print $1}')
		gatk SelectVariants \
	  		-V $i \
 	 		--select-type-to-include SNP \
 	 		--output ${name}_snps.vcf 
 		gatk SelectVariants \
  			-V $i \
  			--select-type-to-include INDEL \
  			--output ${name}_indels.vcf 
		gatk VariantFiltration \
  			-R ${refgen} \
  			-V ${name}_snps.vcf \
  			--output ${name}_snps_marked.vcf \
  			--filter-expression "DP < 10.0 || FS > 30.0 || QUAL < 30.0 || MQ < 30.0 || QD < 2.0" --filter-name "DP10.FS30.QUAL30.MQ30.QD2"
		gatk SelectVariants \
  			-R ${refgen} \
  			-V ${name}_snps_marked.vcf \
  			--exclude-filtered true \
  			--output ${name}_snps_filtered.vcf  
		gatk VariantFiltration \
  			-R ${refgen} \
  			-V ${name}_indels.vcf \
  			--output ${name}_indels_marked.vcf \
  			--filter-expression "FS > 200.0 || QUAL < 30.0 || ReadPosRankSum > -20.0 || QD < 2.0" --filter-name "FS200.QUAL30.RPRS-20.QD2"
		gatk SelectVariants \
  			-R ${refgen} \
  			-V ${name}_indels_marked.vcf \
			--exclude-filtered true \
  			--output ${name}_indels_filtered.vcf
done

  
rm *marked*

echo "Finished Post-processing"
