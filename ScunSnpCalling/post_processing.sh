#!/bin/bash
#
#
#SBATCH -J scun_snp_post_processing_scun_cat
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=120gb
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/scun_snp_post_processing_scun_cat.%A_%a.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/scun_snp_post_processing_scun_cat.%A_%a.err
#SBATCH --account=coxlab
#SBATCH --partition=standard
#SBATCH --array=1-9000%100


PrimDir="/scratch/dtn2an/workspace/refseq_aligned"
SecDir="/scratch/dtn2an/workspace/refseq_aligned/scun_cat/variant_calls"
ThirDir="/scratch/dtn2an/workspace/refseq_aligned/scun_cat/variant_calls/Genotype_gVCFs"
refgen="/scratch/dtn2an/genomes/refseq_genome/ncbi-genomes-2021-09-15/GCF_019175285.1_SceUnd_v1.1_genomic.fna"

mkdir ${SecDir}/post_processing/

##Declare intervals for array
interval_list="/scratch/dtn2an/workspace/refseq_aligned/interval_list.txt"
myint=$(cat ${interval_list} | sed "${SLURM_ARRAY_TASK_ID}q;d")

##Separate snps and indels

cd ${ThirDir}

module load gatk

gatk SelectVariants \
  -V ${myint}_jc.vcf \
  --select-type-to-include SNP \
  --output ${SecDir}/post_processing/${myint}_scun_cat_snps.vcf 
  
gatk SelectVariants \
  -V ${myint}_jc.vcf \
  --select-type-to-include INDEL \
  --output ${SecDir}/post_processing/${myint}_scun_cat_indels.vcf 

#rm ${ThirDir}/*_jc.vcf*


## Filter variants based off recommended process for RNAseq data

cd ${SecDir}/post_processing/

gatk VariantFiltration \
  -R ${refgen} \
  -V ${myint}_scun_cat_snps.vcf \
  --output ${myint}_scun_cat_snps_marked.vcf \
  --filter-expression "DP < 10.0 || FS > 30.0 || QUAL < 30.0 || MQ < 30.0 || QD < 2.0" --filter-name "DP10.FS30.QUAL30.MQ30.QD2"

gatk SelectVariants \
  -R ${refgen} \
  -V ${myint}_scun_cat_snps_marked.vcf \
  --exclude-filtered true \
  --output ${myint}_scun_cat_snps_filtered.vcf
  
gatk VariantFiltration \
  -R ${refgen} \
  -V ${myint}_scun_cat_indels.vcf \
  --output ${myint}_scun_cat_indels_marked.vcf \
  --filter-expression "FS > 200.0 || QUAL < 30.0 || ReadPosRankSum > -20.0 || QD < 2.0" --filter-name "FS200.QUAL30.RPRS-20.QD2"

gatk SelectVariants \
  -R ${refgen} \
  -V ${myint}_scun_cat_indels_marked.vcf \
  --exclude-filtered true \
  --output ${myint}_scun_cat_indels_filtered.vcf
  
##rm *marked*

echo "Finished Post-processing"
