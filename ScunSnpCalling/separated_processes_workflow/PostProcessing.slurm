#!/bin/bash
#
#
#SBATCH -J PostProcessing
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=120gb
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/PostProcessing.%A_%a.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/PostProcessing.%A_%a.err
#SBATCH --account=coxlab
#SBATCH --partition=standard
#SBATCH --array=1-9000%100

joint_calls_dir=""
refgen=""

##modules may change according to sys admin
module load gcc/11.4.0 openmpi/4.1.4
module load python/3.11.4
module load gatk/4.2.3.0


cd ${joint_calls_dir}




##Declare intervals for array
interval_list=""
myint=$(cat ${interval_list} | sed "${SLURM_ARRAY_TASK_ID}q;d")

##Separate snps and indels


cd ${joint_calls_dir}
mkdir PostProcessing

gatk SelectVariants \
  -V ${myint}_jc.vcf \
  --select-type-to-include SNP \
  --output PostProcessing/${myint}_snps.vcf 
  
gatk SelectVariants \
  -V ${myint}_jc.vcf \
  --select-type-to-include INDEL \
  --output PostProcessing/${myint}_indels.vcf 


## Filter variants based off recommended process for RNAseq data

cd PostProcessing

gatk VariantFiltration \
  -R ${refgen} \
  -V ${myint}_snps.vcf \
  --output ${myint}_snps_marked.vcf \
  --filter-expression "DP < 10.0 || FS > 30.0 || QUAL < 30.0 || MQ < 30.0 || QD < 2.0" --filter-name "DP10.FS30.QUAL30.MQ30.QD2"

gatk SelectVariants \
  -R ${refgen} \
  -V ${myint}_snps_marked.vcf \
  --exclude-filtered true \
  --output ${myint}_snps_filtered.vcf
  
gatk VariantFiltration \
  -R ${refgen} \
  -V ${myint}_indels.vcf \
  --output ${myint}_indels_marked.vcf \
  --filter-expression "FS > 200.0 || QUAL < 30.0 || ReadPosRankSum > -20.0 || QD < 2.0" --filter-name "FS200.QUAL30.RPRS-20.QD2"

gatk SelectVariants \
  -R ${refgen} \
  -V ${myint}_indels_marked.vcf \
  --exclude-filtered true \
  --output ${myint}_indels_filtered.vcf
  
rm *marked*

echo "Finished Post-processing"
