#!/bin/bash
#
#
#SBATCH -J scun_joint_snp_calling_scun_cat
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=120gb
#SBATCH --output /scratch/dtn2an/workspace/refseq_aligned/stderr/joint_snp_calls_scun_cat.%A_%a.out
#SBATCH --error /scratch/dtn2an/workspace/refseq_aligned/stderr/joint_snp_calls_scun_cat.%A_%a.err
#SBATCH --account=coxlab
#SBATCH --partition=standard
#SBATCH --array=1-9000%100


PrimDir="/scratch/dtn2an/workspace/refseq_aligned"
SecDir="/scratch/dtn2an/workspace/refseq_aligned/scun_cat/variant_calls/"
refgen="/scratch/dtn2an/genomes/refseq_genome/ncbi-genomes-2021-09-15/GCF_019175285.1_SceUnd_v1.1_genomic.fna"

## Now we can run GenomicsDBImport which will make a database of all of our reference called variants and the reference to aid in joint calling. We will then immediately joint call and then delete the genomics database since it takes up so many files

interval_file="/scratch/dtn2an/workspace/refseq_aligned/interval_list.txt"
myint=$(cat ${interval_file} | sed "${SLURM_ARRAY_TASK_ID}q;d")

cd ${SecDir}
module load gatk

gatk GenomicsDBImport --batch-size 7 \
  --reader-threads 10 \
  --variant AF3.vcf \
  --variant AF4.vcf \
  --variant AM5.vcf \
  --variant AM6.vcf \
  --variant AF7.vcf \
  --variant AM8.vcf \
  --variant AM9.vcf \
  --variant AF10.vcf \
  --variant AF11.vcf \
  --variant AM12.vcf \
  --variant AF13.vcf \
  --variant AM14.vcf \
  --variant AM15.vcf \
  --variant YM21.vcf \
  --variant YM22.vcf \
  --variant YM23.vcf \
  --variant YM24.vcf \
  --variant YM25.vcf \
  --variant YF26.vcf \
  --variant YF27.vcf \
  --variant YM28.vcf \
  --variant YF29.vcf \
  --variant YF31.vcf \
  --variant YF32.vcf \
  --variant YF33.vcf \
  --variant HF40.vcf \
  --variant HF41.vcf \
  --variant HF42.vcf \
  --variant HM43.vcf \
  --variant HF44.vcf \
  --variant HM45.vcf \
  --variant HM46.vcf \
  --variant HM47.vcf \
  --variant HM48.vcf \
  --variant HF49.vcf \
  --variant HM50.vcf \
  --variant HF51.vcf \
  --genomicsdb-workspace-path ./${myint}_wd \
  --intervals $myint
  
mkdir ${SecDir}/Genotype_gVCFs

gatk GenotypeGVCFs \
  -R ${refgen} \
  -V gendb://./${myint}_wd \
  -O ${SecDir}/Genotype_gVCFs/${myint}_jc.vcf

rm -r ./${myint}_wd


