#!/bin/bash
#
#
#SBATCH -J scja_joint_snp_calling
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=120gb
#SBATCH --output /scratch/dtn2an/workspace/refseq_aligned/stderr/scja_joint_snp_calls.out
#SBATCH --error /scratch/dtn2an/workspace/refseq_aligned/stderr/scja_joint_snp_calls.err
#SBATCH --account=coxlab
#SBATCH --partition=standard


PrimDir="/scratch/dtn2an/workspace/refseq_aligned"
SecDir="/scratch/dtn2an/workspace/refseq_aligned/haplotypecaller_scja"
refgen="/scratch/dtn2an/genomes/refseq_genome/ncbi-genomes-2021-09-15/GCF_019175285.1_SceUnd_v1.1_genomic.fna"

## Now we can run GenomicsDBImport which will make a database of all of our reference called variants and the reference to aid in joint calling. We will then immediately joint call and then delete the genomics database since it takes up so many files

cd ${SecDir}/chrom1/
module load gatk

gatk GenomicsDBImport --batch-size 7 \
  --reader-threads 10 \
  --sample-name-map ${SecDir}/chrom1/map_file_chrom1.txt \
  -L NC_056522.1 \
  --genomicsdb-workspace-path ./chrom1_wd 
  
mkdir ${SecDir}/chrom1/Genotype_gVCFs

gatk GenotypeGVCFs \
  -R ${refgen} \
  -V gendb://./chrom1_wd \
  -O ${SecDir}/chrom1/Genotype_gVCFs/chrom1_jc.vcf

rm -r ./chrom1_wd


cd ${SecDir}/chrom2/
module load gatk

gatk GenomicsDBImport --batch-size 7 \
  --reader-threads 10 \
  --sample-name-map ${SecDir}/chrom2/map_file_chrom2.txt \
  -L NC_056523.1 \
  --genomicsdb-workspace-path ./chrom2_wd

mkdir ${SecDir}/chrom2/Genotype_gVCFs

gatk GenotypeGVCFs \
  -R ${refgen} \
  -V gendb://./chrom2_wd \
  -O ${SecDir}/chrom2/Genotype_gVCFs/chrom2_jc.vcf

rm -r ./chrom2_wd



cd ${SecDir}/chrom3/
module load gatk

gatk GenomicsDBImport --batch-size 7 \
  --reader-threads 10 \
  --sample-name-map ${SecDir}/chrom3/map_file_chrom3.txt \
  -L NC_056524.1 \
  --genomicsdb-workspace-path ./chrom3_wd

mkdir ${SecDir}/chrom3/Genotype_gVCFs

gatk GenotypeGVCFs \
  -R ${refgen} \
  -V gendb://./chrom3_wd \
  -O ${SecDir}/chrom3/Genotype_gVCFs/chrom3_jc.vcf

rm -r ./chrom3_wd



cd ${SecDir}/chrom4/
module load gatk

gatk GenomicsDBImport --batch-size 7 \
  --reader-threads 10 \
  --sample-name-map ${SecDir}/chrom4/map_file_chrom4.txt \
  -L NC_056525.1 \
  --genomicsdb-workspace-path ./chrom4_wd

mkdir ${SecDir}/chrom4/Genotype_gVCFs

gatk GenotypeGVCFs \
  -R ${refgen} \
  -V gendb://./chrom4_wd \
  -O ${SecDir}/chrom4/Genotype_gVCFs/chrom4_jc.vcf

rm -r ./chrom4_wd




cd ${SecDir}/chrom5/
module load gatk

gatk GenomicsDBImport --batch-size 7 \
  --reader-threads 10 \
  --sample-name-map ${SecDir}/chrom5/map_file_chrom5.txt \
  -L NC_056526.1 \
  --genomicsdb-workspace-path ./chrom5_wd

mkdir ${SecDir}/chrom5/Genotype_gVCFs

gatk GenotypeGVCFs \
  -R ${refgen} \
  -V gendb://./chrom5_wd \
  -O ${SecDir}/chrom5/Genotype_gVCFs/chrom5_jc.vcf

rm -r ./chrom5_wd









cd ${SecDir}/chrom6/
module load gatk

gatk GenomicsDBImport --batch-size 7 \
  --reader-threads 10 \
  --sample-name-map ${SecDir}/chrom6/map_file_chrom6.txt \
  -L NC_056527.1 \
  --genomicsdb-workspace-path ./chrom6_wd

mkdir ${SecDir}/chrom6/Genotype_gVCFs

gatk GenotypeGVCFs \
  -R ${refgen} \
  -V gendb://./chrom6_wd \
  -O ${SecDir}/chrom6/Genotype_gVCFs/chrom6_jc.vcf

rm -r ./chrom6_wd



cd ${SecDir}/chrom7/
module load gatk

gatk GenomicsDBImport --batch-size 7 \
  --reader-threads 10 \
  --sample-name-map ${SecDir}/chrom7/map_file_chrom7.txt \
  -L NC_056528.1 \
  --genomicsdb-workspace-path ./chrom7_wd

mkdir ${SecDir}/chrom7/Genotype_gVCFs

gatk GenotypeGVCFs \
  -R ${refgen} \
  -V gendb://./chrom7_wd \
  -O ${SecDir}/chrom7/Genotype_gVCFs/chrom7_jc.vcf

rm -r ./chrom7_wd







cd ${SecDir}/chrom8/
module load gatk

gatk GenomicsDBImport --batch-size 7 \
  --reader-threads 10 \
  --sample-name-map ${SecDir}/chrom8/map_file_chrom8.txt \
  -L NC_056529.1 \
  --genomicsdb-workspace-path ./chrom8_wd

mkdir ${SecDir}/chrom8/Genotype_gVCFs

gatk GenotypeGVCFs \
  -R ${refgen} \
  -V gendb://./chrom8_wd \
  -O ${SecDir}/chrom8/Genotype_gVCFs/chrom8_jc.vcf

rm -r ./chrom8_wd




cd ${SecDir}/chrom9/
module load gatk

gatk GenomicsDBImport --batch-size 7 \
  --reader-threads 10 \
  --sample-name-map ${SecDir}/chrom9/map_file_chrom9.txt \
  -L NC_056530.1 \
  --genomicsdb-workspace-path ./chrom9_wd

mkdir ${SecDir}/chrom9/Genotype_gVCFs

gatk GenotypeGVCFs \
  -R ${refgen} \
  -V gendb://./chrom9_wd \
  -O ${SecDir}/chrom9/Genotype_gVCFs/chrom9_jc.vcf

rm -r ./chrom9_wd




cd ${SecDir}/chrom10/
module load gatk

gatk GenomicsDBImport --batch-size 7 \
  --reader-threads 10 \
  --sample-name-map ${SecDir}/chrom10/map_file_chrom10.txt \
  -L NC_056531.1 \
  --genomicsdb-workspace-path ./chrom10_wd

mkdir ${SecDir}/chrom10/Genotype_gVCFs

gatk GenotypeGVCFs \
  -R ${refgen} \
  -V gendb://./chrom10_wd \
  -O ${SecDir}/chrom10/Genotype_gVCFs/chrom10_jc.vcf

rm -r ./chrom10_wd




cd ${SecDir}/chrom11/
module load gatk

gatk GenomicsDBImport --batch-size 7 \
  --reader-threads 10 \
  --sample-name-map ${SecDir}/chrom11/map_file_chrom11.txt \
  -L NC_056532.1 \
  --genomicsdb-workspace-path ./chrom11_wd

mkdir ${SecDir}/chrom11/Genotype_gVCFs

gatk GenotypeGVCFs \
  -R ${refgen} \
  -V gendb://./chrom11_wd \
  -O ${SecDir}/chrom11/Genotype_gVCFs/chrom11_jc.vcf

rm -r ./chrom11_wd
