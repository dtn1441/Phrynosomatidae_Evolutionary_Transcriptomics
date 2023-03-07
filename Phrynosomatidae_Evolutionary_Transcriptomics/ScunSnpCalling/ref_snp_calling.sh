#!/bin/bash
#
#
#SBATCH -J scun_ref_snp_calling
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=80gb
#SBATCH --output /scratch/dtn2an/refseq_aligned/stderr/ref_snp_calls.out
#SBATCH --error /scratch/dtn2an/refseq_aligned/stderr/ref_snp_calls.err
#SBATCH --account=coxlab
#SBATCH --partition=standard


##mkdir /scratch/dtn2an/refseq_aligned/snp_calling/

PrimDir="/scratch/dtn2an/refseq_aligned"
SecDir="/scratch/dtn2an/refseq_aligned/snp_calling"
refgen="/scratch/dtn2an/refseq_genome/ncbi-genomes-2021-09-15/GCF_019175285.1_SceUnd_v1.1_genomic.fna"

##mkdir ${SecDir}/liver
##mkdir ${SecDir}/brain
##mkdir ${SecDir}/muscle

##cd ${SecDir}/liver

##mkdir ${SecDir}/liver/haplotypecaller_output/

##module load gatk/4.1.6.0

##for i in *_prehaplotypecaller_sort.bam
##do
  ##name=$(echo $i | awk -F "_prehaplotypecaller_" '{print $1}')
  ##gatk --java-options "-Xmx25g" HaplotypeCaller \
    ##-R ${refgen} \
    ##-I $i \
    ##-O ${SecDir}/liver/haplotypecaller_output/${name}_rd1.vcf \
    ##--standard-min-confidence-threshold-for-calling 30 \
    ##--min-base-quality-score 25 \
    ##-ERC GVCF
##done

##cd ${SecDir}/brain

##mkdir ${SecDir}/brain/haplotypecaller_output/

##module load gatk/4.1.6.0

##for i in *_prehaplotypecaller_sort.bam
##do
##  name=$(echo $i | awk -F "_prehaplotypecaller_" '{print $1}')
##  gatk --java-options "-Xmx25g" HaplotypeCaller \
##    -R ${refgen} \
##    -I $i \
##    -O ${SecDir}/brain/haplotypecaller_output/${name}_rd1.vcf \
##    --standard-min-confidence-threshold-for-calling 30 \
##    --min-base-quality-score 25 \
##    -ERC GVCF
##done

cd ${SecDir}/muscle

mkdir ${SecDir}/muscle/haplotypecaller_output/

module load gatk/4.1.6.0

for i in *_prehaplotypecaller_sort.bam
do
  name=$(echo $i | awk -F "_prehaplotypecaller_" '{print $1}')
  gatk --java-options "-Xmx25g" HaplotypeCaller \
    -R ${refgen} \
    -I $i \
    -O ${SecDir}/muscle/haplotypecaller_output/${name}_rd1.vcf \
    --standard-min-confidence-threshold-for-calling 30 \
    --min-base-quality-score 25 \
    -ERC GVCF
done

cd ${PrimDir}

grep ">" ${refgen} | head -9000 > interval_list_header_too_long.txt
cat interval_list_header_too_long.txt | awk -F " " '{print $1}' > interval_list_headers_w_carrot.txt
cat interval_list_headers_w_carrot.txt | sed 's/^.//' > interval_list.txt

