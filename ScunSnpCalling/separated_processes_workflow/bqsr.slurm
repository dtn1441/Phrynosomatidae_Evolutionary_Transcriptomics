#!/bin/bash
#
#
#SBATCH -J bqsr
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 0-50:00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/bqsr.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/bqsr.err
#SBATCH --account=coxlab
#SBATCH --partition=standard

pre_processing_dir=""

ref_fasta=""
known_sites=""

cd ${pre_processing_dir}

##module may change according to sys admin
module load gcc/11.4.0 openmpi/4.1.4 python/3.11.4
module load gatk
module load samtools

mkdir ../bqsr

for i in *_prehaplo_sort.bam
do
        name=$(echo $i | awk -F "_prehaplo" '{print $1}')
        gatk BaseRecalibrator \
              -R ${ref_fasta} \
              -I $i \
              --known-sites ${known_sites} \
              -O ../bqsr/${name}_recal_data.txt
        gatk ApplyBQSR \
                -R ${ref_fasta} \
                -I $i \
                --bqsr-recal-file bqsr/${name}_recal_data.txt \
                -O ../bqsr/${name}_bqsr.bam
done


cd ../bqsr

for i in *_bqsr.bam
do
  	    name=$(echo $i | awk -F "_bqsr" '{print $1}')
        samtools sort $i > ${name}_PostRecalPreHaploSort.bam
        samtools index ${name}_PostRecalPreHaploSort.bam
done
