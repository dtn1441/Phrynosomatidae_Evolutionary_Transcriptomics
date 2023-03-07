#!/bin/bash
#
#
#SBATCH -J align
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH --cpus-per-task 20
#SBATCH -t 72:00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/align.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/align.err
#SBATCH --account=coxlab
#SBATCH --partition=standard

module load gcc/9.2.0
module load bwa/0.7.17

cd /scratch/dtn2an/genomes/refseq_genome/ncbi-genomes-2021-09-15/
ref_loc=/scratch/dtn2an/genomes/refseq_genome/ncbi-genomes-2021-09-15/

bwa index -p bwa_index GCF_019175285.1_SceUnd_v1.1_genomic.fna


cd /scratch/dtn2an/workspace/refseq_aligned/quality_test/fastqs/lanes_merged

for i in *_R1.fq
do
	name=$(echo $i | awk -F "_R" '{print $1}')
	echo ${name}
	bwa mem -R ${ref_loc}/bwa_index -P $i ${name}_R2.fq > ${name}.sam 
done
