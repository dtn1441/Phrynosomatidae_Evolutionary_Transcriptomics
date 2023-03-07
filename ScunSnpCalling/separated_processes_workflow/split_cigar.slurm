#!/bin/bash
#
#
#SBATCH -J split_cigar_scun_cat
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 0-10:00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/split_cigar_scun_cat.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/split_cigar_scun_cat.err
#SBATCH --account=coxlab
#SBATCH --partition=standard


pre_processing_dir=""


cd ${pre_processing_dir}

module load gatk/4.2.3.0

for i in *_removed_dups.bam
do
	name=$(echo $i | awk -F "_paired" '{print $1}')
	gatk SplitNCigarReads \
		-I $i \
		-O ${name}_splitNcigars.bam \
		-R /scratch/dtn2an/genomes/refseq_genome/ncbi-genomes-2021-09-15/GCF_019175285.1_SceUnd_v1.1_genomic.fna
done




