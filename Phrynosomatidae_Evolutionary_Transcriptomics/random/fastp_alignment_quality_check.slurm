#!/bin/bash
#
#
#SBATCH -J quality_check_fastp
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH --cpus-per-task 20
#SBATCH -t 72:00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/quality_check_fastp.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/quality_check_fastp.err
#SBATCH --account=coxlab
#SBATCH --partition=standard

cd /scratch/dtn2an/workspace/refseq_aligned/quality_test/fastqs/lanes_merged

for i in *_R1_lanes_merged.fastq
do
	name=$(echo $i | awk -F "_R" '{print $1}')
	echo ${name}
	/scratch/dtn2an/tools/fastp -i $i -I ${name}_R2_lanes_merged.fastq -o ${name}_R1.fq -O ${name}_R2.fq --trim_poly_g --trim_poly_x -c -y -3 --cut_tail_window_size 5 --cut_tail_mean_quality 20 -l 36 -q 25
done
