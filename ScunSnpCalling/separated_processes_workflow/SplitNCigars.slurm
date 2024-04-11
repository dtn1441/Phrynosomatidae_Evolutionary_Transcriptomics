#!/bin/bash
#
#
#SBATCH -J SplitNCigars
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 0-10:00:00
#SBATCH --mem 20G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/SplitNCigars.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/SplitNCigars.err
#SBATCH --account=coxlab
#SBATCH --partition=standard


pre_processing_dir=""
ref_fasta=""

cd ${pre_processing_dir}

##modules below may change according to sys admin
module load gcc/11.4.0 openmpi/4.1.4
module load python/3.11.4
module load gatk/4.2.3.0

for i in *_removed_dups.bam
do
	name=$(echo $i | awk -F "_removed" '{print $1}')
	gatk SplitNCigarReads \
		-I $i \
		-O ${name}_splitNcigars.bam \
		-R ${ref_fasta}
done




