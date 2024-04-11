#!/bin/bash
#
#
#SBATCH -J PreRecalJointSnpCalling
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=120gb
#SBATCH --output /scratch/dtn2an/workspace/refseq_aligned/stderr/PreRecalJointSnpCalling.%A_%a.out
#SBATCH --error /scratch/dtn2an/workspace/refseq_aligned/stderr/PreRecalJointSnpCalling.%A_%a.err
#SBATCH --account=coxlab
#SBATCH --partition=standard
#SBATCH --array=1-9000%100


pre_processing_dir=""

refgen=""

## Now we can run GenomicsDBImport which will make a database of all of our reference called variants and the reference to aid in joint calling. We will then immediately joint call and then delete the genomics database since it takes up so many files

interval_file=""
myint=$(cat ${interval_file} | sed "${SLURM_ARRAY_TASK_ID}q;d")

cd ${pre_processing_dir}

##Modules may change according to sys admin
module load gcc/11.4.0 openmpi/4.1.4
module load python/3.11.4
module load gatk/4.2.3.0

mkdir joint_call_dir

gatk GenomicsDBImport --batch-size 7 \
  --sample-name-map scun_sample_map.txt \
  --genomicsdb-workspace-path ./joint_call_dir/${myint}_wd \
  --intervals $myint
  
cd joint_call_dir

gatk GenotypeGVCFs \
  -R ${refgen} \
  -V gendb://./${myint}_wd \
  -O ${myint}_jc.vcf

rm -r ./${myint}_wd


