#!/bin/bash
#
#
#SBATCH -J scun_joint_snp_calling_merged
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=120gb
#SBATCH --output /scratch/dtn2an/workspace/refseq_aligned/stderr/joint_snp_calls_r1.%A_%a.out
#SBATCH --error /scratch/dtn2an/workspace/refseq_aligned/stderr/joint_snp_calls_r1.%A_%a.err
#SBATCH --account=coxlab
#SBATCH --partition=standard
#SBATCH --array=1-9167%100


refgen="/scratch/dtn2an/genomes/refseq_genome/ncbi-genomes-2021-09-15/GCF_019175285.1_SceUnd_v1.1_genomic.fna"

## Now we can run GenomicsDBImport which will make a database of all of our reference called variants and the reference to aid in joint calling. We will then immediately joint call and then delete the ge$

interval_file="/scratch/dtn2an/workspace/refseq_aligned/interval_windowed_list.txt"
myint=$(cat ${interval_file} | sed "${SLURM_ARRAY_TASK_ID}q;d")

cd /scratch/dtn2an/workspace/refseq_aligned/scun_cat/variant_calls/

module load gcc/11.4.0 openmpi/4.1.4
module load python/3.11.4
module load gatk/4.2.3.0

gatk GenomicsDBImport --batch-size 7 \
  --reader-threads 10 \
  --sample-name-map /scratch/dtn2an/workspace/refseq_aligned/vcf_list.txt \
  --genomicsdb-workspace-path ./${myint}_wd \
  --intervals ${myint}


gatk GenotypeGVCFs \
  -R ${refgen} \
  -V gendb://./${myint}_wd \
  --intervals ${myint} \
  --include-non-variant-sites TRUE \
  -O Genotype_gVCFs_new/${myint}_jc_rd1.vcf

