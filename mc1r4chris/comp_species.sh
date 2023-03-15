#!/bin/bash
#
#
#SBATCH -J stringtie
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 0-5:00:00
#SBATCH --mem 10G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/mc1r_ORF.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/mc1r_ORF.err
#SBATCH --account=coxlab
#SBATCH --partition=standard



module load vcftools
module load samtools
module load gcc/9.2.0
export PATH=$PATH:/scratch/dtn2an/tools/TransDecoder-TransDecoder-v5.7.0
export PATH=$PATH:/scratch/dtn2an/tools/cdhit/


scun_vcfs="/scratch/dtn2an/workspace/refseq_aligned/scun_cat/variant_calls/post_processing/main_chroms"
scja_vcfs="/scratch/dtn2an/workspace/refseq_aligned/scja_snps/"
old_mc1r_workspace="/scratch/dtn2an/workspace/refseq_aligned/scun_cat/variant_calls/post_processing/main_chroms/mc1r4chris"
new_mc1r_workspace="/scratch/dtn2an/workspace/refseq_aligned/mc1r4chris"









##extract TUBB3 from genome based on coordinates found in .gtf

##vcftools --vcf ${scun_vcfs}/NC_056529.1_scun_cat_snps_filtered.vcf --bed ${old_mc1r_workspace}/tubb3.bed --out ${old_mc1r_workspace}/scun_tubb3 --recode --keep-INFO-all
##vcftools --vcf ${scja_vcfs}/chrom8_snps_filtered.vcf --bed ${old_mc1r_workspace}/tubb3.bed --out ${old_mc1r_workspace}/scja_tubb3 --recode --keep-INFO-all

##There were no variable sites in this gene
##Now we can compare TUBB3 from multiple species pretty easily since there is no intraspecific variation (at least in the ones we actually produced vcfs for)
##We can use plate 4646 data that I previously sent through the stringtie workflow and input those into transdecoder




##Finding open reading frames
cd ${new_mc1r_workspace}

for p in *.fa
do
	name=$(echo $p | awk -F _ '{print $1}')
	echo $name
	TransDecoder.LongOrfs -t $p -O $name
done

for q in *.fa
do
	name=$(echo $q | awk -F "_transcripts.fa" '{print $1}')
	echo $name
	TransDecoder.Predict -t $q -O $name
	done

mkdir ${new_mc1r_workspace}/cdhit_output
	
for i in *longest_orfs*
do
	prefix=$(echo $i | awk -F _ '{print $1}')
	suffix=$(echo $i | awk -F . '{print $2}')
	cd-hit -i $i -o ./cdhit_output/${prefix}_cdhit_longest_orfs.${suffix} -c 0.95 -d 0
done

