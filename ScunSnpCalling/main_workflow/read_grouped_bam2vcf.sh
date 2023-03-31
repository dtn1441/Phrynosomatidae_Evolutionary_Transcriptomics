## Welcome to 2 years of work that boils down to an unsettlingly small amount of code
## This file will take you from alignment .bam data to variant called .vcf data
## This is not the end all of filtering for quality control, we will need to check the .vcf's afterwards to see if there are bizarre patterns
## From this point, it is a short jump to most popgen statistics

##################################################################################################################################################################################



############################################
## This is the only part that needs input ##
############################################


##provide direct path to reference
ref_genome=/scratch/dtn2an/genomes/refseq_genome/ncbi-genomes-2021-09-15/GCF_019175285.1_SceUnd_v1.1_genomic.fna

##provide direct path to directory with bams
bam_directory=/scratch/dtn2an/workspace/refseq_aligned/scun_merge_runs/H_finished/test4bam2vcf

##provide a brief code for the species, it will be at the start of every file
species=SCUN

##provide path to species specific interval list derived from reference genome
interval_list=/scratch/dtn2an/workspace/refseq_aligned/interval_list.txt












######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

## Start of automatic run ##




module load samtools
module load gatk

cd ${bam_directory}

##RNAseq data is plagued with technical duplicates. This would mess up our allele count later on
##for i in *.bam
##do
##	name=$(echo $i | awk -F . '{print $1}')
##	echo working on ${name}
##	samtools sort $i > ${name}_sorted.bam
##	gatk MarkDuplicates \
##		-I ${name}_sorted.bam \
##		-O ${name}_removed_dups.bam \
##		-M ${species}_marked_dups_metrics \
##		-REMOVE_DUPLICATES true \
##		-READ_NAME_REGEX null
##done

##echo finished duplicate removal

##Splice sites cause problems, so it splits the reads where there are many NNNNNs into spliced pieces 
##for j in *_removed_dups.bam
##do
##	name=$(echo $j | awk -F "_removed" '{print $1}')
##	echo splitting ${name}
##	gatk SplitNCigarReads \
##		-I $j \
##		-O ${name}_splitNcigars.bam \
##		-R ${ref_genome}
##done

##echo finished splitting cigars

##variant calling requires sorted data and an index file that lets it quickly know where to go next
##for k in *_splitNcigars.bam
##do
##	name=$(echo $k | awk -F "_split" '{print $1}')
##	echo sortdexing ${name}
##	samtools sort $k > ${name}_prehaplotypecaller_sort.bam
##	samtools index ${name}_prehaplotypecaller_sort.bam
##done

##echo finished sortdexing


##first we have to find all sites that differ between the individual and its reference genome. It does a de novo assembly of the reads around variable sites because these reads can cause problems in the original alignment
##for l in *_prehaplotypecaller_sort.bam
##do
##	name=$(echo $l | awk -F "_prehaplotype" '{print $1}')
##	echo ${name} goes into haplotypecaller
##	gatk --java-options "-Xmx50g" HaplotypeCaller \
##		-R ${ref_genome} \
##		-I $l \
##		-O ${name}.vcf \
##		--standard-min-confidence-threshold-for-calling 30 \
##		--min-base-quality-score 25 \
##		-ERC GVCF
##done


##echo haplotypecaller finished

##clear up space
rm *_sorted*
rm *_splitNcigars*
rm *_dups_removed*

##This phase creates a database for all the variants across all the samples and calls snps between the individuals (this is the official snp calling part of the workflow)
##We can create a file that tells GATK where to look for the vcfs to include, so we don't have to manually enter the information. The file must be | sample name | hard path | in column format.

ls *.vcf | awk -F . '{print $1}' > sample_list.txt
ls -d "$PWD"/*.vcf > full_path_list.txt

paste sample_list.txt full_path_list.txt > ${species}_map_file.txt

myint=$(cat ${interval_list} | sed "${SLURM_ARRAY_TASK_ID}q;d")




gatk GenomicsDBImport --batch-size 10 \
	--reader-threads 10 \
	--genomicsdb-workspace-path ./${myint}_wd \
	--intervals ${myint} \
	--sample-name-map ${species}_map_file.txt

mkdir ${species}_GenotypeGVCFs \

gatk GenotypeGVCFs \
	-R ${ref_genome} \
	-V gendb://./${myint}_wd \
	-O ${bam_directory}/${species}_GenotypeGVCFs/${myint}_jc.vcf \
	-all-sites ##keep all sites so pixy can calculate nucleotide diversity more robustly

rm -r ./${myint}_wd


## we don't need indels, and need to filter our variants to avoid variants that are unlikely to be real. Finally we exlcude all variants that didn't pass, and we have a .vcf file per region of the genome (chromosome/scaffold).
cd ${bam_directory}/${species}_GenotypeGVCFs/

gatk SelectVariants \
	-V ${myint}_jc.vcf \
	--select-type-to-include SNP \
	--output ${myint}_snps.vcf

gatk VariantFiltration \
	-R ${ref_genome} \
	-V ${myint}_snps.vcf \
	--output ${myint}_snps_marked.vcf \
	--filter-expression "DP < 10 || FS > 30.0 || QUAL< 30 || MQ < 30.0 \\ QD < 2.0" --filter-name "DP10.FS30.QUAL30.MQ30.QD2"

gatk SelectVariants \
	-R ${ref_genome} \
	-V ${myint}_snps_marked.vcf \
	--exclude-filtered true \
	--output ${species}_${myint}_snps_filtered.vcf
