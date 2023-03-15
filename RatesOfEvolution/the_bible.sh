#!/bin/bash
#
#
#SBATCH -J the_bible
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 0-72:00:00
#SBATCH --mem 200G
#SBATCH -o /scratch/dtn2an/sex_conflict_evolution/the_bible.out
#SBATCH -e /scratch/dtn2an/sex_conflict_evolution/the_bible.err
#SBATCH --account=coxlab
#SBATCH --partition=standard

#######################################
########### INPUT STEPS ###############
#######################################

##make the head directory
mkdir /scratch/dtn2an/sex_conflict_evolution/!DIR_NAME!

##create variables so you don't have to touch anything below

##input##
species_number=1,2,3,4,5 ##keep only the integer that describes your data
sp1=species 1 header
sp2=species 2 header
sp3=species 3 header
sp4=species 4 header
sp5=species 5 header

head_dir=/scratch/dtn2an/sex_conflict_evolution/!DIR_NAME!

##lead to the directory that contains each species's fastq files
pre_trim_sp1=path to fastqs
pre_trim_sp2=path to fastqs
pre_trim_sp3=path to fastqs
pre_trim_sp4=path to fastqs
pre_trim_sp5=path to fastqs

## if two species share a reference genome, then make sure they are in the sp1 and sp2 position, if three, then in the sp1 sp2 and sp3 positions.

reference_sp1=/scratch/dtn2an/path_to_reference_genome+annotation ## must already be a hisat indexed genome+annotation
indx_sp1_directory=name of index directory inside of reference directory
indx_name_sp1=actual header of each index file
reference_genome_fasta_sp1=/scratch/dtn2an/path_to_genome/genome.fasta
reference_gff3_sp1=name of .gff3 for reference file

reference_sp2=/scratch/dtn2an/path_to_reference_genome+annotation ## must already be a hisat indexed genome+annotation
indx_sp2_directory=name of index directory inside of reference directory
indx_name_sp2=actual header of each index file
reference_genome_fasta_sp2=/scratch/dtn2an/path_to_genome/genome.fasta
reference_gff3_sp2=name of .gff3 for reference file

reference_sp3=/scratch/dtn2an/path_to_reference_genome+annotation ## must already be a hisat indexed genome+annotation
indx_sp3_directory=name of index directory inside of reference directory
indx_name_sp3=actual header of each index file
reference_genome_fasta_sp3=/scratch/dtn2an/path_to_genome/genome.fasta
reference_gff3_sp3=name of .gff3 for reference file

reference_sp4=/scratch/dtn2an/path_to_reference_genome+annotation ## must already be a hisat indexed genome+annotation
indx_sp4_directory=name of index directory inside of reference directory
indx_name_sp4=actual header of each index file
reference_genome_fasta_sp4=/scratch/dtn2an/path_to_genome/genome.fasta
reference_gff3_sp4=name of .gff3 for reference file

reference_sp5=/scratch/dtn2an/path_to_reference_genome+annotation ## must already be a hisat indexed genome+annotation
indx_sp5_directory=name of index directory inside of reference directory 
indx_name_sp5=actual header of each index file
reference_genome_fasta_sp5=/scratch/dtn2an/path_to_genome/genome.fasta
reference_gff3_sp5=name of .gff3 for reference file


share_genome=1,2,3 ##keep only the integer that describes your data


## species that share a reference genome (up to three) must be in sp1, sp2, and sp3 positions. If only two then use sp1 and sp2.
## a value of one means that none share genome


## rename your samples retaining an ID and the Illumina lane/paired-end information if applicable. eg Scj10_L1_R1.fastq Scj10_L1_R2.fastq Scj10_L2_R1.fastq Scj10_L2_R2.fastq


mkdir ${head_dir}/${sp1}
mkdir ${head_dir}/${sp2}
mkdir ${head_dir}/${sp3}
mkdir ${head_dir}/${sp4}
mkdir ${head_dir}/${sp5}

	

if [[ ${share_genome} == 1] && (${species_number} == 1 || ${species_number} == 2 || ${species_number} == 3 || ${species_number} == 4 || ${species_number} == 5) ]]
then

	###########################################################################
	###################### AUTONOMOUS STEPS ##############################
	###########################################################################

	##Trim_Galore step
	
	##sp1
	mkdir ${head_dir}/${sp1}/trimmed/
	cd ${pre_trim_sp1}/

	module load trimgalore/0.6.4

	for i in *.fastq
		do
			if [[ $i == *"_R1"* || $i == *"_R2"* ]] 	
        	then
        		if [[ $i == *"_R1"* ]]
        		then
        			trim_name=$(echo $i | awk -F _ '{print $1}')
        			echo ${trim_name}
        			trim_galore --fastqc --paired -o ${head_dir}/${sp1}/trimmed ${pre_trim_sp1}/${trim_name}_R1.fastq ${pre_trim_sp1}/${trim_name}_R2.fastq
				else
					:
				fi
			else
				trim_name=$(echo $i | awk -F ".fastq" '{print $1}')
				echo ${trim_name}
				trim_galore --fastqc -o ${head_dir}/${sp1}/trimmed ${pre_trim_sp1}/$i
			fi
		done

	cd ${head_dir}/${sp1}/trimmed
	
	for j in *.fq
        do
        	detrim_name=$(echo $j | awk -F "_trimmed" '{print $1}')
        	echo ${detrim_name}
        	if [[ $j == *"trimmed"* ]]
            then
            	mv ${detrim_name}_trimmed.fq ${detrim_name}.fq
            else
            	:
            fi
    	 done


	for j in *.fq
        do
        	change_name=$(echo $j | awk -F _ '{print $1}')
            echo ${change_name}
            if [[ $j == *"_val_1"* || $j == *"_val_2"* ]]
            then
                if [[ $j == *"_val_1"* ]]
                then
                    mv ${change_name}_R1_val_1.fq ${change_name}_R1.fq
                else
                    mv ${change_name}_R2_val_2.fq ${change_name}_R2.fq
                fi
            else
                :
            fi
        done


	##Hisat2 alignment step
	mkdir ${head_dir}/${sp1}/hisat_aligned/
	cd ${head_dir}/${sp1}/trimmed/

	module load gcc/7.1.0
	module load hisat2/2.1.0


	for j in *.fq
		do
	  		align_name=$(echo $j | awk -F _ '{print $1}')
	        echo ${align_name}
	       
	        if [[ $j == *"_R1"* || $j == *"_R2"* ]]
	        then 
	        	if [[ $j == *"_R1"* ]]	
	        	then
	        		hisat2 -q -x ${reference_sp1}/${indx_sp1_directory}/${indx_name_sp1} -S ${head_dir}/${sp1}/hisat_aligned/${align_name}.sam -1 ${head_dir}/${sp1}/trimmed/${align_name}_R1.fq -2 ${head_dir}/${sp1}/trimmed/${align_name}_R2.fq
	        	else
	        		:
	        	fi
	        else
	        	hisat2 -q -x ${reference_sp1}/${indx_sp1_directory}/${indx_name_sp1} -S ${head_dir}/${sp1}/hisat_aligned/${align_name}.sam -U ${head_dir}/${sp1}/trimmed/${align_name}.fq 
			fi
		done
					 
	##convert to bam and sort

	cd ${head_dir}/${sp1}/hisat_aligned/

	module load samtools/1.10

	for k in *.sam
		do
			convert2bam_name=$(echo $k | awk -F ".sam" '{print $1}')
			echo ${convert2bam_name}
			samtools view -S -b -h $k > ${convert2bam_name}.bam 
		done
	
	for k in *.bam
		do
			sort_name=$(echo $k | awk -F ".bam" '{print $1}')
			echo ${sort_name}	
			samtools sort ${sort_name}.bam -o ${sort_name}_sorted.bam
		done
		
	## Stringtie Transcript Creation

	cd ${head_dir}/${sp1}/hisat_aligned/
	mkdir ${head_dir}/${sp1}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for m in *_sorted.bam
		do
	        stringtie_r1_name=$(echo $m | awk -F "_sorted.bam" '{print $1}')
	        echo ${stringtie_r1_name}
	        stringtie ${head_dir}/${sp1}/hisat_aligned/$m -o ${head_dir}/${sp1}/stringtie/${stringtie_r1_name}.gtf -A -C -G ${reference_sp1}/${reference_gff3_sp1} -l ${sp1}_transcripts
		done

	cd ${head_dir}/${sp1}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	stringtie ${head_dir}/${sp1}/stringtie/* --merge -o ${head_dir}/${sp1}/stringtie/merged.gtf -G ${reference_sp1}/${reference_gff3_sp1} -l merged_transcripts

	mkdir ${head_dir}/${sp1}/stringtie/stringtie_final_output/
	cd ${head_dir}/${sp1}/hisat_aligned/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for n in *_sorted.bam
		do
 	       stringtie_r3_name=$(echo $n | awk -F "_sorted.bam" '{print $1}')
 	       echo ${stringtie_r3_name}
 	       stringtie ${head_dir}/${sp1}/hisat_aligned/$n -o ${head_dir}/${sp1}/stringtie/stringtie_final_output/${stringtie_r3_name}.gtf -A -C -G ${head_dir}/${sp1}/stringtie/merged.gtf
		done

	## Converting transcripts to fasta

	export PATH=$PATH:/scratch/dtn2an/gffread

	cd ${head_dir}/${sp1}/stringtie/stringtie_final_output/
	mkdir ${head_dir}/${sp1}/fastasfromtranscriptome/
	
	for o in *.gtf
		do
 	 		getfasta_name=$(echo $o | awk -F ".gtf" '{print $1}')
 	        echo ${getfasta_name}
 	        gffread -w ${head_dir}/${sp1}/fastasfromtranscriptome/${getfasta_name}_transcripts.fa -g ${reference_sp1}/${reference_genome_fasta_sp1} ${head_dir}/${sp1}/stringtie/stringtie_final_output/$o
		done
	
	##TransDecoder Step

	mkdir ${head_dir}/${sp1}/TransDecoder_${sp1}_output/
	cd ${head_dir}/${sp1}/fastasfromtranscriptome/

	export PATH=$PATH:/scratch/dtn2an/TransDecoder

	for p in *.fa
	do
       	 name=$(echo $p | awk -F "_transcripts.fa" '{print $1}')
       	 echo $name
       	 TransDecoder.LongOrfs -t ${head_dir}/${sp1}/fastasfromtranscriptome/$p -O ${head_dir}/${sp1}/TransDecoder_${sp1}_output/
	done

	cd ${head_dir}/${sp1}/TransDecoder_${sp1}_output/

	for q in *.fa
	do
  		name=$(echo $q | awk -F "_transcripts.fa" '{print $1}')
        echo $name
        TransDecoder.Predict -t ${head_dir}/${sp1}/fastasfromtranscriptome/$q -O ${head_dir}/${sp1}/TransDecoder_${sp1}_output/
	done

	mkdir ${head_dir}/${sp1}/CD_HIT_output
	
	export PATH=$PATH:/scratch/dtn2an/cd-hit-v4.8.1-2019-0228

	./cd-hit -i ${head_dir}/${sp1}/TransDecoder_${sp1}_output/longest_orfs.pep -o ${head_dir}/${sp1}/CD_HIT_output/cdhit_longest_orfs.pep -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp1}/TransDecoder_${sp1}_output/longest_orfs.gff3 -o ${head_dir}/${sp1}/CD_HIT_output/cdhit_longest_orfs.gff3 -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp1}/TransDecoder_${sp1}_output/longest_orfs.cds -o ${head_dir}/${sp1}/CD_HIT_output/cdhit_longest_orfs.cds -c 0.95 -d 0

	echo "species one, not shared done"
else
	echo "species one, not shared skipped"
fi


##########################################################################################################################################

if [[ ${share_genome} == 1 && (${species_number} == 2 || ${species_number} == 3 || ${species_number} == 4 || ${species_number} == 5) ]]
then

	###########################################################################
	###################### AUTONOMOUS STEPS ##############################
	###########################################################################

	##Trim_Galore step
	
	##sp2
	mkdir ${head_dir}/${sp2}/trimmed/
	cd /${pre_trim_sp2}/

	module load trimgalore/0.6.4

	for i in *.fastq
		do
			if [[ $i == *"_R1"* || $i == *"_R2"* ]] 	
        	then
        		if [[ $i == *"_R1"* ]]
        		then
        			trim_name=$(echo $i | awk -F _ '{print $1}')
        			echo ${trim_name}
        			trim_galore --fastqc --paired -o ${head_dir}/${sp2}/trimmed ${pre_trim_sp2}/${trim_name}_R1.fastq ${pre_trim_sp2}/${trim_name}_R2.fastq
				else
					:
				fi
			else
				trim_name=$(echo $i | awk -F ".fastq" '{print $1}')
				echo ${trim_name}
				trim_galore --fastqc -o ${head_dir}/${sp2}/trimmed ${pre_trim_sp2}/$i
			fi
		done

	cd ${head_dir}/${sp2}/trimmed

	for j in *.fq
        do
        	detrim_name=$(echo $j | awk -F "_trimmed" '{print $1}')
        	echo ${detrim_name}
        	if [[ $j == *"trimmed"* ]]
            then
            	mv ${detrim_name}_trimmed.fq ${detrim_name}.fq
            else
            	:
            fi
    	done


	for j in *.fq
        do
        	change_name=$(echo $j | awk -F _ '{print $1}')
            echo ${change_name}
            if [[ $j == *"_val_1"* || $j == *"_val_2"* ]]
            then
                if [[ $j == *"_val_1"* ]]
                then
                    mv ${change_name}_R1_val_1.fq ${change_name}_R1.fq
                else
                    mv ${change_name}_R2_val_2.fq ${change_name}_R2.fq
                fi
            else
                :
            fi
        done
			
	##Hisat2 alignment step
	mkdir ${head_dir}/${sp2}/hisat_aligned/
	cd ${head_dir}/${sp2}/trimmed/

	module load gcc/7.1.0
	module load hisat2/2.1.0

	for j in *.fq
		do
	  		align_name=$(echo $j | awk -F _ '{print $1}')
	        echo ${align_name}
	       
	        if [[ $j == *"_R1"* || $j == *"_R2"* ]]
	        then 
	        	if [[ $j == *"_R1"* ]]
	        	then	
	        		hisat2 -q -x ${reference_sp2}/${indx_sp2_directory}/${indx_name_sp2} -S ${head_dir}/${sp2}/hisat_aligned/${align_name}.sam -1 ${head_dir}/${sp2}/trimmed/${align_name}_R1.fq -2 ${head_dir}/${sp2}/trimmed/${align_name}_R2.fq
				else
					:
				fi	        
	        else
	        	hisat2 -q -x ${reference_sp2}/${indx_sp2_directory}/${indx_name_sp2} -S ${head_dir}/${sp2}/hisat_aligned/${align_name}.sam -U ${head_dir}/${sp2}/trimmed/${align_name}.fq 
			fi
		done

	##convert to bam and sort

	cd ${head_dir}/${sp2}/hisat_aligned/

	module load samtools/1.10

	for k in *.sam
		do
			convert2bam_name=$(echo $k | awk -F ".sam" '{print $1}')
			echo ${convert2bam_name}
			samtools view -S -b -h $k > ${convert2bam_name}.bam 
		done
	
	for k in *.bam
		do
			sort_name=$(echo $k | awk -F ".bam" '{print $1}')
			echo ${sort_name}	
			samtools sort ${sort_name}.bam -o ${sort_name}_sorted.bam
		done

	## Stringtie Transcript Creation

	cd ${head_dir}/${sp2}/hisat_aligned/
	mkdir ${head_dir}/${sp2}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for m in *_sorted.bam
		do
	        stringtie_r1_name=$(echo $m | awk -F "_sorted.bam" '{print $1}')
	        echo ${stringtie_r1_name}
	        stringtie ${head_dir}/${sp2}/hisat_aligned/$m -o ${head_dir}/${sp2}/stringtie/${stringtie_r1_name}.gtf -A -C -G ${reference_sp2}/${reference_gff3_sp2} -l ${sp2}_transcripts
		done

	cd ${head_dir}/${sp2}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	stringtie ${head_dir}/${sp2}/stringtie/* --merge -o ${head_dir}/${sp2}/stringtie/merged.gtf -G ${reference_sp2}/${reference_gff3_sp2} -l merged_transcripts

	mkdir ${head_dir}/${sp2}/stringtie/stringtie_final_output/
	cd ${head_dir}/${sp2}/hisat_aligned/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for n in *_sorted.bam
		do
 	       stringtie_r3_name=$(echo $n | awk -F "_sorted.bam" '{print $1}')
 	       echo ${stringtie_r3_name}
 	       stringtie ${head_dir}/${sp2}/hisat_aligned/$n -o ${head_dir}/${sp2}/stringtie/stringtie_final_output/${stringtie_r3_name}.gtf -A -C -G ${head_dir}/${sp2}/stringtie/merged.gtf
		done

	## Converting transcripts to fasta

	export PATH=$PATH:/scratch/dtn2an/gffread

	cd ${head_dir}/${sp2}/stringtie/stringtie_final_output/
	mkdir ${head_dir}/${sp2}/fastasfromtranscriptome/
	
	for o in *.gtf
		do
 	 		getfasta_name=$(echo $o | awk -F ".gtf" '{print $1}')
 	        echo ${getfasta_name}
 	        gffread -w ${head_dir}/${sp2}/fastasfromtranscriptome/${getfasta_name}_transcripts.fa -g ${reference_sp2}/${reference_genome_fasta_sp2} ${head_dir}/${sp2}/stringtie/stringtie_final_output/$o
		done

	##TransDecoder Step

	mkdir ${head_dir}/${sp2}/TransDecoder_${sp2}_output/
	cd ${head_dir}/${sp2}/fastasfromtranscriptome/

	export PATH=$PATH:/scratch/dtn2an/TransDecoder

	for p in *.fa
		do
       	 	name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		 echo $name
       	 	TransDecoder.LongOrfs -t ${head_dir}/${sp2}/fastasfromtranscriptome/$p -O ${head_dir}/${sp2}/TransDecoder_${sp2}_output/
		done

	cd ${head_dir}/${sp2}/TransDecoder_${sp2}_output/

	for q in *.fa
		do
  		 name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       	 echo $name
       	 TransDecoder.Predict -t ${head_dir}/${sp2}/fastasfromtranscriptome/$q -O ${head_dir}/${sp2}/TransDecoder_${sp2}_output/
		done
	
	mkdir ${head_dir}/${sp1}/CD_HIT_output
	
	export PATH=$PATH:/scratch/dtn2an/cd-hit-v4.8.1-2019-0228

	./cd-hit -i ${head_dir}/${sp2}/TransDecoder_${sp2}_output/longest_orfs.pep -o ${head_dir}/${sp2}/CD_HIT_output/cdhit_longest_orfs.pep -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp2}/TransDecoder_${sp2}_output/longest_orfs.gff3 -o ${head_dir}/${sp2}/CD_HIT_output/cdhit_longest_orfs.gff3 -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp2}/TransDecoder_${sp2}_output/longest_orfs.cds -o ${head_dir}/${sp2}/CD_HIT_output/cdhit_longest_orfs.cds -c 0.95 -d 0

	echo "species two, not shared done"
else
	echo "species two, not shared skipped"
fi


##########################################################################################################################################

if [[ ${share_genome} == 1 && (${species_number} == 3 || ${species_number}== 4 || ${species_number} == 5) ]]
then

	###########################################################################
	###################### AUTONOMOUS STEPS ##############################
	###########################################################################

	##Trim_Galore step
	
	##sp3
	mkdir ${head_dir}/${sp3}/trimmed/
	cd /${pre_trim_sp3}/

	module load trimgalore/0.6.4

	for i in *.fastq
		do
			if [[ $i == *"_R1"* || $i == *"_R2"* ]] 	
        	then
        		if [[ $i == *"_R1"* ]]
        		then
        			trim_name=$(echo $i | awk -F _ '{print $1}')
        			echo ${trim_name}
        			trim_galore --fastqc --paired -o ${head_dir}/${sp3}/trimmed ${pre_trim_sp3}/${trim_name}_R1.fastq ${pre_trim_sp3}/${trim_name}_R2.fastq
				else
					:
				fi
			else
				trim_name=$(echo $i | awk -F ".fastq" '{print $1}')
				echo ${trim_name}
				trim_galore --fastqc -o ${head_dir}/${sp3}/trimmed ${pre_trim_sp3}/$i
			fi
		done
	
	cd ${head_dir}/${sp3}/trimmed
		
	for j in *.fq
        do
        	detrim_name=$(echo $j | awk -F "_trimmed" '{print $1}')
        	echo ${detrim_name}
        	if [[ $j == *"trimmed"* ]]
            then
            	mv ${detrim_name}_trimmed.fq ${detrim_name}.fq
            else
            	:
            fi
    	 done


	for j in *.fq
        do
        	change_name=$(echo $j | awk -F _ '{print $1}')
            echo ${change_name}
            if [[ $j == *"_val_1"* || $j == *"_val_2"* ]]
            then
                if [[ $j == *"_val_1"* ]]
                then
                    mv ${change_name}_R1_val_1.fq ${change_name}_R1.fq
                else
                    mv ${change_name}_R2_val_2.fq ${change_name}_R2.fq
                fi
            else
                :
            fi
        done
	
	##Hisat2 alignment step
	mkdir ${head_dir}/${sp3}/hisat_aligned/
	cd ${head_dir}/${sp3}/trimmed/

	module load gcc/7.1.0
	module load hisat2/2.1.0

	for j in *.fq
		do
	  		align_name=$(echo $j | awk -F _ '{print $1}')
	        echo ${align_name}
	       
	        if [[ $j == *"_R1"* || $j == *"_R2"* ]]
	        then 	
	        	if [[ $j == *"_R1"* ]]
	        	then
	        		hisat2 -q -x ${reference_sp3}/${indx_sp3_directory}/${indx_name_sp3} -S ${head_dir}/${sp3}/hisat_aligned/${align_name}.sam -1 ${head_dir}/${sp3}/trimmed/${align_name}_R1.fq -2 ${head_dir}/${sp3}/trimmed/${align_name}_R2.fq
	        	else
	        		:
	        	fi
	        else
	        	hisat2 -q -x ${reference_sp3}/${indx_sp3_directory}/${indx_name_sp3} -S ${head_dir}/${sp3}/hisat_aligned/${align_name}.sam -U ${head_dir}/${sp3}/trimmed/${align_name}.fq 
			fi
		done
	
	##convert to bam and sort

	cd ${head_dir}/${sp3}/hisat_aligned/

	module load samtools/1.10

	for k in *.sam
		do
			convert2bam_name=$(echo $k | awk -F ".sam" '{print $1}')
			echo ${convert2bam_name}
			samtools view -S -b -h $k > ${convert2bam_name}.bam 
		done
	
	for k in *.bam
		do
			sort_name=$(echo $k | awk -F ".bam" '{print $1}')
			echo ${sort_name}	
			samtools sort ${sort_name}.bam -o ${sort_name}_sorted.bam
		done

	## Stringtie Transcript Creation

	cd ${head_dir}/${sp3}/hisat_aligned/
	mkdir ${head_dir}/${sp3}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for m in *_sorted.bam
		do
	        stringtie_r1_name=$(echo $m | awk -F "_sorted.bam" '{print $1}')
	        echo ${stringtie_r1_name}
	        stringtie ${head_dir}/${sp3}/hisat_aligned/$m -o ${head_dir}/${sp3}/stringtie/${stringtie_r1_name}.gtf -A -C -G ${reference_sp3}/${reference_gff3_sp3} -l ${sp3}_transcripts
		done

	cd ${head_dir}/${sp3}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	stringtie ${head_dir}/${sp3}/stringtie/* --merge -o ${head_dir}/${sp3}/stringtie/merged.gtf -G ${reference_sp3}/${reference_gff3_sp3} -l merged_transcripts

	mkdir ${head_dir}/${sp3}/stringtie/stringtie_final_output/
	cd ${head_dir}/${sp3}/hisat_aligned/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for n in *_sorted.bam
		do
 	   	   stringtie_r3_name=$(echo $n | awk -F "_sorted.bam" '{print $1}')
 	       echo ${stringtie_r3_name}
 	       stringtie ${head_dir}/${sp3}/hisat_aligned/$n -o ${head_dir}/${sp3}/stringtie/stringtie_final_output/${stringtie_r3_name}.gtf -A -C -G ${head_dir}/${sp3}/stringtie/merged.gtf
		done

	## Converting transcripts to fasta

	export PATH=$PATH:/scratch/dtn2an/gffread

	cd ${head_dir}/${sp3}/stringtie/stringtie_final_output/
	mkdir ${head_dir}/${sp3}/fastasfromtranscriptome/
	
	for o in *.gtf
		do
 	 		getfasta_name=$(echo $o | awk -F ".gtf" '{print $1}')
 	        echo ${getfasta_name}
 	        gffread -w ${head_dir}/${sp3}/fastasfromtranscriptome/${getfasta_name}_transcripts.fa -g ${reference_sp3}/${reference_genome_fasta_sp3} ${head_dir}/${sp3}/stringtie/stringtie_final_output/$o
		done
		
	##TransDecoder Step

	mkdir ${head_dir}/${sp3}/TransDecoder_${sp3}_output/
	cd ${head_dir}/${sp3}/fastasfromtranscriptome/

	export PATH=$PATH:/scratch/dtn2an/TransDecoder

	for p in *.fa
		do
       	 	name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		 echo $name
       	 	TransDecoder.LongOrfs -t ${head_dir}/${sp3}/fastasfromtranscriptome/$p -O ${head_dir}/${sp3}/TransDecoder_${sp3}_output/
		done

	cd ${head_dir}/${sp3}/TransDecoder_${sp3}_output/

	for q in *.fa
		do
  		 	name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       	 	echo $name
       		TransDecoder.Predict -t ${head_dir}/${sp3}/fastasfromtranscriptome/$q -O ${head_dir}/${sp3}/TransDecoder_${sp3}_output/
		done
	
	mkdir ${head_dir}/${sp3}/CD_HIT_output
	
	export PATH=$PATH:/scratch/dtn2an/cd-hit-v4.8.1-2019-0228

	./cd-hit -i ${head_dir}/${sp3}/TransDecoder_${sp3}_output/longest_orfs.pep -o ${head_dir}/${sp3}/CD_HIT_output/cdhit_longest_orfs.pep -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp3}/TransDecoder_${sp3}_output/longest_orfs.gff3 -o ${head_dir}/${sp3}/CD_HIT_output/cdhit_longest_orfs.gff3 -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp3}/TransDecoder_${sp3}_output/longest_orfs.cds -o ${head_dir}/${sp3}/CD_HIT_output/cdhit_longest_orfs.cds -c 0.95 -d 0

	echo "species three, not shared done"
	
else
	echo "species three, not shared skipped"		
fi


##########################################################################################################################################

if [[ ${share_genome} == 1 && (${species_number} == 4 || ${species_number} == 5) ]]
then

	###########################################################################
	###################### AUTONOMOUS STEPS ##############################
	###########################################################################

	##Trim_Galore step
	
	##sp4
	mkdir ${head_dir}/${sp4}/trimmed/
	cd /${pre_trim_sp4}/

	module load trimgalore/0.6.4

	for i in *.fastq
		do
			if [[ $i == *"_R1"* || $i == *"_R2"* ]] 	
        	then
        		if [[ $i == *"_R1"* ]]
        		then
        			trim_name=$(echo $i | awk -F _ '{print $1}')
        			echo ${trim_name}
        			trim_galore --fastqc --paired -o ${head_dir}/${sp4}/trimmed ${pre_trim_sp4}/${trim_name}_R1.fastq ${pre_trim_sp4}/${trim_name}_R2.fastq
				else
					:
				fi
			else
				trim_name=$(echo $i | awk -F ".fastq" '{print $1}')
				echo ${trim_name}
				trim_galore --fastqc -o ${head_dir}/${sp4}/trimmed ${pre_trim_sp4}/$i
			fi
		done
	
	cd ${head_dir}/${sp4}/trimmed
		
	for j in *.fq
        do
        	detrim_name=$(echo $j | awk -F "_trimmed" '{print $1}')
        	echo ${detrim_name}
        	if [[ $j == *"trimmed"* ]]
            then
            	mv ${detrim_name}_trimmed.fq ${detrim_name}.fq
            else
            	:
            fi
    	 done

	for j in *.fq
        do
        	change_name=$(echo $j | awk -F _ '{print $1}')
            echo ${change_name}
            if [[ $j == *"_val_1"* || $j == *"_val_2"* ]]
            then
                if [[ $j == *"_val_1"* ]]
                then
                    mv ${change_name}_R1_val_1.fq ${change_name}_R1.fq
                else
                    mv ${change_name}_R2_val_2.fq ${change_name}_R2.fq
                fi
            else
                :
            fi
        done
        
	##Hisat2 alignment step
	mkdir ${head_dir}/${sp4}/hisat_aligned/
	cd ${head_dir}/${sp4}/trimmed/

	module load gcc/7.1.0
	module load hisat2/2.1.0

	for j in *.fq
		do
	  		align_name=$(echo $j | awk -F _ '{print $1}')
	        echo ${align_name}
	       
	        if [[ $j == *"_R1"* || $j == *"_R2"* ]]
	        then 	
	        	if [[ $j == *"_R1"* ]]
	        	then
	        		hisat2 -q -x ${reference_sp4}/${indx_sp4_directory}/${indx_name_sp4} -S ${head_dir}/${sp4}/hisat_aligned/${align_name}.sam -1 ${head_dir}/${sp4}/trimmed/${align_name}_R1.fq -2 ${head_dir}/${sp4}/trimmed/${align_name}_R2.fq
	        	else
	        		:
	        	fi
	        else
	        	hisat2 -q -x ${reference_sp4}/${indx_sp4_directory}/${indx_name_sp4} -S ${head_dir}/${sp4}/hisat_aligned/${align_name}.sam -U ${head_dir}/${sp4}/trimmed/${align_name}.fq 
			fi
		done

	##convert to bam and sort

	cd ${head_dir}/${sp4}/hisat_aligned/


	module load samtools/1.10

	for k in *.sam
		do
			convert2bam_name=$(echo $k | awk -F ".sam" '{print $1}')
			echo ${convert2bam_name}
			samtools view -S -b -h $k > ${convert2bam_name}.bam 
		done
	
	for k in *.bam
		do
			sort_name=$(echo $k | awk -F ".bam" '{print $1}')
			echo ${sort_name}	
			samtools sort ${sort_name}.bam -o ${sort_name}_sorted.bam
		done

	## Stringtie Transcript Creation

	cd ${head_dir}/${sp4}/hisat_aligned/
	mkdir ${head_dir}/${sp4}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for m in *_sorted.bam
		do
	        stringtie_r1_name=$(echo $m | awk -F "_sorted.bam" '{print $1}')
	        echo ${stringtie_r1_name}
	        stringtie ${head_dir}/${sp4}/hisat_aligned/$m -o ${head_dir}/${sp4}/stringtie/${stringtie_r1_name}.gtf -A -C -G ${reference_sp4}/${reference_gff3_sp4} -l ${sp4}_transcripts
		done

	cd ${head_dir}/${sp4}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	stringtie ${head_dir}/${sp4}/stringtie/* --merge -o ${head_dir}/${sp4}/stringtie/merged.gtf -G ${reference_sp4}/${reference_gff3_sp4} -l merged_transcripts

	mkdir ${head_dir}/${sp4}/stringtie/stringtie_final_output/
	cd ${head_dir}/${sp4}/hisat_aligned/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for n in *_sorted.bam
		do
 	       stringtie_r3_name=$(echo $n | awk -F "_sorted.bam" '{print $1}')
 	       echo ${stringtie_r3_name}
 	       stringtie ${head_dir}/${sp4}/hisat_aligned/$n -o ${head_dir}/${sp4}/stringtie/stringtie_final_output/${stringtie_r3_name}.gtf -A -C -G ${head_dir}/${sp4}/stringtie/merged.gtf
		done

	## Converting transcripts to fasta

	export PATH=$PATH:/scratch/dtn2an/gffread

	cd ${head_dir}/${sp4}/stringtie/stringtie_final_output/
	mkdir ${head_dir}/${sp4}/fastasfromtranscriptome/
	
	for o in *.gtf
		do
 	 		getfasta_name=$(echo $o | awk -F ".gtf" '{print $1}')
 	        echo ${getfasta_name}
 	        gffread -w ${head_dir}/${sp4}/${getfasta_name}_transcripts.fa -g ${reference_sp4}/${reference_genome_fasta_sp4} ${head_dir}/${sp4}/stringtie/stringtie_final_output/$o
		done
	##TransDecoder Step

	mkdir ${head_dir}/${sp4}/TransDecoder_${sp4}_output/
	cd ${head_dir}/${sp4}/fastasfromtranscriptome/

	export PATH=$PATH:/scratch/dtn2an/TransDecoder

	for p in *.fa
		do
       	 	name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		echo $name
       	 	TransDecoder.LongOrfs -t ${head_dir}/${sp4}/fastasfromtranscriptome/$p -O ${head_dir}/${sp4}/TransDecoder_${sp4}_output/
		done

	cd ${head_dir}/${sp4}/TransDecoder_${sp4}_output/

	for q in *.fa
		do
  			 name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		 echo $name
       		 TransDecoder.Predict -t ${head_dir}/${sp4}/fastasfromtranscriptome/$q -O ${head_dir}/${sp4}/TransDecoder_${sp4}_output/
		done
	
	mkdir ${head_dir}/${sp4}/CD_HIT_output
	
	export PATH=$PATH:/scratch/dtn2an/cd-hit-v4.8.1-2019-0228

	./cd-hit -i ${head_dir}/${sp4}/TransDecoder_${sp4}_output/longest_orfs.pep -o ${head_dir}/${sp4}/CD_HIT_output/cdhit_longest_orfs.pep -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp4}/TransDecoder_${sp4}_output/longest_orfs.gff3 -o ${head_dir}/${sp4}/CD_HIT_output/cdhit_longest_orfs.gff3 -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp4}/TransDecoder_${sp4}_output/longest_orfs.cds -o ${head_dir}/${sp4}/CD_HIT_output/cdhit_longest_orfs.cds -c 0.95 -d 0
		
	echo "species four, not shared done"
else
	echo "species four, not shared skipped"
fi

##########################################################################################################################################

if [[ ${share_genome} == 1 && (${species_number} == 5) ]]
then

	###########################################################################
	###################### AUTONOMOUS STEPS ##############################
	###########################################################################

	##Trim_Galore step
	
	##sp5
	mkdir ${head_dir}/${sp5}/trimmed/
	cd /${pre_trim_sp5}/

	module load trimgalore/0.6.4

	for i in *.fastq
		do
			if [[ $i == *"_R1"* || $i == *"_R2"* ]] 	
        	then
        		if [[ $i == *"_R1"* ]]
        		then
        			trim_name=$(echo $i | awk -F _ '{print $1}')
        			echo ${trim_name}
        			trim_galore --fastqc --paired -o ${head_dir}/${sp5}/trimmed ${pre_trim_sp5}/${trim_name}_R1.fastq ${pre_trim_sp5}/${trim_name}_R2.fastq
				else
					:
				fi
			else
				trim_name=$(echo $i | awk -F ".fastq" '{print $1}')
				echo ${trim_name}
				trim_galore --fastqc -o ${head_dir}/${sp5}/trimmed ${pre_trim_sp5}/$i
			fi
		done
	
	cd ${head_dir}/${sp5}/trimmed
	
	for j in *.fq
        do
        	detrim_name=$(echo $j | awk -F "_trimmed" '{print $1}')
        	echo ${detrim_name}
        	if [[ $j == *"trimmed"* ]]
            then
            	mv ${detrim_name}_trimmed.fq ${detrim_name}.fq
            else
            	:
            fi
    	 done


	for j in *.fq
        do
        	change_name=$(echo $j | awk -F _ '{print $1}')
            echo ${change_name}
            if [[ $j == *"_val_1"* || $j == *"_val_2"* ]]
            then
                if [[ $j == *"_val_1"* ]]
                then
                    mv ${change_name}_R1_val_1.fq ${change_name}_R1.fq
                else
                    mv ${change_name}_R2_val_2.fq ${change_name}_R2.fq
                fi
            else
                :
            fi
        done

	##Hisat2 alignment step
	mkdir ${head_dir}/${sp5}/hisat_aligned/
	cd ${head_dir}/${sp5}/trimmed/

	module load gcc/7.1.0
	module load hisat2/2.1.0

	for j in *.fq
		do
	  		align_name=$(echo $j | awk -F _ '{print $1}')
	        echo ${align_name}
	       
	        if [[ $j == *"_R1"* || $j == *"_R2"* ]]
	        then 
	        	if [[ $j == *"_R1" ]]
	        	then	
	        		hisat2 -q -x ${reference_sp5}/${indx_sp5_directory}/${indx_name_sp5} -S ${head_dir}/${sp5}/hisat_aligned/${align_name}.sam -1 ${head_dir}/${sp5}/trimmed/${align_name}_R1.fq -2 ${head_dir}/${sp5}/trimmed/${align_name}_R2.fq
	        	else
	        		:
	        	fi
	        else
	        	hisat2 -q -x ${reference_sp5}/${indx_sp5_directory}/${indx_name_sp5} -S ${head_dir}/${sp5}/hisat_aligned/${align_name}.sam -U ${head_dir}/${sp5}/trimmed/${align_name}.fq 
			fi
		done

	##convert to bam and sort

	cd ${head_dir}/${sp5}/hisat_aligned/

	module load samtools/1.10

	for k in *.sam
		do
			convert2bam_name=$(echo $k | awk -F ".sam" '{print $1}')
			echo ${convert2bam_name}
			samtools view -S -b -h $k > ${convert2bam_name}.bam 
		done
	
	for k in *.bam
		do
			sort_name=$(echo $k | awk -F ".bam" '{print $1}')
			echo ${sort_name}	
			samtools sort ${sort_name}.bam -o ${sort_name}_sorted.bam
		done

	## Stringtie Transcript Creation

	cd ${head_dir}/${sp5}/hisat_aligned/
	mkdir ${head_dir}/${sp5}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for m in *_sorted.bam
		do
	        stringtie_r1_name=$(echo $m | awk -F "_sorted.bam" '{print $1}')
	        echo ${stringtie_r1_name}
	        stringtie ${head_dir}/${sp5}/hisat_aligned/$m -o ${head_dir}/${sp5}/stringtie/${stringtie_r1_name}.gtf -A -C -G ${reference_sp5}/${reference_gff3_sp5} -l ${sp5}_transcripts
		done

	cd ${head_dir}/${sp5}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	stringtie ${head_dir}/${sp5}/stringtie/* --merge -o ${head_dir}/${sp5}/stringtie/merged.gtf -G ${reference_sp5}/${reference_gff3_sp5} -l merged_transcripts

	mkdir ${head_dir}/${sp5}/stringtie/stringtie_final_output/
	cd ${head_dir}/${sp5}/hisat_aligned/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for n in *_sorted.bam
		do
 	       stringtie_r3_name=$(echo $n | awk -F "_sorted.bam" '{print $1}')
 	       echo ${stringtie_r3_name}
 	       stringtie ${head_dir}/${sp5}/hisat_aligned/$n -o ${head_dir}/${sp5}/stringtie/stringtie_final_output/${stringtie_r3_name}.gtf -A -C -G ${head_dir}/${sp5}/stringtie/merged.gtf
		done

	## Converting transcripts to fasta

	export PATH=$PATH:/scratch/dtn2an/gffread

	cd ${head_dir}/${sp5}/stringtie/stringtie_final_output/
	mkdir ${head_dir}/${sp5}/fastasfromtranscriptome/
	
	for o in *.gtf
		do
 	 		getfasta_name=$(echo $o | awk -F ".gtf" '{print $1}')
 	        echo ${getfasta_name}
 	        gffread -w ${head_dir}/${sp5}/fastasfromtranscriptome/${getfasta_name}_transcripts.fa -g ${reference_sp5}/${reference_genome_fasta_sp5} ${head_dir}/${sp5}/stringtie/stringtie_final_output/$o
		done
		
	##TransDecoder Step

	mkdir ${head_dir}/${sp5}/TransDecoder_${sp5}_output/
	cd ${head_dir}/${sp5}/fastasfromtranscriptome/

	export PATH=$PATH:/scratch/dtn2an/TransDecoder

	for p in *.fa
		do
       	 	 name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		 echo $name
       	 	 TransDecoder.LongOrfs -t ${head_dir}/${sp5}/fastasfromtranscriptome/$p -O ${head_dir}/${sp5}/TransDecoder_${sp5}_output/
		done

	cd ${head_dir}/${sp5}/TransDecoder_${sp5}_output/

	for q in *.fa
		do
  		 	 name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		 echo $name
       		 TransDecoder.Predict -t ${head_dir}/${sp5}/fastasfromtranscriptome/$q -O ${head_dir}/${sp5}/TransDecoder_${sp5}_output/
		done
		
	mkdir ${head_dir}/${sp5}/CD_HIT_output
	
	export PATH=$PATH:/scratch/dtn2an/cd-hit-v4.8.1-2019-0228

	./cd-hit -i ${head_dir}/${sp5}/TransDecoder_${sp5}_output/longest_orfs.pep -o ${head_dir}/${sp5}/CD_HIT_output/cdhit_longest_orfs.pep -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp5}/TransDecoder_${sp5}_output/longest_orfs.gff3 -o ${head_dir}/${sp5}/CD_HIT_output/cdhit_longest_orfs.gff3 -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp5}/TransDecoder_${sp5}_output/longest_orfs.cds -o ${head_dir}/${sp5}/CD_HIT_output/cdhit_longest_orfs.cds -c 0.95 -d 0
	
	echo "species five, not shared done"
else
	echo "species five, not shared skipped"
fi

##########################################################################################################################################
##########################################################################################################################################

if [[ ${share_genome} == 2 && (${species_number} == 2 || ${species_number} == 3 || ${species_number} == 4 || ${species_number} == 5) ]]
then


	###########################################################################
	###################### AUTONOMOUS STEPS ##############################
	###########################################################################

	##Trim_Galore step
	
	##sp1
	mkdir ${head_dir}/${sp1}/trimmed/
	cd ${pre_trim_sp1}/

	module load trimgalore/0.6.4

	for i in *.fastq
		do
			if [[ $i == *"_R1"* || $i == *"_R2"* ]] 	
        	then
        		if [[ $i == *"_R1"* ]]
        		then
        			trim_name=$(echo $i | awk -F _ '{print $1}')
        			echo ${trim_name}
        			trim_galore --fastqc --paired -o ${head_dir}/${sp1}/trimmed ${pre_trim_sp1}/${trim_name}_R1.fastq ${pre_trim_sp1}/${trim_name}_R2.fastq
				else
					:
				fi
			else
				trim_name=$(echo $i | awk -F ".fastq" '{print $1}')
				echo ${trim_name}
				trim_galore --fastqc -o ${head_dir}/${sp1}/trimmed ${pre_trim_sp1}/$i
			fi
		done
		
	cd ${head_dir}/${sp1}/trimmed

	for j in *.fq
        do
        	detrim_name=$(echo $j | awk -F "_trimmed" '{print $1}')
        	echo ${detrim_name}
        	if [[ $j == *"trimmed"* ]]
            then
            	mv ${detrim_name}_trimmed.fq ${detrim_name}.fq
            else
            	:
            fi
    	 done


	for j in *.fq
        do
        	change_name=$(echo $j | awk -F _ '{print $1}')
            echo ${change_name}
            if [[ $j == *"_val_1"* || $j == *"_val_2"* ]]
            then
                if [[ $j == *"_val_1"* ]]
                then
                    mv ${change_name}_R1_val_1.fq ${change_name}_R1.fq
                else
                    mv ${change_name}_R2_val_2.fq ${change_name}_R2.fq
                fi
            else
                :
            fi
        done
		
	##Hisat2 alignment step
	mkdir ${head_dir}/${sp1}/hisat_aligned/
	cd ${head_dir}/${sp1}/trimmed/

	module load gcc/7.1.0
	module load hisat2/2.1.0

	for j in *.fq
		do
	  		align_name=$(echo $j | awk -F _ '{print $1}')
	        echo ${align_name}
	       
	        if [[ $j == *"_R1"* || $j = *"_R2"* ]]
	        then 	
	        	if [[ $j == *"_R1"* ]]
	        	then
	        		hisat2 -q -x ${reference_sp1}/${indx_sp1_directory}/${indx_name_sp1} -S ${head_dir}/${sp1}/hisat_aligned/${align_name}.sam -1 ${head_dir}/${sp1}/trimmed/${align_name}_R1.fq -2 ${head_dir}/${sp1}/trimmed/${align_name}_R2.fq
	        	else
	        		:
	        	fi
	        else
	        	hisat2 -q -x ${reference_sp1}/${indx_sp1_directory}/${indx_name_sp1} -S ${head_dir}/${sp1}/hisat_aligned/${align_name}.sam -U ${head_dir}/${sp1}/trimmed/${align_name}.fq 
			fi
		done

	##convert to bam and sort

	cd ${head_dir}/${sp1}/hisat_aligned/

	module load samtools/1.10

	for k in *.sam
		do
			convert2bam_name=$(echo $k | awk -F ".sam" '{print $1}')
			echo ${convert2bam_name}
			samtools view -S -b -h $k > ${convert2bam_name}.bam 
		done
	
	for k in *.bam
		do
			sort_name=$(echo $k | awk -F ".bam" '{print $1}')
			echo ${sort_name}	
			samtools sort ${sort_name}.bam -o ${sort_name}_sorted.bam
		done
	
	## Stringtie Transcript Creation

	cd ${head_dir}/${sp1}/hisat_aligned/
	mkdir ${head_dir}/${sp1}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for m in *_sorted.bam
		do
	        stringtie_r1_name=$(echo $m | awk -F "_sorted.bam" '{print $1}')
	        echo ${stringtie_r1_name}
	        stringtie ${head_dir}/${sp1}/hisat_aligned/$m -o ${head_dir}/${sp1}/stringtie/${stringtie_r1_name}.gtf -A -C -G ${reference_sp1}/${reference_gff3_sp1} -l ${sp1}_transcripts
		done
	
	echo "species one, two shared done"
else
	echo "species one, two shared skipped"
fi

##########################################################################################################################################

if [[ ${share_genome} == 2 && (${species_number} == 2 || ${species_number} == 3 || ${species_number} == 4 || ${species_number} == 5) ]]
then

	###########################################################################
	###################### AUTONOMOUS STEPS ##############################
	###########################################################################

	##Trim_Galore step
	
	##sp2
	mkdir ${head_dir}/${sp2}/trimmed/
	cd /${pre_trim_sp2}/

	module load trimgalore/0.6.4

	for i in *.fastq
		do
			if [[ $i == *"_R1"* || $i == *"_R2"* ]] 	
        	then
        		if [[ $i == *"_R1"* ]]
        		then
        			trim_name=$(echo $i | awk -F _ '{print $1}')
        			echo ${trim_name}
        			trim_galore --fastqc --paired -o ${head_dir}/${sp2}/trimmed ${pre_trim_sp2}/${trim_name}_R1.fastq ${pre_trim_sp2}/${trim_name}_R2.fastq
				else
					:
				fi
			else
				trim_name=$(echo $i | awk -F ".fastq" '{print $1}')
				echo ${trim_name}
				trim_galore --fastqc -o ${head_dir}/${sp2}/trimmed ${pre_trim_sp2}/$i
			fi
		done
	
	cd ${head_dir}/${sp2}/trimmed
	
	for j in *.fq
        do
        	detrim_name=$(echo $j | awk -F "_trimmed" '{print $1}')
        	echo ${detrim_name}
        	if [[ $j == *"trimmed"* ]]
            then
            	mv ${detrim_name}_trimmed.fq ${detrim_name}.fq
            else
            	:
            fi
    	 done


	for j in *.fq
        do
        	change_name=$(echo $j | awk -F _ '{print $1}')
            echo ${change_name}
            if [[ $j == *"_val_1"* || $j == *"_val_2"* ]]
            then
                if [[ $j == *"_val_1"* ]]
                then
                    mv ${change_name}_R1_val_1.fq ${change_name}_R1.fq
                else
                    mv ${change_name}_R2_val_2.fq ${change_name}_R2.fq
                fi
            else
                :
            fi
        done
		
	##Hisat2 alignment step
	mkdir ${head_dir}/${sp2}/hisat_aligned/
	cd ${head_dir}/${sp2}/trimmed/

	module load gcc/7.1.0
	module load hisat2/2.1.0
	
	for j in *.fq
		do
	  		align_name=$(echo $j | awk -F _ '{print $1}')
	        echo ${align_name}
	       
	        if [[ $j == *"_R1"* || $j == *"_R2"* ]]
	        then 
	        	if [[ "$j" == *"_R1"* ]]
	        	then	
	        		hisat2 -q -x ${reference_sp2}/${indx_sp2_directory}/${indx_name_sp2} -S ${head_dir}/${sp2}/hisat_aligned/${align_name}.sam -1 ${head_dir}/${sp2}/trimmed/${align_name}_R1.fq -2 ${head_dir}/${sp2}/trimmed/${align_name}_R2.fq
	        	else
	        		:
	        	fi
	        else
	        	hisat2 -q -x ${reference_sp2}/${indx_sp2_directory}/${indx_name_sp2} -S ${head_dir}/${sp2}/hisat_aligned/${align_name}.sam -U ${head_dir}/${sp2}/trimmed/${align_name}.fq 
			fi
		done
		
	##convert to bam and sort

	cd ${head_dir}/${sp2}/hisat_aligned/

	module load samtools/1.10

	for k in *.sam
		do
			convert2bam_name=$(echo $k | awk -F ".sam" '{print $1}')
			echo ${convert2bam_name}
			samtools view -S -b -h $k > ${convert2bam_name}.bam 
		done
	
	for k in *.bam
		do
			sort_name=$(echo $k | awk -F ".bam" '{print $1}')
			echo ${sort_name}	
			samtools sort ${sort_name}.bam -o ${sort_name}_sorted.bam
		done

	## Stringtie Transcript Creation

	cd ${head_dir}/${sp2}/hisat_aligned/
	mkdir ${head_dir}/${sp2}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for m in *_sorted.bam
		do
	        stringtie_r1_name=$(echo $m | awk -F "_sorted.bam" '{print $1}')
	        echo ${stringtie_r1_name}
	        stringtie ${head_dir}/${sp2}/hisat_aligned/$m -o ${head_dir}/${sp2}/stringtie/${stringtie_r1_name}.gtf -A -C -G ${reference_sp2}/${reference_gff3_sp2} -l ${sp2}_transcripts
		done
		
	mkdir ${head_dir}/${sp1}_${sp2}/
	mkdir ${head_dir}/${sp1}_${sp2}/stringtie/
	
	cp -r ${head_dir}/${sp1}/stringtie/ ${head_dir}/${sp1}_${sp2}/
	cp -r ${head_dir}/${sp2}/stringtie/ ${head_dir}/${sp1}_${sp2}/
	
	cd ${head_dir}/${sp1}_${sp2}/stringtie/
	
	module load gcc/7.1.0
	module load stringtie/2.0.6

	stringtie ${head_dir}/${sp1}_${sp2}/stringtie/* --merge -o ${head_dir}/${sp1}_${sp2}/stringtie/merged.gtf -G ${reference_sp2}/${reference_gff3_sp2} -l merged_transcripts
	
	mkdir ${head_dir}/${sp1}/stringtie/stringtie_final_output/
	cd ${head_dir}/${sp1}/hisat_aligned/
	
	module load gcc/7.1.0
	module load stringtie/2.0.6
	
	for n in *_sorted.bam
		do
			stringtie_r3_name=$(echo $n | awk -F "_sorted.bam" '{print $1}')
			echo ${stringtie_r3_name}
			stringtie ${head_dir}/${sp1}/hisat_aligned/$n -o ${head_dir}/${sp1}/stringtie/stringtie_final_output/${stringtie_r3_name}.gtf -A -C -G ${head_dir}/${sp1}_${sp2}/stringtie/merged.gtf
		done
	
	mkdir ${head_dir}/${sp2}/stringtie/stringtie_final_output/
	cd ${head_dir}/${sp2}/hisat_aligned/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for n in *_sorted.bam
		do
 	       stringtie_r3_name=$(echo $n | awk -F "_sorted.bam" '{print $1}')
 	       echo ${stringtie_r3_name}
 	       stringtie ${head_dir}/${sp2}/hisat_aligned/$n -o ${head_dir}/${sp2}/stringtie/stringtie_final_output/${stringtie_r3_name}.gtf -A -C -G ${head_dir}/${sp1}_${sp2}/stringtie/merged.gtf
		done

	## Converting transcripts to fasta
	export PATH=$PATH:/scratch/dtn2an/gffread
	
	cd ${head_dir}/${sp1}/stringtie/stringtie_final_output/
	mkdir ${head_dir}/${sp1}/fastasfromtranscriptome/
	
	for o in *.gtf
		do
			getfasta_name=$(echo $o | awk -F ".gtf" '{print $1}')
			echo ${getfasta_name}
			gffread -w ${head_dir}/${sp1}/fastasfromtranscriptome/${getfasta_name}_transcripts.fa -g ${reference_sp1}/${reference_genome_fasta_sp1} ${head_dir}/${sp1}/stringtie/stringtie_final_output/$o
		done
	
	export PATH=$PATH:/scratch/dtn2an/gffread

	cd ${head_dir}/${sp2}/stringtie/stringtie_final_output/
	mkdir ${head_dir}/${sp2}/fastasfromtranscriptome/
	
	for o in *.gtf
		do
 	 		getfasta_name=$(echo $o | awk -F ".gtf" '{print $1}')
 	        echo ${getfasta_name}
 	        gffread -w ${head_dir}/${sp2}/fastasfromtranscriptome/${getfasta_name}_transcripts.fa -g ${reference_sp2}/${reference_genome_fasta_sp2} ${head_dir}/${sp2}/stringtie/stringtie_final_output/$o
		done

	##TransDecoder Step

	mkdir ${head_dir}/${sp1}/TransDecoder_${sp1}_output/
	cd ${head_dir}/${sp1}/fastasfromtranscriptome/

	export PATH=$PATH:/scratch/dtn2an/TransDecoder

	for p in *.fa
		do
       	 	 name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		 echo $name
       	 	 TransDecoder.LongOrfs -t ${head_dir}/${sp1}/fastasfromtranscriptome/$p -O ${head_dir}/${sp1}/TransDecoder_${sp1}_output/
		done

	cd ${head_dir}/${sp1}/TransDecoder_${sp1}_output/

	for q in *.fa
		do
  		 	name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       	 	echo $name
       	 	TransDecoder.Predict -t ${head_dir}/${sp1}/fastasfromtranscriptome/$q -O ${head_dir}/${sp1}/TransDecoder_${sp1}_output/
		done
		
	mkdir ${head_dir}/${sp2}/TransDecoder_${sp2}_output/
	cd ${head_dir}/${sp2}/fastasfromtranscriptome/

	export PATH=$PATH:/scratch/dtn2an/TransDecoder

	for p in *.fa
		do
       	 	name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		 echo $name
       	 	TransDecoder.LongOrfs -t ${head_dir}/${sp2}/fastasfromtranscriptome/$p -O ${head_dir}/${sp2}/TransDecoder_${sp2}_output/
		done

	cd ${head_dir}/${sp2}/TransDecoder_${sp2}_output/

	for q in *.fa
		do
  		 name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       	 echo $name
       	 TransDecoder.Predict -t ${head_dir}/${sp2}/fastasfromtranscriptome/$q -O ${head_dir}/${sp2}/TransDecoder_${sp2}_output/
		done
		
	mkdir ${head_dir}/${sp1}/CD_HIT_output
	
	export PATH=$PATH:/scratch/dtn2an/cd-hit-v4.8.1-2019-0228

	./cd-hit -i ${head_dir}/${sp1}/TransDecoder_${sp1}_output/longest_orfs.pep -o ${head_dir}/${sp1}/CD_HIT_output/cdhit_longest_orfs.pep -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp1}/TransDecoder_${sp1}_output/longest_orfs.gff3 -o ${head_dir}/${sp1}/CD_HIT_output/cdhit_longest_orfs.gff3 -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp1}/TransDecoder_${sp1}_output/longest_orfs.cds -o ${head_dir}/${sp1}/CD_HIT_output/cdhit_longest_orfs.cds -c 0.95 -d 0

	mkdir ${head_dir}/${sp2}/CD_HIT_output
	
	export PATH=$PATH:/scratch/dtn2an/cd-hit-v4.8.1-2019-0228

	./cd-hit -i ${head_dir}/${sp2}/TransDecoder_${sp2}_output/longest_orfs.pep -o ${head_dir}/${sp2}/CD_HIT_output/cdhit_longest_orfs.pep -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp2}/TransDecoder_${sp2}_output/longest_orfs.gff3 -o ${head_dir}/${sp2}/CD_HIT_output/cdhit_longest_orfs.gff3 -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp2}/TransDecoder_${sp2}_output/longest_orfs.cds -o ${head_dir}/${sp2}/CD_HIT_output/cdhit_longest_orfs.cds -c 0.95 -d 0
	
	echo "species two, two shared done"
else
	echo "species two, two shared skipped"
fi

##########################################################################################################################################

if [[ ${share_genome} == 2 && (${species_number} == 3 || ${species_number} == 4 || ${species_number} == 5) ]]
then

	###########################################################################
	###################### AUTONOMOUS STEPS ##############################
	###########################################################################

	##Trim_Galore step
	
	##sp3
	mkdir ${head_dir}/${sp3}/trimmed/
	cd /${pre_trim_sp3}/

	module load trimgalore/0.6.4

	for i in *.fastq
		do
			if [[ $i == *"_R1"* || $i == *"_R2"* ]] 	
        	then
        		if [[ $i == *"_R1"* ]]
        		then
        			trim_name=$(echo $i | awk -F _ '{print $1}')
        			echo ${trim_name}
        			trim_galore --fastqc --paired -o ${head_dir}/${sp3}/trimmed ${pre_trim_sp3}/${trim_name}_R1.fastq ${pre_trim_sp3}/${trim_name}_R2.fastq
				else
					:
				fi
			else
				trim_name=$(echo $i | awk -F ".fastq" '{print $1}')
				echo ${trim_name}
				trim_galore --fastqc -o ${head_dir}/${sp3}/trimmed ${pre_trim_sp3}/$i
			fi
		done
	
	cd ${head_dir}/${sp3}/trimmed
	
	for j in *.fq
        do
        	detrim_name=$(echo $j | awk -F "_trimmed" '{print $1}')
        	echo ${detrim_name}
        	if [[ $j == *"trimmed"* ]]
            then
            	mv ${detrim_name}_trimmed.fq ${detrim_name}.fq
            else
            	:
            fi
    	 done

	for j in *.fq
        do
        	change_name=$(echo $j | awk -F _ '{print $1}')
            echo ${change_name}
            if [[ $j == *"_val_1"* || $j == *"_val_2"* ]]
            then
                if [[ $j == *"_val_1"* ]]
                then
                    mv ${change_name}_R1_val_1.fq ${change_name}_R1.fq
                else
                    mv ${change_name}_R2_val_2.fq ${change_name}_R2.fq
                fi
            else
                :
            fi
        done	
		
	##Hisat2 alignment step
	mkdir ${head_dir}/${sp3}/hisat_aligned/
	cd ${head_dir}/${sp3}/trimmed/

	module load gcc/7.1.0
	module load hisat2/2.1.0
	
	for j in *.fq
		do
	  		align_name=$(echo $j | awk -F _ '{print $1}')
	        echo ${align_name}
	       
	        if [[ $j == *"_R1"* || $j == *"_R2"* ]]
	        then 
	        	if [[ $j == *"_R1"* ]]
	        	then	
	        		hisat2 -q -x ${reference_sp3}/${indx_sp3_directory}/${indx_name_sp3} -S ${head_dir}/${sp3}/hisat_aligned/${align_name}.sam -1 ${head_dir}/${sp3}/trimmed/${align_name}_R1.fq -2 ${head_dir}/${sp3}/trimmed/${align_name}_R2.fq
	        	else
	        		:
	        	fi
	        else
	        	hisat2 -q -x ${reference_sp3}/${indx_sp3_directory}/${indx_name_sp3} -S ${head_dir}/${sp3}/hisat_aligned/${align_name}.sam -U ${head_dir}/${sp3}/trimmed/${align_name}.fq 
			fi
		done
		
	##convert to bam and sort

	cd ${head_dir}/${sp3}/hisat_aligned/

	module load samtools/1.10

	for k in *.sam
		do
			convert2bam_name=$(echo $k | awk -F ".sam" '{print $1}')
			echo ${convert2bam_name}
			samtools view -S -b -h $k > ${convert2bam_name}.bam 
		done
	
	for k in *.bam
		do
			sort_name=$(echo $k | awk -F ".bam" '{print $1}')
			echo ${sort_name}	
			samtools sort ${sort_name}.bam -o ${sort_name}_sorted.bam
		done
		
	## Stringtie Transcript Creation

	cd ${head_dir}/${sp3}/hisat_aligned/
	mkdir ${head_dir}/${sp3}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for m in *_sorted.bam
		do
	        stringtie_r1_name=$(echo $m | awk -F "_sorted.bam" '{print $1}')
	        echo ${stringtie_r1_name}
	        stringtie ${head_dir}/${sp3}/hisat_aligned/$m -o ${head_dir}/${sp3}/stringtie/${stringtie_r1_name}.gtf -A -C -G ${reference_sp3}/${reference_gff3_sp3} -l ${sp3}_transcripts
		done

	cd ${head_dir}/${sp3}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	stringtie ${head_dir}/${sp3}/stringtie/* --merge -o merged.gtf -G ${reference_sp3}/${reference_gff3_sp3} -l merged_transcripts

	mkdir ${head_dir}/${sp3}/stringtie/stringtie_final_output/
	cd ${head_dir}/${sp3}/hisat_aligned/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for n in *_sorted.bam
		do
 	   	   stringtie_r3_name=$(echo $n | awk -F "_sorted.bam" '{print $1}')
 	       echo ${stringtie_r3_name}
 	       stringtie ${head_dir}/${sp3}/hisat_aligned/$n -o ${head_dir}/${sp3}/stringtie/stringtie_final_output/${stringtie_r3_name}.gtf -A -C -G ${head_dir}/${sp3}/stringtie/merged.gtf
		done

	## Converting transcripts to fasta

	export PATH=$PATH:/scratch/dtn2an/gffread

	cd ${head_dir}/${sp3}/stringtie/stringtie_final_output/
	mkdir ${head_dir}/${sp3}/fastasfromtranscriptome/
	
	for o in *.gtf
		do
 	 		getfasta_name=$(echo $o | awk -F ".gtf" '{print $1}')
 	        echo ${getfasta_name}
 	        gffread -w ${head_dir}/${sp3}/fastasfromtranscriptome/${getfasta_name}_transcripts.fa -g ${reference_sp3}/${reference_genome_fasta_sp3} ${head_dir}/${sp3}/stringtie/stringtie_final_output/$o
		done
		
	##TransDecoder Step

	mkdir ${head_dir}/${sp3}/TransDecoder_${sp3}_output/
	cd ${head_dir}/${sp3}/fastasfromtranscriptome/

	export PATH=$PATH:/scratch/dtn2an/TransDecoder

	for p in *.fa
		do
       	 	name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		 echo $name
       	 	TransDecoder.LongOrfs -t ${head_dir}/${sp3}/fastasfromtranscriptome/$p -O ${head_dir}/${sp3}/TransDecoder_${sp3}_output/
		done

	cd ${head_dir}/${sp3}/TransDecoder_${sp3}_output/

	for q in *.fa
		do
  		 	name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       	 	echo $name
       		TransDecoder.Predict -t ${head_dir}/${sp3}/fastasfromtranscriptome/$q -O ${head_dir}/${sp3}/TransDecoder_${sp3}_output/
		done

	mkdir ${head_dir}/${sp3}/CD_HIT_output
	
	export PATH=$PATH:/scratch/dtn2an/cd-hit-v4.8.1-2019-0228

	./cd-hit -i ${head_dir}/${sp3}/TransDecoder_${sp3}_output/longest_orfs.pep -o ${head_dir}/${sp3}/CD_HIT_output/cdhit_longest_orfs.pep -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp3}/TransDecoder_${sp3}_output/longest_orfs.gff3 -o ${head_dir}/${sp3}/CD_HIT_output/cdhit_longest_orfs.gff3 -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp3}/TransDecoder_${sp3}_output/longest_orfs.cds -o ${head_dir}/${sp3}/CD_HIT_output/cdhit_longest_orfs.cds -c 0.95 -d 0
	
	echo "species three, two shared done"
else
	echo "species three, two shared skipped"				
fi

##########################################################################################################################################

if [[ ${share_genome} == 2 && (${species_number} == 4 || ${species_number} == 5) ]]
then

	###########################################################################
	###################### AUTONOMOUS STEPS ##############################
	###########################################################################

	##Trim_Galore step
	
	##sp4
	mkdir ${head_dir}/${sp4}/trimmed/
	cd /${pre_trim_sp4}/

	module load trimgalore/0.6.4

	for i in *.fastq
		do
			if [[ $i == *"_R1"* || $i == *"_R2"* ]] 	
        	then
        		if [[ $i == *"_R1"* ]]
        		then
        			trim_name=$(echo $i | awk -F _ '{print $1}')
        			echo ${trim_name}
        			trim_galore --fastqc --paired -o ${head_dir}/${sp4}/trimmed ${pre_trim_sp4}/${trim_name}_R1.fastq ${pre_trim_sp4}/${trim_name}_R2.fastq
				else
					:
				fi
			else
				trim_name=$(echo $i | awk -F ".fastq" '{print $1}')
				echo ${trim_name}
				trim_galore --fastqc -o ${head_dir}/${sp4}/trimmed ${pre_trim_sp4}/$i
			fi
		done
	
	cd ${head_dir}/${sp4}/trimmed
	
	for j in *.fq
        do
        	detrim_name=$(echo $j | awk -F "_trimmed" '{print $1}')
        	echo ${detrim_name}
        	if [[ $j == *"trimmed"* ]]
            then
            	mv ${detrim_name}_trimmed.fq ${detrim_name}.fq
            else
            	:
            fi
    	 done


	for j in *.fq
        do
        	change_name=$(echo $j | awk -F _ '{print $1}')
            echo ${change_name}
            if [[ $j == *"_val_1"* || $j == *"_val_2"* ]]
            then
                if [[ $j == *"_val_1"* ]]
                then
                    mv ${change_name}_R1_val_1.fq ${change_name}_R1.fq
                else
                    mv ${change_name}_R2_val_2.fq ${change_name}_R2.fq
                fi
            else
                :
            fi
        done
		
	##Hisat2 alignment step
	mkdir ${head_dir}/${sp4}/hisat_aligned/
	cd ${head_dir}/${sp4}/trimmed/

	module load gcc/7.1.0
	module load hisat2/2.1.0

	for j in *.fq
		do
	  		align_name=$(echo $j | awk -F _ '{print $1}')
	        echo ${align_name}
	       
	        if [[ $j == *"_R1"* || $j == *"_R2"* ]]
	        then 	
	        	if [[ $j == *"_R1"* ]]
	        	then
	        		hisat2 -q -x ${reference_sp4}/${indx_sp4_directory}/${indx_name_sp4} -S ${head_dir}/${sp4}/hisat_aligned/${align_name}.sam -1 ${head_dir}/${sp4}/trimmed/${align_name}_R1.fq -2 ${head_dir}/${sp4}/trimmed/${align_name}_R2.fq
				else
					:
				fi	        
	        else
	        	hisat2 -q -x ${reference_sp4}/${indx_sp4_directory}/${indx_name_sp4} -S ${head_dir}/${sp4}/hisat_aligned/${align_name}.sam -U ${head_dir}/${sp4}/trimmed/${align_name}.fq 
			fi
		done
	
	##convert to bam and sort

	cd ${head_dir}/${sp4}/hisat_aligned/

	module load samtools/1.10

	for k in *.sam
		do
			convert2bam_name=$(echo $k | awk -F ".sam" '{print $1}')
			echo ${convert2bam_name}
			samtools view -S -b -h $k > ${convert2bam_name}.bam 
		done
	
	for k in *.bam
		do
			sort_name=$(echo $k | awk -F ".bam" '{print $1}')
			echo ${sort_name}	
			samtools sort ${sort_name}.bam -o ${sort_name}_sorted.bam
		done

	## Stringtie Transcript Creation

	cd ${head_dir}/${sp4}/hisat_aligned/
	mkdir ${head_dir}/${sp4}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for m in *_sorted.bam
		do
	        stringtie_r1_name=$(echo $m | awk -F "_sorted.bam" '{print $1}')
	        echo ${stringtie_r1_name}
	        stringtie ${head_dir}/${sp4}/hisat_aligned/$m -o ${head_dir}/${sp4}/stringtie/${stringtie_r1_name}.gtf -A -C -G ${reference_sp4}/${reference_gff3_sp4} -l ${sp4}_transcripts
		done

	cd ${head_dir}/${sp4}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	stringtie ${head_dir}/${sp4}/stringtie/* --merge -o merged.gtf -G ${reference_sp4}/${reference_gff3_sp4} -l merged_transcripts

	mkdir ${head_dir}/${sp4}/stringtie/stringtie_final_output/
	cd ${head_dir}/${sp4}/hisat_aligned/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for n in *_sorted.bam
		do
 	       stringtie_r3_name=$(echo $n | awk -F "_sorted.bam" '{print $1}')
 	       echo ${stringtie_r3_name}
 	       stringtie ${head_dir}/${sp4}/hisat_aligned/$n -o ${head_dir}/${sp4}/stringtie/stringtie_final_output/${stringtie_r3_name}.gtf -A -C -G ${head_dir}/${sp4}/stringtie/merged.gtf
		done

	## Converting transcripts to fasta

	export PATH=$PATH:/scratch/dtn2an/gffread

	cd ${head_dir}/${sp4}/stringtie/stringtie_final_output/
	mkdir ${head_dir}/${sp4}/fastasfromtranscriptome/
	
	for o in *.gtf
		do
 	 		getfasta_name=$(echo $o | awk -F ".gtf" '{print $1}')
 	        echo ${getfasta_name}
 	        gffread -w ${head_dir}/${sp4}/${getfasta_name}_transcripts.fa -g ${reference_sp4}/${reference_genome_fasta_sp4} ${head_dir}/${sp4}/stringtie/stringtie_final_output/$o
		done
	
	##TransDecoder Step

	mkdir ${head_dir}/${sp4}/TransDecoder_${sp4}_output/
	cd ${head_dir}/${sp4}/fastasfromtranscriptome/

	export PATH=$PATH:/scratch/dtn2an/TransDecoder

	for p in *.fa
		do
       	 	name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		 echo $name
       	 	TransDecoder.LongOrfs -t ${head_dir}/${sp4}/fastasfromtranscriptome/$p -O ${head_dir}/${sp4}/TransDecoder_${sp4}_output/
		done

	cd ${head_dir}/${sp4}/TransDecoder_${sp4}_output/

	for q in *.fa
		do
  			 name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		 echo $name
       		 TransDecoder.Predict -t ${head_dir}/${sp4}/fastasfromtranscriptome/$q -O ${head_dir}/${sp4}/TransDecoder_${sp4}_output/
		done
	
	mkdir ${head_dir}/${sp4}/CD_HIT_output
	
	export PATH=$PATH:/scratch/dtn2an/cd-hit-v4.8.1-2019-0228

	./cd-hit -i ${head_dir}/${sp4}/TransDecoder_${sp4}_output/longest_orfs.pep -o ${head_dir}/${sp4}/CD_HIT_output/cdhit_longest_orfs.pep -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp4}/TransDecoder_${sp4}_output/longest_orfs.gff3 -o ${head_dir}/${sp4}/CD_HIT_output/cdhit_longest_orfs.gff3 -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp4}/TransDecoder_${sp4}_output/longest_orfs.cds -o ${head_dir}/${sp4}/CD_HIT_output/cdhit_longest_orfs.cds -c 0.95 -d 0

	echo "species four, two shared done"
else
	echo "species four, two shared skipped"
fi

##########################################################################################################################################

if [[ ${share_genome} == 2 && (${species_number} == 5) ]]
then

	###########################################################################
	###################### AUTONOMOUS STEPS ##############################
	###########################################################################

	##Trim_Galore step
	
	##sp5
	mkdir ${head_dir}/${sp5}/trimmed/
	cd /${pre_trim_sp5}/

	module load trimgalore/0.6.4

	for i in *.fastq
		do
			if [[ $i == *"_R1"* || $i == *"_R2"* ]] 	
        	then
        		if [[ $i == *"_R1"* ]]
        		then
        			trim_name=$(echo $i | awk -F _ '{print $1}')
        			echo ${trim_name}
        			trim_galore --fastqc --paired -o ${head_dir}/${sp5}/trimmed ${pre_trim_sp5}/${trim_name}_R1.fastq ${pre_trim_sp5}/${trim_name}_R2.fastq
				else
					:
				fi
			else
				trim_name=$(echo $i | awk -F ".fastq" '{print $1}')
				echo ${trim_name}
				trim_galore --fastqc -o ${head_dir}/${sp5}/trimmed ${pre_trim_sp5}/$i
			fi
		done
	
	cd ${head_dir}/${sp5}/trimmed
	
	for j in *.fq
        do
        	detrim_name=$(echo $j | awk -F "_trimmed" '{print $1}')
        	echo ${detrim_name}
        	if [[ $j == *"trimmed"* ]]
            then
            	mv ${detrim_name}_trimmed.fq ${detrim_name}.fq
            else
            	:
            fi
    	 done


	for j in *.fq
        do
        	change_name=$(echo $j | awk -F _ '{print $1}')
            echo ${change_name}
            if [[ $j == *"_val_1"* || $j == *"_val_2"* ]]
            then
                if [[ $j == *"_val_1"* ]]
                then
                    mv ${change_name}_R1_val_1.fq ${change_name}_R1.fq
                else
                    mv ${change_name}_R2_val_2.fq ${change_name}_R2.fq
                fi
            else
                :
            fi
        done
		
	##Hisat2 alignment step
	mkdir ${head_dir}/${sp5}/hisat_aligned/
	cd ${head_dir}/${sp5}/trimmed/

	module load gcc/7.1.0
	module load hisat2/2.1.0

	for j in *.fq
		do
	  		align_name=$(echo $j | awk -F _ '{print $1}')
	        echo ${align_name}
	       
	        if [[ $j == *"_R1"* || $j == *"_R2"* ]]
	        then 	
	        	if [[ $j == *"_R1"* ]]
	        	then
	        		hisat2 -q -x ${reference_sp5}/${indx_sp5_directory}/${indx_name_sp5} -S ${head_dir}/${sp5}/hisat_aligned/${align_name}.sam -1 ${head_dir}/${sp5}/trimmed/${align_name}_R1.fq -2 ${head_dir}/${sp5}/trimmed/${align_name}_R2.fq
	        	else
	        		:
	        	fi
	        else
	        	hisat2 -q -x ${reference_sp5}/${indx_sp5_directory}/${indx_name_sp5} -S ${head_dir}/${sp5}/hisat_aligned/${align_name}.sam -U ${head_dir}/${sp5}/trimmed/${align_name}.fq 
			fi
		done


	##convert to bam and sort

	cd ${head_dir}/${sp5}/hisat_aligned/

	module load samtools/1.10

	for k in *.sam
		do
			convert2bam_name=$(echo $k | awk -F ".sam" '{print $1}')
			echo ${convert2bam_name}
			samtools view -S -b -h $k > ${convert2bam_name}.bam 
		done
	
	for k in *.bam
		do
			sort_name=$(echo $k | awk -F ".bam" '{print $1}')
			echo ${sort_name}	
			samtools sort ${sort_name}.bam -o ${sort_name}_sorted.bam
		done
		
	## Stringtie Transcript Creation

	cd ${head_dir}/${sp5}/hisat_aligned/
	mkdir ${head_dir}/${sp5}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for m in *_sorted.bam
		do

	        stringtie_r1_name=$(echo $m | awk -F "_sorted.bam" '{print $1}')
	        echo ${stringtie_r1_name}
	        stringtie ${head_dir}/${sp5}/hisat_aligned/$m -o ${head_dir}/${sp5}/stringtie/${stringtie_r1_name}.gtf -A -C -G ${reference_sp5}/${reference_gff3_sp5} -l ${sp5}_transcripts
		done

	cd ${head_dir}/${sp5}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	stringtie ${head_dir}/${sp5}/stringtie/* --merge -o merged.gtf -G ${reference_sp5}/${reference_gff3_sp5} -l merged_transcripts

	mkdir ${head_dir}/${sp5}/stringtie/stringtie_final_output/
	cd ${head_dir}/${sp5}/hisat_aligned/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for n in *_sorted.bam
		do
 	       stringtie_r3_name=$(echo $n | awk -F "_sorted.bam" '{print $1}')
 	       echo ${stringtie_r3_name}
 	       stringtie ${head_dir}/${sp5}/hisat_aligned/$n -o ${head_dir}/${sp5}/stringtie/stringtie_final_output/${stringtie_r3_name}.gtf -A -C -G ${head_dir}/${sp5}/stringtie/merged.gtf
		done

	## Converting transcripts to fasta

	export PATH=$PATH:/scratch/dtn2an/gffread

	cd ${head_dir}/${sp5}/stringtie/stringtie_final_output/
	mkdir ${head_dir}/${sp5}/fastasfromtranscriptome/
	
	for o in *.gtf
		do
 	 		getfasta_name=$(echo $o | awk -F ".gtf" '{print $1}')
 	        echo ${getfasta_name}
 	        gffread -w ${head_dir}/${sp5}/fastasfromtranscriptome/${getfasta_name}_transcripts.fa -g ${reference}/${reference_genome_fasta} ${head_dir}/${sp5}/stringtie/stringtie_final_output/$o
		done
		
	##TransDecoder Step

	mkdir ${head_dir}/${sp5}/TransDecoder_${sp5}_output/
	cd ${head_dir}/${sp5}/fastasfromtranscriptome/

	export PATH=$PATH:/scratch/dtn2an/TransDecoder

	for p in *.fa
		do
       	 	 name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		 echo $name
       	 	 TransDecoder.LongOrfs -t ${head_dir}/${sp5}/fastasfromtranscriptome/$p -O ${head_dir}/${sp5}/TransDecoder_${sp5}_output/
		done

	cd ${head_dir}/${sp5}/TransDecoder_${sp5}_output/

	for q in *.fa
		do
  		 	 name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		 echo $name
       		 TransDecoder.Predict -t ${head_dir}/${sp5}/fastasfromtranscriptome/$q -O ${head_dir}/${sp5}/TransDecoder_${sp5}_output/
		done
	
	mkdir ${head_dir}/${sp5}/CD_HIT_output
	
	export PATH=$PATH:/scratch/dtn2an/cd-hit-v4.8.1-2019-0228

	./cd-hit -i ${head_dir}/${sp5}/TransDecoder_${sp5}_output/longest_orfs.pep -o ${head_dir}/${sp5}/CD_HIT_output/cdhit_longest_orfs.pep -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp5}/TransDecoder_${sp5}_output/longest_orfs.gff3 -o ${head_dir}/${sp5}/CD_HIT_output/cdhit_longest_orfs.gff3 -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp5}/TransDecoder_${sp5}_output/longest_orfs.cds -o ${head_dir}/${sp5}/CD_HIT_output/cdhit_longest_orfs.cds -c 0.95 -d 0

	echo "species five, two shared done"
	
else
	echo "species five, two shared skipped"
fi

##########################################################################################################################################
##########################################################################################################################################

if [[ ${share_genome} == 3 && (${species_number} == 3 || ${species_number} == 4 || ${species_number} == 5) ]]
then


	###########################################################################
	###################### AUTONOMOUS STEPS ##############################
	###########################################################################

	##Trim_Galore step
	
	##sp1
	mkdir ${head_dir}/${sp1}/trimmed/
	cd ${pre_trim_sp1}/

	module load trimgalore/0.6.4

	for i in *.fastq
		do
			if [[ $i == *"_R1"* || $i == *"_R2"* ]] 	
        	then
        		if [[ $i == *"_R1"* ]]
        		then
        			trim_name=$(echo $i | awk -F _ '{print $1}')
        			echo ${trim_name}
        			trim_galore --fastqc --paired -o ${head_dir}/${sp1}/trimmed ${pre_trim_sp1}/${trim_name}_R1.fastq ${pre_trim_sp1}/${trim_name}_R2.fastq
				else
					:
				fi
			else
				trim_name=$(echo $i | awk -F ".fastq" '{print $1}')
				echo ${trim_name}
				trim_galore --fastqc -o ${head_dir}/${sp1}/trimmed ${pre_trim_sp1}/$i
			fi
		done
	
	cd ${head_dir}/${sp1}/trimmed
	
	for j in *.fq
        do
        	detrim_name=$(echo $j | awk -F "_trimmed" '{print $1}')
        	echo ${detrim_name}
        	if [[ $j == *"trimmed"* ]]
            then
            	mv ${detrim_name}_trimmed.fq ${detrim_name}.fq
            else
            	:
            fi
    	 done


	for j in *.fq
        do
        	change_name=$(echo $j | awk -F _ '{print $1}')
            echo ${change_name}
            if [[ $j == *"_val_1"* || $j == *"_val_2"* ]]
            then
                if [[ $j == *"_val_1"* ]]
                then
                    mv ${change_name}_R1_val_1.fq ${change_name}_R1.fq
                else
                    mv ${change_name}_R2_val_2.fq ${change_name}_R2.fq
                fi
            else
                :
            fi
        done

	##Hisat2 alignment step
	mkdir ${head_dir}/${sp1}/hisat_aligned/
	cd ${head_dir}/${sp1}/trimmed/

	module load gcc/7.1.0
	module load hisat2/2.1.0

	for j in *.fq
		do
	  		align_name=$(echo $j | awk -F _ '{print $1}')
	        echo ${align_name}
	       
	        if [[ $j == *"_R1"* || $j == *"_R2"* ]]
	        then 	
	        	if [[ $j == *"_R1"* ]]
	        	then
	        		hisat2 -q -x ${reference_sp1}/${indx_sp1_directory}/${indx_name_sp1} -S ${head_dir}/${sp1}/hisat_aligned/${align_name}.sam -1 ${head_dir}/${sp1}/trimmed/${align_name}_R1.fq -2 ${head_dir}/${sp1}/trimmed/${align_name}_R2.fq
	        	else
	        		:
	        	fi
	        else
	        	hisat2 -q -x ${reference_sp1}/${indx_sp1_directory}/${indx_name_sp1} -S ${head_dir}/${sp1}/hisat_aligned/${align_name}.sam -U ${head_dir}/${sp1}/trimmed/${align_name}.fq 
			fi
		done
	

	##convert to bam and sort

	cd ${head_dir}/${sp1}/hisat_aligned/

	module load samtools/1.10

	for k in *.sam
		do
			convert2bam_name=$(echo $k | awk -F ".sam" '{print $1}')
			echo ${convert2bam_name}
			samtools view -S -b -h $k > ${convert2bam_name}.bam 
		done
	
	for k in *.bam
		do
			sort_name=$(echo $k | awk -F ".bam" '{print $1}')
			echo ${sort_name}	
			samtools sort ${sort_name}.bam -o ${sort_name}_sorted.bam
		done


	## Stringtie Transcript Creation

	cd ${head_dir}/${sp1}/hisat_aligned/
	mkdir ${head_dir}/${sp1}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for m in *_sorted.bam
		do
	        stringtie_r1_name=$(echo $m | awk -F "_sorted.bam" '{print $1}')
	        echo ${stringtie_r1_name}
	        stringtie ${head_dir}/${sp1}/hisat_aligned/$m -o ${head_dir}/${sp1}/stringtie/${stringtie_r1_name}.gtf -A -C -G ${reference_sp1}/${reference_gff3_sp1} -l ${sp1}_transcripts
		done
		
	echo "species one, three shared done"
else
	echo "species one, three shared skipped"
fi

##########################################################################################################################################

if [[ ${share_genome} == 3 && (${species_number} == 3 || ${species_number} == 4 || ${species_number} == 5) ]]
then

	###########################################################################
	###################### AUTONOMOUS STEPS ##############################
	###########################################################################

	##Trim_Galore step
	
	##sp2
	mkdir ${head_dir}/${sp2}/trimmed/
	cd /${pre_trim_sp2}/

	module load trimgalore/0.6.4

	for i in *.fastq
		do
			if [[ $i == *"_R1"* || $i == *"_R2"* ]] 	
        	then
        		if [[ $i == *"_R1"* ]]
        		then
        			trim_name=$(echo $i | awk -F _ '{print $1}')
        			echo ${trim_name}
        			trim_galore --fastqc --paired -o ${head_dir}/${sp2}/trimmed ${pre_trim_sp2}/${trim_name}_R1.fastq ${pre_trim_sp2}/${trim_name}_R2.fastq
				else
					:
				fi
			else
				trim_name=$(echo $i | awk -F ".fastq" '{print $1}')
				echo ${trim_name}
				trim_galore --fastqc -o ${head_dir}/${sp2}/trimmed ${pre_trim_sp2}/$i
			fi
		done

	cd ${head_dir}/${sp2}/trimmed

	for j in *.fq
        do
        	detrim_name=$(echo $j | awk -F "_trimmed" '{print $1}')
        	echo ${detrim_name}
        	if [[ $j == *"trimmed"* ]]
            then
            	mv ${detrim_name}_trimmed.fq ${detrim_name}.fq
            else
            	:
            fi
    	 done


	for j in *.fq
        do
        	change_name=$(echo $j | awk -F _ '{print $1}')
            echo ${change_name}
            if [[ $j == *"_val_1"* || $j == *"_val_2"* ]]
            then
                if [[ $j == *"_val_1"* ]]
                then
                    mv ${change_name}_R1_val_1.fq ${change_name}_R1.fq
                else
                    mv ${change_name}_R2_val_2.fq ${change_name}_R2.fq
                fi
            else
                :
            fi
        done

	##Hisat2 alignment step
	mkdir ${head_dir}/${sp2}/hisat_aligned/
	cd ${head_dir}/${sp2}/trimmed/

	module load gcc/7.1.0
	module load hisat2/2.1.0
	
	for j in *.fq
		do
	  		align_name=$(echo $j | awk -F _ '{print $1}')
	        echo ${align_name}
	       
	        if [[ $j == *"_R1"* || $j == *"_R2"* ]]
	        then 
	        	if [[ $j == *"_R1"* ]]
	        	then	
	        		hisat2 -q -x ${reference_sp2}/${indx_sp2_directory}/${indx_name_sp2} -S ${head_dir}/${sp2}/hisat_aligned/${align_name}.sam -1 ${head_dir}/${sp2}/trimmed/${align_name}_R1.fq -2 ${head_dir}/${sp2}/trimmed/${align_name}_R2.fq
	       		else
	       			:
	       		fi
	        else
	        	hisat2 -q -x ${reference_sp2}/${indx_sp2_directory}/${indx_name_sp2} -S ${head_dir}/${sp2}/hisat_aligned/${align_name}.sam -U ${head_dir}/${sp2}/trimmed/${align_name}.fq 
			fi
		done


	##convert to bam and sort

	cd ${head_dir}/${sp2}/hisat_aligned/

	module load samtools/1.10

	for k in *.sam
		do
			convert2bam_name=$(echo $k | awk -F ".sam" '{print $1}')
			echo ${convert2bam_name}
			samtools view -S -b -h $k > ${convert2bam_name}.bam 
		done
	
	for k in *.bam
		do
			sort_name=$(echo $k | awk -F ".bam" '{print $1}')
			echo ${sort_name}	
			samtools sort ${sort_name}.bam -o ${sort_name}_sorted.bam
		done
		
	## Stringtie Transcript Creation

	cd ${head_dir}/${sp2}/hisat_aligned/
	mkdir ${head_dir}/${sp2}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for m in *_sorted.bam
		do
	        stringtie_r1_name=$(echo $m | awk -F "_sorted.bam" '{print $1}')
	        echo ${stringtie_r1_name}
	        stringtie ${head_dir}/${sp2}/hisat_aligned/$m -o ${head_dir}/${sp2}/stringtie/${stringtie_r1_name}.gtf -A -C -G ${reference_sp2}/${reference_gff3_sp2} -l ${sp2}_transcripts
		done
	
	echo "species two, three shared done"
else
	echo "species two, three shared skipped"
fi

##########################################################################################################################################

if [[ ${share_genome} == 3 && (${species_number} == 3 || ${species_number} == 4 || ${species_number} == 5) ]]
then 

	###########################################################################
	###################### AUTONOMOUS STEPS ##############################
	###########################################################################

	##Trim_Galore step
	
	##sp3
	mkdir ${head_dir}/${sp3}/trimmed/
	cd /${pre_trim_sp3}/

	module load trimgalore/0.6.4

	for i in *.fastq
		do
			if [[ $i == *"_R1"* || $i == *"_R2"* ]] 	
        	then
        		if [[ $i == *"_R1"* ]]
        		then
        			trim_name=$(echo $i | awk -F _ '{print $1}')
        			echo ${trim_name}
        			trim_galore --fastqc --paired -o ${head_dir}/${sp3}/trimmed ${pre_trim_sp3}/${trim_name}_R1.fastq ${pre_trim_sp3}/${trim_name}_R2.fastq
				else
					:
				fi
			else
				trim_name=$(echo $i | awk -F ".fastq" '{print $1}')
				echo ${trim_name}
				trim_galore --fastqc -o ${head_dir}/${sp3}/trimmed ${pre_trim_sp3}/$i
			fi
		done
	
	cd ${head_dir}/${sp3}/trimmed

	for j in *.fq
        do
        	detrim_name=$(echo $j | awk -F "_trimmed" '{print $1}')
        	echo ${detrim_name}
        	if [[ $j == *"trimmed"* ]]
            then
            	mv ${detrim_name}_trimmed.fq ${detrim_name}.fq
            else
            	:
            fi
    	 done


	for j in *.fq
        do
        	change_name=$(echo $j | awk -F _ '{print $1}')
            echo ${change_name}
            if [[ $j == *"_val_1"* || $j == *"_val_2"* ]]
            then
                if [[ $j == *"_val_1"* ]]
                then
                    mv ${change_name}_R1_val_1.fq ${change_name}_R1.fq
                else
                    mv ${change_name}_R2_val_2.fq ${change_name}_R2.fq
                fi
            else
                :
            fi
        done

	##Hisat2 alignment step
	mkdir ${head_dir}/${sp3}/hisat_aligned/
	cd ${head_dir}/${sp3}/trimmed/

	module load gcc/7.1.0
	module load hisat2/2.1.0

	for j in *.fq
		do
	  		align_name=$(echo $j | awk -F _ '{print $1}')
	        echo ${align_name}
	       
	        if [[ $j == *"_R1"* || $j == *"_R2"* ]]
	        then 
	        	if [[ $j == *"_R1"* ]]
	        	then	
	        		hisat2 -q -x ${reference_sp3}/${indx_sp3_directory}/${indx_name_sp3} -S ${head_dir}/${sp3}/hisat_aligned/${align_name}.sam -1 ${head_dir}/${sp3}/trimmed/${align_name}_R1.fq -2 ${head_dir}/${sp3}/trimmed/${align_name}_R2.fq
	        	else
	        		:
	        	fi
	        else
	        	hisat2 -q -x ${reference_sp3}/${indx_sp3_directory}/${indx_name_sp3} -S ${head_dir}/${sp3}/hisat_aligned/${align_name}.sam -U ${head_dir}/${sp3}/trimmed/${align_name}.fq 
			fi
		done

	##convert to bam and sort

	cd ${head_dir}/${sp3}/hisat_aligned/

	module load samtools/1.10

	for k in *.sam
		do
			convert2bam_name=$(echo $k | awk -F ".sam" '{print $1}')
			echo ${convert2bam_name}
			samtools view -S -b -h $k > ${convert2bam_name}.bam 
		done
	
	for k in *.bam
		do
			sort_name=$(echo $k | awk -F ".bam" '{print $1}')
			echo ${sort_name}	
			samtools sort ${sort_name}.bam -o ${sort_name}_sorted.bam
		done
	

	## Stringtie Transcript Creation

	cd ${head_dir}/${sp3}/hisat_aligned/
	mkdir ${head_dir}/${sp3}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for m in *_sorted.bam
		do
	        stringtie_r1_name=$(echo $m | awk -F "_sorted.bam" '{print $1}')
	        echo ${stringtie_r1_name}
	        stringtie ${head_dir}/${sp3}/hisat_aligned/$m -o ${head_dir}/${sp3}/stringtie/${stringtie_r1_name}.gtf -A -C -G ${reference_sp3}/${reference_gff3_sp3} -l ${sp3}_transcripts
		done

	mkdir ${head_dir}/${sp1}_${sp2}_${sp3}/
	mkdir ${head_dir}/${sp1}_${sp2}_${sp3}/stringtie/
	
	cp -r ${head_dir}/${sp1}/stringtie/ ${head_dir}/${sp1}_${sp2}_${sp3}/
	cp -r ${head_dir}/${sp2}/stringtie/ ${head_dir}/${sp1}_${sp2}_${sp3}/
	cp -r ${head_dir}/${sp3}/stringtie/ ${head_dir}/${sp1}_${sp2}_${sp3}/
	
	cd ${head_dir}/${sp1}_${sp2}_${sp3}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	stringtie ${head_dir}/${sp1}_${sp2}_${sp3}/stringtie/* --merge -o ${head_dir}/${sp1}_${sp2}_${sp3}/merged.gtf -G ${reference_sp3}/${reference_gff3_sp3} -l merged_transcripts
	
	mkdir ${head_dir}/${sp1}/stringtie/stringtie_final_output/
	cd ${head_dir}/${sp1}/hisat_aligned/
	
	module load gcc/7.1.0
	module load stringtie/2.0.6
	
	for n in *_sorted.bam
		do
			stringtie_r3_name=$(echo $n | awk -F "_sorted.bam" '{print $1}')
			echo ${stringtie_r3_name}
			stringtie ${head_dir}/${sp1}/hisat_aligned/$n -o ${head_dir}/${sp1}/stringtie/stringtie_final_output/${stringtie_r3_name}.gtf -A -C -G ${head_dir}/${sp1}_${sp2}_${sp3}/stringtie/merged.gtf
		done
	mkdir ${head_dir}/${sp2}/stringtie/stringtie_final_output/
	cd ${head_dir}/${sp2}/hisat_aligned/
	
	module load gcc/7.1.0
	module load stringtie/2.0.6
	
	for n in *_sorted.bam
		do
			stringtie_r3_name=$(echo $n | awk -F "_sorted.bam" '{print $1}')
			echo ${stringtie_r3_name}
			stringtie ${head_dir}/${sp2}/hisat_aligned/$n -o ${head_dir}/${sp2}/stringtie/stringtie_final_output/${stringtie_r3_name}.gtf -A -C -G ${head_dir}/${sp1}_${sp2}_${sp3}/stringtie/merged.gtf
		done
	mkdir ${head_dir}/${sp3}/stringtie/stringtie_final_output/
	cd ${head_dir}/${sp3}/hisat_aligned/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for n in *_sorted.bam
		do
 	   	   stringtie_r3_name=$(echo $n | awk -F "_sorted.bam" '{print $1}')
 	       echo ${stringtie_r3_name}
 	       stringtie ${head_dir}/${sp3}/hisat_aligned/$n -o ${head_dir}/${sp3}/stringtie/stringtie_final_output/${stringtie_r3_name}.gtf -A -C -G ${head_dir}/${sp1}_${sp2}_${sp3}/stringtie/merged.gtf
		done

	## Converting transcripts to fasta
	export PATH=$PATH:/scratch/dtn2an/gffread
	
	cd ${head_dir}/${sp1}/stringtie_stringtie_final_output/
	mkdir ${head_dir}/${sp1}/fastasfromtranscriptome/
	
	for o in *.gtf
		do
			getfasta_name=$(echo $o | awk -F ".gtf" '{print $1}')
			echo ${getfasta_name}
			gffread -w ${head_dir}/${sp1}/fastasfromtranscriptome/${getfasta_name}_transcripts.fa -g ${reference_sp1}/${reference_genome_fasta_sp1} ${head_dir}/${sp1}/stringtie/stringtie_final_output/$o
		done
	export PATH=$PATH:/scratch/dtn2an/gffread
	
	cd ${head_dir}/${sp2}/stringtie/stringtie_final_output/
	mkdir ${head_dir}/${sp2}/fastasfromtranscriptome/
	
	for o in *.gtf
		do
			getfasta_name=$(echo $o | awk -F ".gtf" '{print $1}')
			echo ${getfasta_name}
			gffread -w ${head_dir}/${sp2}/fastasfromtranscriptome/${getfasta_name}_transcripts.fa -g ${reference_sp2}/${reference_genome_fasta_sp2} ${head_dir}/${sp2}/stringtie/stringtie_final_output/$o
		done	
	export PATH=$PATH:/scratch/dtn2an/gffread

	cd ${head_dir}/${sp3}/stringtie/stringtie_final_output/
	mkdir ${head_dir}/${sp3}/fastasfromtranscriptome/
	
	for o in *.gtf
		do
 	 		getfasta_name=$(echo $o | awk -F ".gtf" '{print $1}')
 	        echo ${getfasta_name}
 	        gffread -w ${head_dir}/${sp3}/fastasfromtranscriptome/${getfasta_name}_transcripts.fa -g ${reference_sp3}/${reference_genome_fasta_sp3} ${head_dir}/${sp3}/stringtie/stringtie_final_output/$o
		done
	
	
	##TransDecoder Step

	mkdir ${head_dir}/${sp1}/TransDecoder_${sp1}_output/
	cd ${head_dir}/${sp1}/fastasfromtranscriptome/

	export PATH=$PATH:/scratch/dtn2an/TransDecoder

	for p in *.fa
		do
       	 	name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		 echo $name
       	 	TransDecoder.LongOrfs -t ${head_dir}/${sp1}/fastasfromtranscriptome/$p -O ${head_dir}/${sp1}/TransDecoder_${sp1}_output/
		done

	cd ${head_dir}/${sp1}/TransDecoder_${sp1}_output/

	for q in *.fa
		do
  		 	name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       	 	echo $name
       		TransDecoder.Predict -t ${head_dir}/${sp1}/fastasfromtranscriptome/$q -O ${head_dir}/${sp1}/TransDecoder_${sp1}_output/
		done	
	
	mkdir ${head_dir}/${sp2}/TransDecoder_${sp2}_output/
	cd ${head_dir}/${sp2}/fastasfromtranscriptome/

	export PATH=$PATH:/scratch/dtn2an/TransDecoder

	for p in *.fa
		do
       	 	name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		 echo $name
       	 	TransDecoder.LongOrfs -t ${head_dir}/${sp2}/fastasfromtranscriptome/$p -O ${head_dir}/${sp2}/TransDecoder_${sp2}_output/
		done

	cd ${head_dir}/${sp2}/TransDecoder_${sp2}_output/

	for q in *.fa
		do
  		 	name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       	 	echo $name
       		TransDecoder.Predict -t ${head_dir}/${sp2}/fastasfromtranscriptome/$q -O ${head_dir}/${sp2}/TransDecoder_${sp2}_output/
		done		
				
	mkdir ${head_dir}/${sp3}/TransDecoder_${sp3}_output/
	cd ${head_dir}/${sp3}/fastasfromtranscriptome/

	export PATH=$PATH:/scratch/dtn2an/TransDecoder

	for p in *.fa
		do
       	 	name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		 echo $name
       	 	TransDecoder.LongOrfs -t ${head_dir}/${sp3}/fastasfromtranscriptome/$p -O ${head_dir}/${sp3}/TransDecoder_${sp3}_output/
		done

	cd ${head_dir}/${sp3}/TransDecoder_${sp3}_output/

	for q in *.fa
		do
  		 	name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       	 	echo $name
       		TransDecoder.Predict -t ${head_dir}/${sp3}/fastasfromtranscriptome/$q -O ${head_dir}/${sp3}/TransDecoder_${sp3}_output/
		done
	
	mkdir ${head_dir}/${sp1}/CD_HIT_output
	
	export PATH=$PATH:/scratch/dtn2an/cd-hit-v4.8.1-2019-0228

	./cd-hit -i ${head_dir}/${sp1}/TransDecoder_${sp1}_output/longest_orfs.pep -o ${head_dir}/${sp1}/CD_HIT_output/cdhit_longest_orfs.pep -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp1}/TransDecoder_${sp1}_output/longest_orfs.gff3 -o ${head_dir}/${sp1}/CD_HIT_output/cdhit_longest_orfs.gff3 -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp1}/TransDecoder_${sp1}_output/longest_orfs.cds -o ${head_dir}/${sp1}/CD_HIT_output/cdhit_longest_orfs.cds -c 0.95 -d 0

	mkdir ${head_dir}/${sp2}/CD_HIT_output
	
	export PATH=$PATH:/scratch/dtn2an/cd-hit-v4.8.1-2019-0228

	./cd-hit -i ${head_dir}/${sp2}/TransDecoder_${sp2}_output/longest_orfs.pep -o ${head_dir}/${sp2}/CD_HIT_output/cdhit_longest_orfs.pep -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp2}/TransDecoder_${sp2}_output/longest_orfs.gff3 -o ${head_dir}/${sp2}/CD_HIT_output/cdhit_longest_orfs.gff3 -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp2}/TransDecoder_${sp2}_output/longest_orfs.cds -o ${head_dir}/${sp2}/CD_HIT_output/cdhit_longest_orfs.cds -c 0.95 -d 0

	mkdir ${head_dir}/${sp3}/CD_HIT_output
	
	export PATH=$PATH:/scratch/dtn2an/cd-hit-v4.8.1-2019-0228

	./cd-hit -i ${head_dir}/${sp3}/TransDecoder_${sp3}_output/longest_orfs.pep -o ${head_dir}/${sp3}/CD_HIT_output/cdhit_longest_orfs.pep -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp3}/TransDecoder_${sp3}_output/longest_orfs.gff3 -o ${head_dir}/${sp3}/CD_HIT_output/cdhit_longest_orfs.gff3 -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp3}/TransDecoder_${sp3}_output/longest_orfs.cds -o ${head_dir}/${sp3}/CD_HIT_output/cdhit_longest_orfs.cds -c 0.95 -d 0
		
	echo "species three, three shared done"
else
	echo "species three, three shared skipped"	
fi

##########################################################################################################################################

if [[ ${share_genome} == 3 && (${species_number} == 4 || ${species_number} == 5) ]]
then

	###########################################################################
	###################### AUTONOMOUS STEPS ##############################
	###########################################################################

	##Trim_Galore step
	
	##sp4
	mkdir ${head_dir}/${sp4}/trimmed/
	cd /${pre_trim_sp4}/

	module load trimgalore/0.6.4

	for i in *.fastq
		do
			if [[ $i == *"_R1"* || $i == *"_R2"* ]] 	
        	then
        		if [[ $i == *"_R1"* ]]
        		then
        			trim_name=$(echo $i | awk -F _ '{print $1}')
        			echo ${trim_name}
        			trim_galore --fastqc --paired -o ${head_dir}/${sp4}/trimmed ${pre_trim_sp4}/${trim_name}_R1.fastq ${pre_trim_sp4}/${trim_name}_R2.fastq
				else
					:
				fi
			else
				trim_name=$(echo $i | awk -F ".fastq" '{print $1}')
				echo ${trim_name}
				trim_galore --fastqc -o ${head_dir}/${sp4}/trimmed ${pre_trim_sp4}/$i
			fi
		done

	cd ${head_dir}/${sp4}/trimmed

	for j in *.fq
        do
        	detrim_name=$(echo $j | awk -F "_trimmed" '{print $1}')
        	echo ${detrim_name}
        	if [[ $j == *"trimmed"* ]]
            then
            	mv ${detrim_name}_trimmed.fq ${detrim_name}.fq
            else
            	:
            fi
    	 done

	for j in *.fq
        do
        	change_name=$(echo $j | awk -F _ '{print $1}')
            echo ${change_name}
            if [[ $j == *"_val_1"* || $j == *"_val_2"* ]]
            then
                if [[ $j == *"_val_1"* ]]
                then
                    mv ${change_name}_R1_val_1.fq ${change_name}_R1.fq
                else
                    mv ${change_name}_R2_val_2.fq ${change_name}_R2.fq
                fi
            else
                :
            fi
        done

	##Hisat2 alignment step
	mkdir ${head_dir}/${sp4}/hisat_aligned/
	cd ${head_dir}/${sp4}/trimmed/

	module load gcc/7.1.0
	module load hisat2/2.1.0

	for j in *.fq
		do
	  		align_name=$(echo $j | awk -F _ '{print $1}')
	        echo ${align_name}
	       
	        if [[ "$j" == *"_R1"* || $j == *"_R2"* ]]
	        then 	
	        	if [[ $j == *"_R1"* ]]
	        	then
	        		hisat2 -q -x ${reference_sp4}/${indx_sp4_directory}/${indx_name_sp4} -S ${head_dir}/${sp4}/hisat_aligned/${align_name}.sam -1 ${head_dir}/${sp4}/trimmed/${align_name}_R1.fq -2 ${head_dir}/${sp4}/trimmed/${align_name}_R2.fq
	       		else
	       			:
	       		fi
	        else
	        	hisat2 -q -x ${reference_sp4}/${indx_sp4_directory}/${indx_name_sp4} -S ${head_dir}/${sp4}/hisat_aligned/${align_name}.sam -U ${head_dir}/${sp4}/trimmed/${align_name}.fq 
			fi
		done

	##convert to bam and sort

	cd ${head_dir}/${sp4}/hisat_aligned/

	module load samtools/1.10

	for k in *.sam
		do
			convert2bam_name=$(echo $k | awk -F ".sam" '{print $1}')
			echo ${convert2bam_name}
			samtools view -S -b -h $k > ${convert2bam_name}.bam 
		done
	
	for k in *.bam
		do
			sort_name=$(echo $k | awk -F ".bam" '{print $1}')
			echo ${sort_name}	
			samtools sort ${sort_name}.bam -o ${sort_name}_sorted.bam
		done

	## Stringtie Transcript Creation

	cd ${head_dir}/${sp4}/hisat_aligned/
	mkdir ${head_dir}/${sp4}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for m in *_sorted.bam
		do
	        stringtie_r1_name=$(echo $m | awk -F "_sorted.bam" '{print $1}')
	        echo ${stringtie_r1_name}
	        stringtie ${head_dir}/${sp4}/hisat_aligned/$m -o ${head_dir}/${sp4}/stringtie/${stringtie_r1_name}.gtf -A -C -G ${reference_sp4}/${reference_gff3_sp4} -l ${sp4}_transcripts
		done

	cd ${head_dir}/${sp4}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	stringtie ${head_dir}/${sp4}/stringtie/* --merge -o merged.gtf -G ${reference_sp4}/${reference_gff3_sp4} -l merged_transcripts

	mkdir ${head_dir}/${sp4}/stringtie/stringtie_final_output/
	cd ${head_dir}/${sp4}/hisat_aligned/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for n in *_sorted.bam
		do
 	       stringtie_r3_name=$(echo $n | awk -F "_sorted.bam" '{print $1}')
 	       echo ${stringtie_r3_name}
 	       stringtie ${head_dir}/${sp4}/hisat_aligned/$n -o ${head_dir}/${sp4}/stringtie/stringtie_final_output/${stringtie_r3_name}.gtf -A -C -G ${head_dir}/${sp4}/stringtie/merged.gtf
		done

	## Converting transcripts to fasta

	export PATH=$PATH:/scratch/dtn2an/gffread

	cd ${head_dir}/${sp4}/stringtie/stringtie_final_output/
	mkdir ${head_dir}/${sp4}/fastasfromtranscriptome/
	
	for o in *.gtf
		do
 	 		getfasta_name=$(echo $o | awk -F ".gtf" '{print $1}')
 	        echo ${getfasta_name}
 	        gffread -w ${head_dir}/${sp4}/${getfasta_name}_transcripts.fa -g ${reference_sp4}/${reference_genome_fasta_sp4} ${head_dir}/${sp4}/stringtie/stringtie_final_output/$o
		done
	
	##TransDecoder Step

	mkdir ${head_dir}/${sp4}/TransDecoder_${sp4}_output/
	cd ${head_dir}/${sp4}/fastasfromtranscriptome/

	export PATH=$PATH:/scratch/dtn2an/TransDecoder

	for p in *.fa
		do
       	 	name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		 echo $name
       	 	TransDecoder.LongOrfs -t ${head_dir}/${sp4}/fastasfromtranscriptome/$p -O ${head_dir}/${sp4}/TransDecoder_${sp4}_output/
		done

	cd ${head_dir}/${sp4}/TransDecoder_${sp4}_output/

	for q in *.fa
		do
  			 name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		 echo $name
       		 TransDecoder.Predict -t ${head_dir}/${sp4}/fastasfromtranscriptome/$q -O ${head_dir}/${sp4}/TransDecoder_${sp4}_output/
		done
	
	mkdir ${head_dir}/${sp4}/CD_HIT_output
	
	export PATH=$PATH:/scratch/dtn2an/cd-hit-v4.8.1-2019-0228

	./cd-hit -i ${head_dir}/${sp4}/TransDecoder_${sp4}_output/longest_orfs.pep -o ${head_dir}/${sp4}/CD_HIT_output/cdhit_longest_orfs.pep -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp4}/TransDecoder_${sp4}_output/longest_orfs.gff3 -o ${head_dir}/${sp4}/CD_HIT_output/cdhit_longest_orfs.gff3 -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp4}/TransDecoder_${sp4}_output/longest_orfs.cds -o ${head_dir}/${sp4}/CD_HIT_output/cdhit_longest_orfs.cds -c 0.95 -d 0
		
	echo "species four, three shared done"
else
	echo "species four, three shared skipped"		
fi

##########################################################################################################################################

if [[ ${share_genome} == 3 && (${species_number} == 5) ]]
then 

	###########################################################################
	###################### AUTONOMOUS STEPS ##############################
	###########################################################################

	##Trim_Galore step
	
	##sp5
	mkdir ${head_dir}/${sp5}/trimmed/
	cd /${pre_trim_sp5}/

	module load trimgalore/0.6.4

	for i in *.fastq
		do
			if [[ $i == *"_R1"* || $i == *"_R2"* ]] 	
        	then
        		if [[ $i == *"_R1"* ]]
        		then
        			trim_name=$(echo $i | awk -F _ '{print $1}')
        			echo ${trim_name}
        			trim_galore --fastqc --paired -o ${head_dir}/${sp5}/trimmed ${pre_trim_sp5}/${trim_name}_R1.fastq ${pre_trim_sp5}/${trim_name}_R2.fastq
				else
					:
				fi
			else
				trim_name=$(echo $i | awk -F ".fastq" '{print $1}')
				echo ${trim_name}
				trim_galore --fastqc -o ${head_dir}/${sp5}/trimmed ${pre_trim_sp5}/$i
			fi
		done

	cd ${head_dir}/${sp5}/trimmed

	for j in *.fq
        do
        	detrim_name=$(echo $j | awk -F "_trimmed" '{print $1}')
        	echo ${detrim_name}
        	if [[ $j == *"trimmed"* ]]
            then
            	mv ${detrim_name}_trimmed.fq ${detrim_name}.fq
            else
            	:
            fi
    	 done


	for j in *.fq
        do
        	change_name=$(echo $j | awk -F _ '{print $1}')
            echo ${change_name}
            if [[ $j == *"_val_1"* || $j == *"_val_2"* ]]
            then
                if [[ $j == *"_val_1"* ]]
                then
                    mv ${change_name}_R1_val_1.fq ${change_name}_R1.fq
                else
                    mv ${change_name}_R2_val_2.fq ${change_name}_R2.fq
                fi
            else
                :
            fi
        done

	##Hisat2 alignment step
	mkdir ${head_dir}/${sp5}/hisat_aligned/
	cd ${head_dir}/${sp5}/trimmed/

	module load gcc/7.1.0
	module load hisat2/2.1.0

	for j in *.fq
		do
	  		align_name=$(echo $j | awk -F _ '{print $1}')
	        echo ${align_name}
	       
	        if [[ $j == *"_R1"* || $j == *"_R2"* ]]
	        then 	
	    		if [[ "$j" == *"_R1"* ]]
	    		then
	        		hisat2 -q -x ${reference_sp5}/${indx_sp5_directory}/${indx_name_sp5} -S ${head_dir}/${sp5}/hisat_aligned/${align_name}.sam -1 ${head_dir}/${sp5}/trimmed/${align_name}_R1.fq -2 ${head_dir}/${sp5}/trimmed/${align_name}_R2.fq
				else
					:
				fi	        
	        else
	        	hisat2 -q -x ${reference_sp5}/${indx_sp5_directory}/${indx_name_sp5} -S ${head_dir}/${sp5}/hisat_aligned/${align_name}.sam -U ${head_dir}/${sp5}/trimmed/${align_name}.fq 
			fi
		done

	
	##convert to bam and sort

	cd ${head_dir}/${sp5}/hisat_aligned/

	module load samtools/1.10

	for k in *.sam
		do
			convert2bam_name=$(echo $k | awk -F ".sam" '{print $1}')
			echo ${convert2bam_name}
			samtools view -S -b -h $k > ${convert2bam_name}.bam 
		done
	
	for k in *.bam
		do
			sort_name=$(echo $k | awk -F ".bam" '{print $1}')
			echo ${sort_name}	
			samtools sort ${sort_name}.bam -o ${sort_name}_sorted.bam
		done

	## Stringtie Transcript Creation

	cd ${head_dir}/${sp5}/hisat_aligned/
	mkdir ${head_dir}/${sp5}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for m in *_sorted.bam
		do
	        stringtie_r1_name=$(echo $m | awk -F "_sorted.bam" '{print $1}')
	        echo ${stringtie_r1_name}
	        stringtie ${head_dir}/${sp5}/hisat_aligned/$m -o ${head_dir}/${sp5}/stringtie/${stringtie_r1_name}.gtf -A -C -G ${reference_sp5}/${reference_gff3_sp5} -l ${sp5}_transcripts
		done

	cd ${head_dir}/${sp5}/stringtie/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	stringtie ${head_dir}/${sp5}/stringtie/* --merge -o merged.gtf -G ${reference_sp5}/${reference_gff3_sp5} -l merged_transcripts

	mkdir ${head_dir}/${sp5}/stringtie/stringtie_final_output/
	cd ${head_dir}/${sp5}/hisat_aligned/

	module load gcc/7.1.0
	module load stringtie/2.0.6

	for n in *_sorted.bam
		do
		   stringtie_r3_name=$(echo $n | awk -F "_sorted.bam" '{print $1}')
 	       echo ${stringtie_r3_name}
 	       stringtie ${head_dir}/${sp5}/hisat_aligned/$n -o ${head_dir}/${sp5}/stringtie/stringtie_final_output/${stringtie_r3_name}.gtf -A -C -G ${head_dir}/${sp5}/stringtie/merged.gtf
		done

	## Converting transcripts to fasta

	export PATH=$PATH:/scratch/dtn2an/gffread

	cd ${head_dir}/${sp5}/stringtie/stringtie_final_output/
	mkdir ${head_dir}/${sp5}/fastasfromtranscriptome/
	
	for o in *.gtf
		do
 	 		getfasta_name=$(echo $o | awk -F ".gtf" '{print $1}')
 	        echo ${getfasta_name}
 	        gffread -w ${head_dir}/${sp5}/fastasfromtranscriptome/${getfasta_name}_transcripts.fa -g ${reference_sp5}/${reference_genome_fasta_sp5} ${head_dir}/${sp5}/stringtie/stringtie_final_output/$o
		done
		
	##TransDecoder Step

	mkdir ${head_dir}/${sp5}/TransDecoder_${sp5}_output/
	cd ${head_dir}/${sp5}/fastasfromtranscriptome/

	export PATH=$PATH:/scratch/dtn2an/TransDecoder

	for p in *.fa
		do
       	 	name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       		 echo $name
       	 	TransDecoder.LongOrfs -t ${head_dir}/${sp5}/fastasfromtranscriptome/$p -O ${head_dir}/${sp5}/TransDecoder_${sp5}_output/
		done

	cd ${head_dir}/${sp5}/TransDecoder_${sp5}_output/

	for q in *.fa
		do
  		 name=$(echo $i | awk -F "_transcripts.fa" '{print $1}')
       	 echo $name
       	 TransDecoder.Predict -t ${head_dir}/${sp5}/fastasfromtranscriptome/$q -O ${head_dir}/${sp5}/TransDecoder_${sp5}_output/
		done
		
	mkdir ${head_dir}/${sp5}/CD_HIT_output
	
	export PATH=$PATH:/scratch/dtn2an/cd-hit-v4.8.1-2019-0228

	./cd-hit -i ${head_dir}/${sp5}/TransDecoder_${sp5}_output/longest_orfs.pep -o ${head_dir}/${sp5}/CD_HIT_output/cdhit_longest_orfs.pep -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp5}/TransDecoder_${sp5}_output/longest_orfs.gff3 -o ${head_dir}/${sp5}/CD_HIT_output/cdhit_longest_orfs.gff3 -c 0.95 -d 0
	./cd-hit -i ${head_dir}/${sp5}/TransDecoder_${sp5}_output/longest_orfs.cds -o ${head_dir}/${sp5}/CD_HIT_output/cdhit_longest_orfs.cds -c 0.95 -d 0
		
	echo "species five, three shared done"
else
	echo "species five, three shared skipped"
fi

##########################################################################################################################################
##########################################################################################################################################

