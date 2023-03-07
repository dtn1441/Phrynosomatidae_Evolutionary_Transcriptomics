#!/bin/bash
#
#
#SBATCH -J stringtie
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 0-20:00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/stringtie.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/stringtie.err
#SBATCH --account=coxlab
#SBATCH --partition=standard

cd /scratch/dtn2an/workspace/refseq_aligned/4646_merged/

module load gcc/7.1.0
module load stringtie/2.0.6
module load cufflinks/2.2.1

refgen_gff=/scratch/dtn2an/genomes/refseq_genome/ncbi-genomes-2021-09-15/GCF_019175285.1_SceUnd_v1.1_genomic.gff
refgen_fasta=/scratch/dtn2an/genomes/refseq_genome/ncbi-genomes-2021-09-15/GCF_019175285.1_SceUnd_v1.1_genomic.fna


##first round creates gtf files from an individual's transcripts given the reference genome

for i in *.bam
do
	stringtie_r1_name=$(echo $i | awk -F "_lanes" '{print $1}')
	species_name=$(echo $i | awk -F "[0-9]" '{print $1}')
	echo ${stringtie_r1_name}
	stringtie $i -o ${stringtie_r1_name}.gtf -A -C -G ${refgen_gff} -l ${species_name}_transcripts
done

mkdir cote
mkdir sccl
mkdir scja
mkdir scma
mkdir scme
mkdir scsl
mkdir scun
mkdir scvi
mkdir uror
mkdir utst

for i in *.gtf
do
	name=$(echo $i | awk -F . '{print $1}')
	if [[ ${name} == COTE* ]]
	then
		mv $i ./cote/
	elif [[ ${name} == SCCL* ]]
	then
		mv $i ./sccl/
	elif [[ ${name} == SCJA* ]]
    then
            	mv $i ./scja/
	elif [[ ${name} == SCMA* ]]
	then
		mv $i ./scma/
	elif [[ ${name} == SCME* ]]
	then
		mv $i ./scme/
	elif [[ ${name} == SCSL* ]]
	then
		mv $i ./scsl/
	elif [[ ${name} == SCUN* ]]
	then
            	mv $i ./scun/
	elif [[ ${name} == SCVI* ]]
	then
		mv $i ./scvi/
	elif [[ ${name} == UROR* ]]
	then
		mv $i ./uror/
	elif [[ ${name} == UTST* ]]
	then
		mv $i ./utst/
	else
		:
	fi
done



## Second part merges transcripts from the same species to create a full transcriptome


cd ./cote/

stringtie ./* --merge -o cote_merged.gtf -G ${refgen_gff} -l COTE_transcripts

cd ../sccl/

stringtie ./* --merge -o sccl_merged.gtf -G ${refgen_gff} -l SCCL_transcripts

cd ../scja/

stringtie ./* --merge -o scja_merged.gtf -G ${refgen_gff} -l SCJA_transcripts

cd ../scma/

stringtie ./* --merge -o scma_merged.gtf -G ${refgen_gff} -l SCMA_transcripts

cd ../scme/

stringtie ./* --merge -o scme_merged.gtf -G ${refgen_gff} -l SCME_transcripts

cd ../scsl/

stringtie ./* --merge -o scsl_merged.gtf -G ${refgen_gff} -l SCSL_transcripts

cd ../scun/

stringtie ./* --merge -o scun_merged.gtf -G ${refgen_gff} -l SCUN_transcripts

cd ../scvi/

stringtie ./* --merge -o scvi_merged.gtf -G ${refgen_gff} -l SCVI_transcripts

cd ../uror/

stringtie ./* --merge -o uror_merged.gtf -G ${refgen_gff} -l UROR_transcripts

cd ../utst/

stringtie ./* --merge -o utst_merged.gtf -G ${refgen_gff} -l UTST_transcripts

cd ..


## Third part aligns the reads to the new transcritpome we've created


mkdir stringtie_finished


for i in *_merged.bam
do
	stringtie_r3_name=$(echo $i | awk -F "_merged" '{print $1}')
	species_name=$(echo $i | awk -F "[0-9]" '{print $1}')
 	echo ${stringtie_r3_name}
	if [[ ${species_name} == COTE ]]
 	then
		stringtie $i -o ./stringtie_finished/${stringtie_r3_name}.gtf -A -C -G ./cote/cote_merged.gtf
	elif [[ ${species_name} == SCCL ]]
	then
		stringtie $i -o ./stringtie_finished/${stringtie_r3_name}.gtf -A -C -G ./sccl/sccl_merged.gtf
	elif [[ ${species_name} == SCJA ]]
	then
               stringtie $i -o ./stringtie_finished/${stringtie_r3_name}.gtf -A -C -G ./scja/scja_merged.gtf
	elif [[ ${species_name} == SCMA ]]
	then
		stringtie $i -o ./stringtie_finished/${stringtie_r3_name}.gtf -A -C -G ./scma/scma_merged.gtf
	elif [[ ${species_name} == SCME ]]
	then
		stringtie $i -o ./stringtie_finished/${stringtie_r3_name}.gtf -A -C -G ./scme/scme_merged.gtf
	elif [[ ${species_name} == SCSL ]]
	then
		stringtie $i -o ./stringtie_finished/${stringtie_r3_name}.gtf -A -C -G ./scsl/scsl_merged.gtf
	elif [[ ${species_name} == SCUN ]]
    then    
	    	stringtie $i -o ./stringtie_finished/${stringtie_r3_name}.gtf -A -C -G ./scun/scun_merged.gtf
	elif [[ ${species_name} == SCVI ]]
    then        
		stringtie $i -o ./stringtie_finished/${stringtie_r3_name}.gtf -A -C -G ./scvi/scvi_merged.gtf
	elif [[ ${species_name} == UROR ]]
	then
		stringtie $i -o ./stringtie_finished/${stringtie_r3_name}.gtf -A -C -G ./uror/uror_merged.gtf
	elif [[ ${species_name} == UTST ]]
	then
		stringtie $i -o ./stringtie_finished/${stringtie_r3_name}.gtf -A -C -G ./utst/utst_merged.gtf
	else
		:
	fi
done




##Convert to fasta

cd stringtie_finished/

for i in *.gtf
do
	getfasta_name=$(echo $i | awk -F ".gtf" '{print $1}')
 	echo ${getfasta_name}
 	gffread -w ${getfasta_name}_transcripts.fa -g ${refgen_fasta} $i
done
