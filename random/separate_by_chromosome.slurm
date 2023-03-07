#!/bin/bash
#
#
#SBATCH -J separate_by_chromosome
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 0-40:00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/separate_by_chromsome.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/separate_by_chromosome.err
#SBATCH --account=coxlab
#SBATCH --partition=standard

cd /scratch/dtn2an/workspace/refseq_aligned/merge_runs/tissue_merged/mpileups/


for i in *.mpileup
do
        name=$(echo $i | awk -F	. '{print $1}')
        echo $name
        cat $i | grep "NC_056525.1" > ./finished/${name}_chrom4.mpileup
done

for i in *.mpileup
do
        name=$(echo $i | awk -F	. '{print $1}')
        echo $name
        cat $i | grep "NC_056526.1" > ./finished/${name}_chrom5.mpileup
done

for i in *.mpileup
do
        name=$(echo $i | awk -F	. '{print $1}')
        echo $name
        cat $i | grep "NC_056527.1" > ./finished/${name}_chrom6.mpileup
done

for i in *.mpileup
do
        name=$(echo $i | awk -F	. '{print $1}')
        echo $name
        cat $i | grep "NC_056528.1" > ./finished/${name}_chrom7.mpileup
done

for i in *.mpileup
do
        name=$(echo $i | awk -F	. '{print $1}')
        echo $name
        cat $i | grep "NC_056529.1" > ./finished/${name}_chrom8.mpileup
done

for i in *.mpileup
do
        name=$(echo $i | awk -F	. '{print $1}')
        echo $name
        cat $i | grep "NC_056530.1" > ./finished/${name}_chrom9.mpileup
done

for i in *.mpileup
do
        name=$(echo $i | awk -F	. '{print $1}')
        echo $name
        cat $i | grep "NC_056531.1" > ./finished/${name}_chrom10_X.mpileup
done

for i in *.mpileup
do
        name=$(echo $i | awk -F	. '{print $1}')
        echo $name
        cat $i | grep "NC_056532.1" > ./finished/${name}_chrom11.mpileup
done

