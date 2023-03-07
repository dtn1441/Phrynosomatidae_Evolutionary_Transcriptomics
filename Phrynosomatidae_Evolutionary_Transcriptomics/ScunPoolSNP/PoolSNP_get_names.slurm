#!/bin/bash
#
#
#SBATCH -J PoolSNP_get_names
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 0-40:00:00
#SBATCH --mem 50G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/PoolSNP_get_names.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/PoolSNP_get_nameserr
#SBATCH --account=coxlab
#SBATCH --partition=standard

cd /scratch/dtn2an/workspace/refseq_aligned/merge_runs/tissue_merged/mpileups/finished/PoolSNP/x_chromosomes/allele_freqs/

for i in *.frq
do
	name=$(echo $i | awk -F _ '{print $1}')
	gawk ' 
		NR==FNR { f1[$1][$2] = 1; next }
    		$1 in f1 {
        		for (val in f1[$1]) {
            			if (0+$2 <= 0+val && 0+val <= 0+$3) {
                			print val, $4
            							}
       							 }
    						}
	' $i X_chrom_gene_names_coordinates.bed > ${name}_gene_names.tmp
	echo names for $i are done
done

