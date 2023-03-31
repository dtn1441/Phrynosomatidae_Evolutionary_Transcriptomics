## first I had to make a conda environment for pixy by using
## "conda create --name pixy"

conda activate pixy

conda install --yes -c conda-forge pixy
conda install --yes -c bioconda htslib


cd /scratch/dtn2an/workspace/refseq_aligned/scun_cat/variant_calls/post_processing/main_chroms

pop="scun_pop_ordered.txt"
scun_bed="/scratch/dtn2an/"

##compress and index our vcfs
for i in *_scun_cat_snps_filtered.vcf
do
	bgzip $i
	tabix ${i}.gz
done


##calculate pi from each chromosome

for i in *.gz
do
	name=$(echo $i | awk -F ".vcf" '{print $1}')
	chrom=$(echo $i | awk -F "_scun" '{print $1}')
	pixy --stats pi \
	--vcf $i \
	--populations ${pop} \
	--window_size 10000 \
	--chromosomes chrom
done
##This doesn't work because pixy doesn't let you specify output, so it's all just replacing itself .... lame




pixy --stats pi \
--vcf NC_056522.1_scun_cat_snps_filtered.vcf.gz \
--populations ${pop} \
--window_size 10000 
