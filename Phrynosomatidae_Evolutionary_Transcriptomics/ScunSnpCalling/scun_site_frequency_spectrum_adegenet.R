library("adegenet")
library("pegas")


setwd("/scratch/dtn2an/workspace/refseq_aligned/scun_cat/variant_calls/post_processing/main_chroms/")

chrom1.vcf <- read.vcf("NC_056522.1_scun_cat_snps_filtered.vcf", which.loci = 1:3e8)
chrom2.vcf <- read.vcf("NC_056523.1_scun_cat_snps_filtered.vcf", which.loci = 1:3e8)
chrom3.vcf <- read.vcf("NC_056524.1_scun_cat_snps_filtered.vcf", which.loci = 1:2e8)
chrom4.vcf <- read.vcf("NC_056525.1_scun_cat_snps_filtered.vcf", which.loci = 1:2e8)
chrom5.vcf <- read.vcf("NC_056526.1_scun_cat_snps_filtered.vcf", which.loci = 1:1.9e8)
chrom6.vcf <- read.vcf("NC_056527.1_scun_cat_snps_filtered.vcf", which.loci = 1:1.7e8)
chrom7.vcf <- read.vcf("NC_056528.1_scun_cat_snps_filtered.vcf", which.loci = 1:5.1e7)
chrom8.vcf <- read.vcf("NC_056529.1_scun_cat_snps_filtered.vcf", which.loci = 1:3.9e7)
chrom9.vcf <- read.vcf("NC_056530.1_scun_cat_snps_filtered.vcf", which.loci = 1:3.7e7)
chromX.vcf <- read.vcf("NC_056531.1_scun_cat_snps_filtered.vcf", which.loci = 1:1.6e7)
chrom11.vcf <- read.vcf("NC_056532.1_scun_cat_snps_filtered.vcf", which.loci = 1:1.2e7)







chrom1_sfs <- site.spectrum(chrom1.vcf, folded=TRUE, ancestral=NULL, col="blue", main="Chrom1_SFS")
chrom2_sfs <- site.spectrum(chrom2.vcf, folded=TRUE, ancestral=NULL, col="blue", main="Chrom2_SFS")
chrom3_sfs <- site.spectrum(chrom3.vcf, folded=TRUE, ancestral=NULL, col="blue", main="Chrom3_SFS")
chrom4_sfs <- site.spectrum(chrom4.vcf, folded=TRUE, ancestral=NULL, col="blue", main="Chrom4_SFS")
chrom5_sfs <- site.spectrum(chrom5.vcf, folded=TRUE, ancestral=NULL, col="blue", main="Chrom5_SFS")
chrom6_sfs <- site.spectrum(chrom6.vcf, folded=TRUE, ancestral=NULL, col="blue", main="Chrom6_SFS")
chrom7_sfs <- site.spectrum(chrom7.vcf, folded=TRUE, ancestral=NULL, col="blue", main="Chrom7_SFS")
chrom8_sfs <- site.spectrum(chrom8.vcf, folded=TRUE, ancestral=NULL, col="blue", main="Chrom8_SFS")
chrom9_sfs <- site.spectrum(chrom9.vcf, folded=TRUE, ancestral=NULL, col="blue", main="Chrom9_SFS")
chromX_sfs <- site.spectrum(chromX.vcf, folded=TRUE, ancestral=NULL, col="blue", main="ChromX_SFS")
chrom11_sfs <- site.spectrum(chrom11.vcf, folded=TRUE, ancestral=NULL, col="blue", main="Chrom11_SFS")



readRDS(chrom1_sfs, file="/scratch/dtn2an/workspace/refseq_aligned/scun_cat/post_processing/main_chroms/chrom1_sfs.txt")
readRDS(chrom2_sfs, file="/scratch/dtn2an/workspace/refseq_aligned/scun_cat/post_processing/main_chroms/chrom2_sfs.txt")
readRDS(chrom3_sfs, file="/scratch/dtn2an/workspace/refseq_aligned/scun_cat/post_processing/main_chroms/chrom3_sfs.txt")
readRDS(chrom4_sfs, file="/scratch/dtn2an/workspace/refseq_aligned/scun_cat/post_processing/main_chroms/chrom4_sfs.txt")
readRDS(chrom5_sfs, file="/scratch/dtn2an/workspace/refseq_aligned/scun_cat/post_processing/main_chroms/chrom5_sfs.txt")
readRDS(chrom6_sfs, file="/scratch/dtn2an/workspace/refseq_aligned/scun_cat/post_processing/main_chroms/chrom6_sfs.txt")
readRDS(chrom7_sfs, file="/scratch/dtn2an/workspace/refseq_aligned/scun_cat/post_processing/main_chroms/chrom7_sfs.txt")
readRDS(chrom8_sfs, file="/scratch/dtn2an/workspace/refseq_aligned/scun_cat/post_processing/main_chroms/chrom8_sfs.txt")
readRDS(chrom9_sfs, file="/scratch/dtn2an/workspace/refseq_aligned/scun_cat/post_processing/main_chroms/chrom9_sfs.txt")
readRDS(chromX_sfs, file="/scratch/dtn2an/workspace/refseq_aligned/scun_cat/post_processing/main_chroms/chromX_sfs.txt")
readRDS(chrom11_sfs, file="/scratch/dtn2an/workspace/refseq_aligned/scun_cat/post_processing/main_chroms/chrom11_sfs.txt")

