if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("SNPRelate")
BiocManager::install("gdsfmt")

install.packages("vcfR", repos = "http://cran.us.r-project.org")

library("SNPRelate")
library("gdsfmt")
library("vcfR")


setwd("/scratch/dtn2an/workspace/refseq_aligned/merge_runs/tissue_merged/temp/")

AF10_vcf <- "AF10_PoolSNP.vcf.vcf"
snpgdsVCF2GDS_R(AF10_vcf, "AF10_PoolSNP.gds", method="copy.num.of.ref")

fout <- filte("test.txt", "wt")
apply.gdsn(AF10_PoolSNP.gds, 1, FUN=cat, as.is="none", file=fout)
close(fout)

readLines("test.txt")
