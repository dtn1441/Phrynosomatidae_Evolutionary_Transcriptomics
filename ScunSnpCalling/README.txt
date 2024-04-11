#SNP Calling from Reference-Mapped Bams

Aligned reads must go through pre-processing, snp calling, trimming to high confidence short list, recalibratation, and finally a second round of snp calling post recalibration

## Pre-processing of Reference-Mapped Bams




bcftools view -i "F_MISSING < 0.01 && DP >= 50" all_scaffs.vcf -o all_scaffs_1_percent_missingness_50DP_filter.vcf
