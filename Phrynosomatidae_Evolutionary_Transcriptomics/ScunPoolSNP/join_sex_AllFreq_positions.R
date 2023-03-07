library("plyr")
library("tidyr")
library("dplyr")
library("ggplot2")
library("magrittr")
library("tidyverse")

setwd("/scratch/dtn2an/workspace/refseq_aligned/merge_runs/tissue_merged/mpileups/finished/PoolSNP/x_chromosomes/allele_freqs/")


AF10_X_freq <- read.delim("AF10_chrom10_X_PoolSNP.frq")
AF11_X_freq <- read.delim("AF11_chrom10_X_PoolSNP.frq")
AF13_X_freq <- read.delim("AF13_chrom10_X_PoolSNP.frq")
AF4_X_freq <- read.delim("AF4_chrom10_X_PoolSNP.frq")
AF7_X_freq <- read.delim("AF7_chrom10_X_PoolSNP.frq")
AM12_X_freq <- read.delim("AM12_chrom10_X_PoolSNP.frq")
AM14_X_freq <-  read.delim("AM14_chrom10_X_PoolSNP.frq")
AM15_X_freq <- read.delim("AM15_chrom10_X_PoolSNP.frq")
AM5_X_freq <- read.delim("AM5_chrom10_X_PoolSNP.frq")
AM6_X_freq <- read.delim("AM6_chrom10_X_PoolSNP.frq")
AM8_X_freq <- read.delim("AM8_chrom10_X_PoolSNP.frq")
AM9_X_freq <- read.delim("AM9_chrom10_X_PoolSNP.frq")
HF40_X_freq <- read.delim("HF40_chrom10_X_PoolSNP.frq")
HF41_X_freq <- read.delim("HF41_chrom10_X_PoolSNP.frq")
HF42_X_freq <- read.delim("HF42_chrom10_X_PoolSNP.frq")
HF44_X_freq <- read.delim("HF44_chrom10_X_PoolSNP.frq")
HF49_X_freq <- read.delim("HF49_chrom10_X_PoolSNP.frq")
HM43_X_freq <- read.delim("HM43_chrom10_X_PoolSNP.frq")
HM45_X_freq <- read.delim("HM45_chrom10_X_PoolSNP.frq")
HM46_X_freq <- read.delim("HM46_chrom10_X_PoolSNP.frq")
HM47_X_freq <- read.delim("HM47_chrom10_X_PoolSNP.frq")
HM48_X_freq <- read.delim("HM48_chrom10_X_PoolSNP.frq")
YF26_X_freq <- read.delim("YF26_chrom10_X_PoolSNP.frq")
YF27_X_freq <- read.delim("YF27_chrom10_X_PoolSNP.frq")
YF29_X_freq <- read.delim("YF29_chrom10_X_PoolSNP.frq")
YF31_X_freq <- read.delim("YF31_chrom10_X_PoolSNP.frq")
YF32_X_freq <- read.delim("YF32_chrom10_X_PoolSNP.frq")
YF33_X_freq <- read.delim("YF33_chrom10_X_PoolSNP.frq")
YM21_X_freq <- read.delim("YM21_chrom10_X_PoolSNP.frq")
YM22_X_freq <- read.delim("YM22_chrom10_X_PoolSNP.frq")
YM23_X_freq <- read.delim("YM23_chrom10_X_PoolSNP.frq")
YM24_X_freq <- read.delim("YM24_chrom10_X_PoolSNP.frq")
YM25_X_freq <- read.delim("YM25_chrom10_X_PoolSNP.frq")
YM28_X_freq <- read.delim("YM28_chrom10_X_PoolSNP.frq")



AF10_X_freq %>%
  select(POS, N_ALLELES) -> AF10_X_locfreq
AF11_X_freq %>%
  select(POS, N_ALLELES) -> AF11_X_locfreq
AF13_X_freq %>%
  select(POS, N_ALLELES) -> AF13_X_locfreq
AF4_X_freq %>%
  select(POS, N_ALLELES) -> AF4_X_locfreq
AF7_X_freq %>%
  select(POS, N_ALLELES) -> AF7_X_locfreq
AM12_X_freq %>%
  select(POS, N_ALLELES) -> AM12_X_locfreq
AM14_X_freq %>%
  select(POS, N_ALLELES) -> AM14_X_locfreq
AM15_X_freq %>%
  select(POS, N_ALLELES) -> AM15_X_locfreq
AM5_X_freq %>%
  select(POS, N_ALLELES) -> AM5_X_locfreq
AM6_X_freq %>%
  select(POS, N_ALLELES) -> AM6_X_locfreq
AM8_X_freq %>%
  select(POS, N_ALLELES) -> AM8_X_locfreq
AM9_X_freq %>%
  select(POS, N_ALLELES) -> AM9_X_locfreq
HF40_X_freq %>%
  select(POS, N_ALLELES) -> HF40_X_locfreq
HF41_X_freq %>%
  select(POS, N_ALLELES) -> HF41_X_locfreq
HF42_X_freq %>%
  select(POS, N_ALLELES) -> HF42_X_locfreq
HF44_X_freq %>%
  select(POS, N_ALLELES) -> HF44_X_locfreq
HF49_X_freq %>%
  select(POS, N_ALLELES) -> HF49_X_locfreq
HM43_X_freq %>%
  select(POS, N_ALLELES) -> HM43_X_locfreq
HM45_X_freq %>%
  select(POS, N_ALLELES) -> HM45_X_locfreq
HM46_X_freq %>%
  select(POS, N_ALLELES) -> HM46_X_locfreq
HM47_X_freq %>%
  select(POS, N_ALLELES) -> HM47_X_locfreq
HM48_X_freq %>%
  select(POS, N_ALLELES) -> HM48_X_locfreq
YF26_X_freq %>%
  select(POS, N_ALLELES) -> YF26_X_locfreq
YF27_X_freq %>%
  select(POS, N_ALLELES) -> YF27_X_locfreq
YF29_X_freq %>%
  select(POS, N_ALLELES) -> YF29_X_locfreq
YF31_X_freq %>%
  select(POS, N_ALLELES) -> YF31_X_locfreq
YF32_X_freq %>%
  select(POS, N_ALLELES) -> YF32_X_locfreq
YF33_X_freq %>%
  select(POS, N_ALLELES) -> YF33_X_locfreq
YM21_X_freq %>%
  select(POS, N_ALLELES) -> YM21_X_locfreq
YM22_X_freq %>%
  select(POS, N_ALLELES) -> YM22_X_locfreq
YM23_X_freq %>%
  select(POS, N_ALLELES) -> YM23_X_locfreq
YM24_X_freq %>%
  select(POS, N_ALLELES) -> YM24_X_locfreq
YM25_X_freq %>%
  select(POS, N_ALLELES) -> YM25_X_locfreq
YM28_X_freq %>%
  select(POS, N_ALLELES) -> YM28_X_locfreq


female_list <- list(AF10_X_locfreq, AF11_X_locfreq, AF13_X_locfreq, AF4_X_locfreq, AF7_X_locfreq, HF40_X_locfreq, HF41_X_locfreq, HF42_X_locfreq, HF44_X_locfreq, HF49_X_locfreq, YF26_X_locfreq, YF27_X_locfreq, YF29_X_locfreq, YF31_X_locfreq, YF32_X_locfreq, YF33_X_locfreq)
female_list %>% reduce(full_join, by='POS') -> female_X_allele_freqs
write.csv(female_X_allele_freqs, "./female_X_allele_freqs.csv")

male_list <- list(AM12_X_locfreq, AM14_X_locfreq, AM15_X_locfreq, AM5_X_locfreq, AM6_X_locfreq, AM8_X_locfreq, AM9_X_locfreq, HM43_X_locfreq, HM45_X_locfreq, HM46_X_locfreq, HM47_X_locfreq, HM48_X_locfreq, YM21_X_locfreq, YM22_X_locfreq, YM23_X_locfreq, YM24_X_locfreq, YM25_X_locfreq, YM28_X_locfreq)
male_list %>% reduce(full_join, by='POS') -> male_X_allele_freqs
write.csv(male_X_allele_freqs, "./male_X_allele_freqs.csv")

