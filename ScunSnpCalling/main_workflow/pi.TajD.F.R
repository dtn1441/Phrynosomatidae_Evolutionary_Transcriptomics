library("edgeR")
library("dplyr")
library("limma")
library("car")
library("stringi")
library("magrittr")

as.numeracter <- function(x){
  return(as.numeric(as.character(x)))
}
`%nin%` <- Negate(`%in%`)

###################
##### EdgeR #######
###################

liver.brain.muscle.undulatus.all.ages <- read.csv("/scratch/dtn2an/workspace/refseq_aligned/3906_P1P2_Consolidated_ReadCounts_ScunOnly_LiMuBr.csv", header=TRUE, row.names=1, stringsAsFactors = FALSE)

##View(liver.brain.muscle.undulatus.all.ages)


## remove Scun50(AFB - M), Scun130(SMB - R), Scun64(AMM - T), Scun104(SMM - X), Scun144(SFM - W), Scun153(JFM - U), Scun159(JFM - U), Scun48(AML - B)

## Remove gene "NEB" as it got all screwed up by Excel
liver.brain.muscle.undulatus.all.ages <- liver.brain.muscle.undulatus.all.ages[rownames(liver.brain.muscle.undulatus.all.ages) != "NEB",]

groups_cpm_scun <- c("A",          "A",         "A",          "A",            "A",          "B",                  "B",         "B",        "B",        "B",        "B",         "B",          "C",            "C",          "C",            "C",            "C",               "C",              "D",           "D",            "D",              "D",                "D",             "D",            "E",            "E",           "E",               "E",              "E",             "E",               "F",              "F",             "F",              "F",               "F",               "F",                "M",              "M",             "M",             "M",              "M",            "N",           "N",            "N",              "N",              "N",          "N",             "N",              "O",                   "O",            "O",           "O",              "O",              "O",                  "P",              "P",          "P",             "P",           "P",           "P",           "Q",             "Q",            "Q",              "Q",              "Q",                 "Q",                 "R",                 "R",              "R",                 "R",             "R",                "S",            "S",               "S",                "S",                "S",           "S",                 "T",               "T",             "T",               "T",              "T",                   "T",                  "U",                  "U",                     "U",                  "U",                    "V",                 "V",                 "V",                  "V",                 "V",                  "V",               "W",                   "W",              "W",                 "W",               "W",                   "X",            "X",              "X",             "X",          "X")

sex_groups_scun <- c("Female",    "Female",    "Female",    "Female",       "Female",      "Male",              "Male",      "Male",     "Male",     "Male",     "Male",      "Male",       "Female",       "Female",     "Female",       "Female",       "Female",          "Female",         "Male",        "Male",         "Male",            "Male",            "Male",           "Male",        "Female",        "Female",      "Female",         "Female",        "Female",          "Female",         "Male",         "Male",            "Male",            "Male",           "Male",             "Male",           "Female",         "Female",         "Female",        "Female",        "Female",        "Male",       "Male",          "Male",           "Male",            "Male",       "Male",         "Male",          "Female",              "Female",       "Female",        "Female",        "Female",         "Female",             "Male",          "Male",        "Male",          "Male",       "Male",         "Male",       "Female",         "Female",      "Female",          "Female",         "Female",           "Female",              "Male",             "Male",           "Male",               "Male",            "Male",          "Female",       "Female",          "Female",            "Female",           "Female",       "Female",           "Male",             "Male",        "Male",            "Male",            "Male",                "Male",             "Female",             "Female",                 "Female",             "Female",              "Male",               "Male",              "Male",              "Male",               "Male",               "Male",            "Female",             "Female",         "Female",            "Female",           "Female",             "Male",          "Male",            "Male",        "Male",         "Male")

age_groups_scun <-  c("Adult",    "Adult",      "Adult",    "Adult",        "Adult",       "Adult",             "Adult",     "Adult",    "Adult",    "Adult",    "Adult",     "Adult",      "Neonate",      "Neonate",     "Neonate",     "Neonate",       "Neonate",        "Neonate",         "Neonate",    "Neonate",     "Neonate",         "Neonate",          "Neonate",        "Neonate",     "Maturing",     "Maturing",      "Maturing",     "Maturing",       "Maturing",        "Maturing",      "Maturing",     "Maturing",        "Maturing",         "Maturing",       "Maturing",       "Maturing",          "Adult",          "Adult",          "Adult",         "Adult",         "Adult",         "Adult",      "Adult",         "Adult",          "Adult",           "Adult",       "Adult",        "Adult",       "Neonate",            "Neonate",       "Neonate",        "Neonate",       "Neonate",       "Neonate",            "Neonate",      "Neonate",      "Neonate",      "Neonate",     "Neonate",      "Neonate",    "Maturing",       "Maturing",    "Maturing",        "Maturing",       "Maturing",          "Maturing",          "Maturing",         "Maturing",        "Maturing",           "Maturing",        "Maturing",       "Adult",       "Adult",           "Adult",             "Adult",            "Adult",         "Adult",           "Adult",           "Adult",        "Adult",          "Adult",           "Adult",               "Adult",            "Neonate",             "Neonate",                "Neonate",          "Neonate",             "Neonate",             "Neonate",            "Neonate",          "Neonate",             "Neonate",            "Neonate",         "Maturing",             "Maturing",       "Maturing",         "Maturing",         "Maturing",        "Maturing",         "Maturing",        "Maturing",    "Maturing",     "Maturing")

tissue_groups_scun <- c("Liver",  "Liver",     "Liver",     "Liver",        "Liver",       "Liver",             "Liver",     "Liver",    "Liver",    "Liver",    "Liver",     "Liver",       "Liver",        "Liver",       "Liver",       "Liver",       "Liver",             "Liver",         "Liver",       "Liver",      "Liver",            "Liver",           "Liver",          "Liver",       "Liver",        "Liver",         "Liver",        "Liver",          "Liver",            "Liver",         "Liver",          "Liver",          "Liver",           "Liver",          "Liver",          "Liver",             "Brain",        "Brain",           "Brain",         "Brain",         "Brain",         "Brain",       "Brain",        "Brain",           "Brain",           "Brain",       "Brain",        "Brain",         "Brain",             "Brain",        "Brain",           "Brain",         "Brain",         "Brain",              "Brain",         "Brain",       "Brain",        "Brain",      "Brain",        "Brain",      "Brain",          "Brain",       "Brain",           "Brain",          "Brain",             "Brain",             "Brain",            "Brain",          "Brain",              "Brain",          "Brain",           "Muscle",     "Muscle",          "Muscle",            "Muscle",          "Muscle",         "Muscle",          "Muscle",          "Muscle",       "Muscle",        "Muscle",           "Muscle",                "Muscle",            "Muscle",                "Muscle",                "Muscle",           "Muscle",              "Muscle",              "Muscle",             "Muscle",            "Muscle",              "Muscle",            "Muscle",         "Muscle",                 "Muscle",         "Muscle",          "Muscle",           "Muscle",          "Muscle",         "Muscle",           "Muscle",      "Muscle",      "Muscle")

group_info <- as.data.frame(cbind(sex_groups_scun, age_groups_scun, tissue_groups_scun, groups_cpm_scun))



## remove Scun48(AML - A), Scun50(AFB - M), Scun130(SMB - R), Scun64(AMM - T), Scun104(SMM - X), Scun144(SFM - W), Scun153(JFM - U), Scun159(JFM - U)

##Omnibus

scun_omnibus <- liver.brain.muscle.undulatus.all.ages

colnames(scun_omnibus[4])
scun_omnibus_dropped <- scun_omnibus[,-c(4,41,74,85,88,90,104,106)]



scun_df_omnibus <- DGEList(counts=scun_omnibus_dropped, group=groups_cpm_scun)
scun_df_omnibus$samples
nrow(scun_df_omnibus)


## we may want to adjust from defaults
keep.auto_scun_omnibus <- filterByExpr(scun_df_omnibus, group=groups_cpm_scun)
scun_auto_omnibus <- scun_df_omnibus[keep.auto_scun_omnibus,, keep.lib.sizes=FALSE]
nrow(scun_auto_omnibus)




dge_scun_omnibus <- DGEList(scun_auto_omnibus$counts[apply(scun_auto_omnibus, 1, sum) != 0, ],
                            group=groups_cpm_scun)
dge_scun_omnibus$sampleInfo <- scun_auto_omnibus$SampleInfo
head(dge_scun_omnibus$counts)
nrow(dge_scun_omnibus)

dge_scun_omnibus <-calcNormFactors(dge_scun_omnibus)
dge_scun_omnibus$samples

design_auto_scun_omnibus <- model.matrix(~0+groups_cpm_scun, data=dge_scun_omnibus$samples)
design_auto_scun_omnibus

dge_scun_omnibus <-estimateDisp(dge_scun_omnibus, design_auto_scun_omnibus)
dge_scun_omnibus$common.dispersion




scun_groups_liver_only <- factor(c("A", "A", "A", "A", "A", "B","B","B","B","B","B","B","C","C","C","C","C","C","D","D","D","D","D","D","E","E","E","E","E","E","F","F","F","F","F","F")) 

colnames(liver.brain.muscle.undulatus.all.ages[4])

liver_only <- liver.brain.muscle.undulatus.all.ages[1:37]

liver_only %>% 
  dplyr::select(-starts_with("Scun_48")) -> scun_dropped_liver_only


##do scun groups and scun_dropped equal?
NROW(scun_groups_liver_only)
ncol(scun_dropped_liver_only)

scun_df_liver_only <- DGEList(counts=scun_dropped_liver_only, group=scun_groups_liver_only)
scun_df_liver_only$samples
nrow(scun_df_liver_only)


keep.auto_scun_liver_only <- filterByExpr(scun_df_liver_only, group=scun_groups_liver_only)
scun_auto_liver_only <- scun_df_liver_only[keep.auto_scun_liver_only,, keep.lib.sizes=FALSE]
nrow(scun_auto_liver_only)


scun_groups_brain_only <- factor(c("M","M","M","M","M","N","N","N","N","N","N","N","O","O","O","O","O","O","P","P","P","P","P","P","Q","Q","Q","Q","Q","Q", "R","R","R","R","R")) 

colnames(liver.brain.muscle.undulatus.all.ages[74])

brain_only <- liver.brain.muscle.undulatus.all.ages[38:74]

colnames(brain_only[37])

brain_only %>% 
  dplyr::select(-c(starts_with("Scun_50"), starts_with("Scun_130"))) -> scun_dropped_brain_only




##do scun groups and scun_dropped equal?
NROW(scun_groups_brain_only)
ncol(scun_dropped_brain_only)

scun_df_brain_only <- DGEList(counts=scun_dropped_brain_only, group=scun_groups_brain_only)
scun_df_brain_only$samples
nrow(scun_df_brain_only)



keep.auto_scun_brain_only <- filterByExpr(scun_df_brain_only, group=scun_groups_brain_only)
scun_auto_brain_only <- scun_df_brain_only[keep.auto_scun_brain_only,, keep.lib.sizes=FALSE]
nrow(scun_auto_brain_only)



scun_groups_muscle_only <- factor(c("S","S","S","S","S","S","T","T","T","T","T","T", "U","U","U","U","V","V","V","V","V","V","W","W","W","W","W","X","X","X","X","X")) 

colnames(liver.brain.muscle.undulatus.all.ages[85])

muscle_only <- liver.brain.muscle.undulatus.all.ages[75:111]



scun_dropped_muscle_only <- muscle_only[,-c(11,14,16,30,32)]

muscle_only %>% 
  dplyr::select(-c(starts_with("Scun_64"), starts_with("Scun_104"), starts_with("Scun_144"), starts_with("Scun_153"), starts_with("Scun_159"))) -> scun_dropped_muscle_only

##do scun groups and scun_dropped equal?
NROW(scun_groups_muscle_only)
ncol(scun_dropped_muscle_only)

scun_df_muscle_only <- DGEList(counts=scun_dropped_muscle_only, group=scun_groups_muscle_only)
scun_df_muscle_only$samples
nrow(scun_df_muscle_only)


keep.auto_scun_muscle_only <- filterByExpr(scun_df_muscle_only, group=scun_groups_muscle_only)
scun_auto_muscle_only <- scun_df_muscle_only[keep.auto_scun_muscle_only,, keep.lib.sizes=FALSE]
nrow(scun_auto_muscle_only)



dge_scun_liver_only <- DGEList(scun_auto_liver_only$counts[apply(scun_auto_liver_only, 1, sum) != 0, ],
                               group=scun_groups_liver_only)
dge_scun_liver_only$sampleInfo <- scun_auto_liver_only$SampleInfo
head(dge_scun_liver_only$counts)
nrow(dge_scun_liver_only)

dge_scun_liver_only <-calcNormFactors(dge_scun_liver_only)
dge_scun_liver_only$samples

design_auto_scun_liver_only <- model.matrix(~0+scun_groups_liver_only, data=dge_scun_liver_only$samples)
design_auto_scun_liver_only

dge_scun_liver_only <-estimateDisp(dge_scun_liver_only, design_auto_scun_liver_only)
dge_scun_liver_only$common.dispersion



dge_scun_brain_only <- DGEList(scun_auto_brain_only$counts[apply(scun_auto_brain_only, 1, sum) != 0, ],
                               group=scun_groups_brain_only)
dge_scun_brain_only$sampleInfo <- scun_auto_brain_only$SampleInfo
head(dge_scun_brain_only$counts)
nrow(dge_scun_brain_only)

dge_scun_brain_only <-calcNormFactors(dge_scun_brain_only)
dge_scun_brain_only$samples

design_auto_scun_brain_only <- model.matrix(~0+scun_groups_brain_only, data=dge_scun_brain_only$samples)
design_auto_scun_brain_only

dge_scun_brain_only <-estimateDisp(dge_scun_brain_only, design_auto_scun_brain_only)
dge_scun_brain_only$common.dispersion



dge_scun_muscle_only <- DGEList(scun_auto_muscle_only$counts[apply(scun_auto_muscle_only, 1, sum) != 0, ],
                                group=scun_groups_muscle_only)
dge_scun_muscle_only$sampleInfo <- scun_auto_muscle_only$SampleInfo
head(dge_scun_muscle_only$counts)
nrow(dge_scun_muscle_only)

dge_scun_muscle_only <-calcNormFactors(dge_scun_muscle_only)
dge_scun_muscle_only$samples

design_auto_scun_muscle_only <- model.matrix(~0+scun_groups_muscle_only, data=dge_scun_muscle_only$samples)
design_auto_scun_muscle_only

dge_scun_muscle_only <-estimateDisp(dge_scun_muscle_only, design_auto_scun_muscle_only)
dge_scun_muscle_only$common.dispersion



## Gather logCPMs

liver_cpms <- cpm(dge_scun_liver_only, normalized.lib.sizes = TRUE, log=TRUE, prior.count=2)
muscle_cpms <- cpm(dge_scun_muscle_only, normalized.lib.sizes = TRUE, log=TRUE, prior.count=2)
brain_cpms <- cpm(dge_scun_brain_only, normalized.lib.sizes = TRUE, log=TRUE, prior.count=2)

##hist(liver_cpms)
##hist(muscle_cpms)
##hist(brain_cpms)


tliver_cpms <- t(liver_cpms)
tmuscle_cpms <- t(muscle_cpms)
tbrain_cpms <- t(brain_cpms)


liver_groups_cpm_scun <- c("A",          "A",         "A",          "A",            "A",          "B",                  "B",         "B",        "B",        "B",        "B",         "B",          "C",            "C",          "C",            "C",            "C",               "C",              "D",           "D",            "D",              "D",                "D",             "D",            "E",            "E",           "E",               "E",              "E",             "E",               "F",              "F",             "F",              "F",               "F",               "F")

liver_sex_groups_scun <- c("Female",    "Female",    "Female",    "Female",       "Female",      "Male",              "Male",      "Male",     "Male",     "Male",     "Male",      "Male",       "Female",       "Female",     "Female",       "Female",       "Female",          "Female",         "Male",        "Male",         "Male",            "Male",            "Male",           "Male",        "Female",        "Female",      "Female",         "Female",        "Female",          "Female",         "Male",         "Male",            "Male",            "Male",           "Male",             "Male")

liver_age_groups_scun <-  c("Adult",    "Adult",      "Adult",    "Adult",        "Adult",       "Adult",             "Adult",     "Adult",    "Adult",    "Adult",    "Adult",     "Adult",      "Neonate",      "Neonate",     "Neonate",     "Neonate",       "Neonate",        "Neonate",         "Neonate",    "Neonate",     "Neonate",         "Neonate",          "Neonate",        "Neonate",     "Maturing",     "Maturing",      "Maturing",     "Maturing",       "Maturing",        "Maturing",      "Maturing",     "Maturing",        "Maturing",         "Maturing",       "Maturing",       "Maturing")

brain_groups_cpm_scun <- c("M",              "M",             "M",             "M",              "M",            "N",           "N",            "N",              "N",              "N",          "N",             "N",              "O",                   "O",            "O",           "O",              "O",              "O",                  "P",              "P",          "P",             "P",           "P",           "P",           "Q",             "Q",            "Q",              "Q",              "Q",                 "Q",                 "R",                 "R",              "R",                 "R",             "R")

brain_sex_groups_scun <- c("Female",         "Female",         "Female",        "Female",        "Female",        "Male",       "Male",          "Male",           "Male",            "Male",       "Male",         "Male",          "Female",              "Female",       "Female",        "Female",        "Female",         "Female",             "Male",          "Male",        "Male",          "Male",       "Male",         "Male",       "Female",         "Female",      "Female",          "Female",         "Female",           "Female",              "Male",             "Male",           "Male",               "Male",            "Male")

brain_age_groups_scun <-  c("Adult",          "Adult",          "Adult",         "Adult",         "Adult",         "Adult",      "Adult",         "Adult",          "Adult",           "Adult",       "Adult",        "Adult",       "Neonate",            "Neonate",       "Neonate",        "Neonate",       "Neonate",       "Neonate",            "Neonate",      "Neonate",      "Neonate",      "Neonate",     "Neonate",      "Neonate",    "Maturing",       "Maturing",    "Maturing",        "Maturing",       "Maturing",          "Maturing",          "Maturing",         "Maturing",        "Maturing",           "Maturing",        "Maturing")

muscle_groups_cpm_scun <- c("S",            "S",               "S",                "S",                "S",           "S",                 "T",               "T",             "T",               "T",              "T",                   "T",                  "U",                  "U",                     "U",                  "U",                    "V",                 "V",                 "V",                  "V",                 "V",                  "V",               "W",                   "W",              "W",                 "W",               "W",                   "X",            "X",              "X",             "X",          "X")

muscle_sex_groups_scun <- c("Female",       "Female",          "Female",            "Female",           "Female",       "Female",           "Male",             "Male",        "Male",            "Male",            "Male",                "Male",             "Female",             "Female",                 "Female",             "Female",              "Male",               "Male",              "Male",              "Male",               "Male",               "Male",            "Female",             "Female",         "Female",            "Female",           "Female",             "Male",          "Male",            "Male",        "Male",         "Male")

muscle_age_groups_scun <-  c("Adult",       "Adult",           "Adult",             "Adult",            "Adult",         "Adult",           "Adult",           "Adult",        "Adult",          "Adult",           "Adult",               "Adult",            "Neonate",             "Neonate",                "Neonate",          "Neonate",             "Neonate",             "Neonate",            "Neonate",          "Neonate",             "Neonate",            "Neonate",         "Maturing",             "Maturing",       "Maturing",         "Maturing",         "Maturing",        "Maturing",         "Maturing",        "Maturing",    "Maturing",     "Maturing")

liver_cpms_w_groups <- cbind(tliver_cpms, liver_groups_cpm_scun, liver_sex_groups_scun, liver_age_groups_scun)
brain_cpms_w_groups <- cbind(tbrain_cpms, brain_groups_cpm_scun, brain_sex_groups_scun, brain_age_groups_scun)
muscle_cpms_w_groups <- cbind(tmuscle_cpms, muscle_groups_cpm_scun, muscle_sex_groups_scun, muscle_age_groups_scun)


tliver_cpms_w_groups <- t(liver_cpms_w_groups)
tbrain_cpms_w_groups <- t(brain_cpms_w_groups)
tmuscle_cpms_w_groups <- t(muscle_cpms_w_groups)





##Remove X chromosome
chromosome_names <- read.csv("/scratch/dtn2an/workspace/refseq_aligned/SceUnd1.1_gene_Chr_locations.csv", header=TRUE)

accepted_chromosomes_names <- c("NC_056522.1","NC_056523.1","NC_056524.1","NC_056525.1","NC_056526.1","NC_056527.1","NC_056528.1","NC_056529.1","NC_056530.1","NC_056531.1","NC_056532.1")

accepted_chromosomes_only <- filter(chromosome_names, grepl(paste(accepted_chromosomes_names, collapse="|"), chromosome))

accepted_chromosomes_only$chromosome <- stri_replace_all_regex(accepted_chromosomes_only$chromosome, pattern = c("NC_056522.1","NC_056523.1","NC_056524.1","NC_056525.1","NC_056526.1","NC_056527.1","NC_056528.1","NC_056529.1","NC_056530.1","NC_056531.1","NC_056532.1"), replacement = c("Chrom_1", "Chrom_2", "Chrom_3", "Chrom_4", "Chrom_5", "Chrom_6", "Chrom_7", "Chrom_8", "Chrom_9", "Chrom_10", "Chrom_11"), vectorize=FALSE)

tliver_cpms_w_groups <- tibble::rownames_to_column(as.data.frame(tliver_cpms_w_groups), "gene_name")
as.data.frame(tliver_cpms_w_groups) %>% 
  left_join(accepted_chromosomes_only %>%
    dplyr::select(chromosome, start, stop, gene_name), by = "gene_name") -> tliver_cpms_w_groups
tliver_cpms_w_groups %<>%
  mutate(
  chromosome = case_when(is.na(chromosome) == TRUE ~ "Scaffolds",
  is.na(chromosome) == FALSE ~ chromosome)
  )

tliver_cpms_w_groups %>%
  filter(chromosome != "Chrom_10" & chromosome != "Scaffolds") -> tliver_cpms_w_groups.noX
rownames(tliver_cpms_w_groups.noX) <- tliver_cpms_w_groups.noX$gene_name
tliver_cpms_w_groups.noX %<>%
  dplyr::select(-c("gene_name", "chromosome", "start", "stop"))

liver_cpms_w_groups.noX <- t(tliver_cpms_w_groups.noX)

tbrain_cpms_w_groups <- tibble::rownames_to_column(as.data.frame(tbrain_cpms_w_groups), "gene_name")
as.data.frame(tbrain_cpms_w_groups) %>% 
  left_join(accepted_chromosomes_only %>%
              dplyr::select(chromosome, start, stop, gene_name), by = "gene_name") -> tbrain_cpms_w_groups
tbrain_cpms_w_groups %<>%
  mutate(
    chromosome = case_when(is.na(chromosome) == TRUE ~ "Scaffolds",
                           is.na(chromosome) == FALSE ~ chromosome)
  )

tbrain_cpms_w_groups %>%
  filter(chromosome != "Chrom_10" & chromosome != "Scaffolds") -> tbrain_cpms_w_groups.noX
rownames(tbrain_cpms_w_groups.noX) <- tbrain_cpms_w_groups.noX$gene_name
tbrain_cpms_w_groups.noX %<>%
  dplyr::select(-c("gene_name", "chromosome", "start", "stop"))

brain_cpms_w_groups.noX <- t(tbrain_cpms_w_groups.noX)


tmuscle_cpms_w_groups <- tibble::rownames_to_column(as.data.frame(tmuscle_cpms_w_groups), "gene_name")
as.data.frame(tmuscle_cpms_w_groups) %>% 
  left_join(accepted_chromosomes_only %>%
              dplyr::select(chromosome, start, stop, gene_name), by = "gene_name") -> tmuscle_cpms_w_groups
tmuscle_cpms_w_groups %<>%
  mutate(
    chromosome = case_when(is.na(chromosome) == TRUE ~ "Scaffolds",
                           is.na(chromosome) == FALSE ~ chromosome)
  )

tmuscle_cpms_w_groups %>%
  filter(chromosome != "Chrom_10" & chromosome != "Scaffolds") -> tmuscle_cpms_w_groups.noX
rownames(tmuscle_cpms_w_groups.noX) <- tmuscle_cpms_w_groups.noX$gene_name
tmuscle_cpms_w_groups.noX %<>%
  dplyr::select(-c("gene_name", "chromosome", "start", "stop"))

muscle_cpms_w_groups.noX <- t(tmuscle_cpms_w_groups.noX)

NCOL(liver_cpms_w_groups.noX)
NCOL(brain_cpms_w_groups.noX)
NCOL(muscle_cpms_w_groups.noX)

liver_cpms_w_groups.noX <- as.data.frame(liver_cpms_w_groups.noX)
brain_cpms_w_groups.noX <- as.data.frame(brain_cpms_w_groups.noX)
muscle_cpms_w_groups.noX <- as.data.frame(muscle_cpms_w_groups.noX)


liver_cpms_w_groups.noX <- mutate_all(liver_cpms_w_groups.noX, function(x) as.numeric(as.character(x)))
brain_cpms_w_groups.noX <- mutate_all(brain_cpms_w_groups.noX, function(x) as.numeric(as.character(x)))
muscle_cpms_w_groups.noX <- mutate_all(muscle_cpms_w_groups.noX, function(x) as.numeric(as.character(x)))


liver_cpms_w_groups.noX <- cbind(liver_cpms_w_groups.noX, liver_groups_cpm_scun, liver_sex_groups_scun, liver_age_groups_scun)
brain_cpms_w_groups.noX <- cbind(brain_cpms_w_groups.noX, brain_groups_cpm_scun, brain_sex_groups_scun, brain_age_groups_scun)
muscle_cpms_w_groups.noX <- cbind(muscle_cpms_w_groups.noX, muscle_groups_cpm_scun, muscle_sex_groups_scun, muscle_age_groups_scun)



liver_model <- c()
liver_anova_summaries <- c()
for (i in 1:(NCOL(liver_cpms_w_groups.noX)-3)){
  liver_model[[i]] <- lm(liver_cpms_w_groups.noX[,i] ~ liver_sex_groups_scun + liver_age_groups_scun + liver_sex_groups_scun*liver_age_groups_scun, data=liver_cpms_w_groups.noX)
  liver_anova_summaries[[i]] <- Anova(liver_model[[i]], type=3)
}


brain_model <- c()
brain_anova_summaries <- c()
for (i in 1:(NCOL(brain_cpms_w_groups.noX)-3)){
  brain_model[[i]] <- lm(brain_cpms_w_groups.noX[,i] ~  brain_sex_groups_scun * brain_age_groups_scun, data=brain_cpms_w_groups.noX)
  brain_anova_summaries[[i]] <- Anova(brain_model[[i]], type=3)
}

muscle_model <- c()
muscle_anova_summaries <- c()
for (i in 1:(NCOL(muscle_cpms_w_groups.noX)-3)){
  muscle_model[[i]] <- lm(muscle_cpms_w_groups.noX[,i] ~ muscle_sex_groups_scun + muscle_age_groups_scun + muscle_sex_groups_scun*muscle_age_groups_scun, data=muscle_cpms_w_groups.noX)
  muscle_anova_summaries[[i]] <- Anova(muscle_model[[i]], type=3)
}



##Liver F's
## extract sex effect for each gene into big list

liver_sex_fvals <- c()
for (i in 1:(NCOL(liver_cpms_w_groups.noX)-3)){
  liver_sex_fvals[i] <- liver_anova_summaries[[i]]$`F value`[2]
}

liver_gene_names.noX <-  as.data.frame(colnames(liver_cpms_w_groups.noX[,1:(NCOL(liver_cpms_w_groups.noX)-3)]))

liver_sex_fvalues <- as.data.frame(cbind(liver_gene_names.noX, liver_sex_fvals))
colnames(liver_sex_fvalues) <- c("gene_names", "fvals")
##View(liver_sex_fvalues)


## extract age effect for each gene into big list

liver_age_fvals <- c()
for (i in 1:(NCOL(liver_cpms_w_groups.noX)-3)){
  liver_age_fvals[i] <- liver_anova_summaries[[i]]$`F value`[3]
}

liver_age_fvalues <- as.data.frame(cbind(liver_gene_names.noX, liver_age_fvals))
colnames(liver_age_fvalues) <- c("gene_names", "fvals")
##View(liver_age_fvalues)


## extract Sex by Age effect for each gene into big list

liver_SexAge_fvals <- c()
for (i in 1:(NCOL(liver_cpms_w_groups.noX)-3)){
  liver_SexAge_fvals[i] <- liver_anova_summaries[[i]]$`F value`[4]
}

liver_SexAge_fvalues <- as.data.frame(cbind(liver_gene_names.noX, liver_SexAge_fvals))
colnames(liver_SexAge_fvalues) <- c("gene_names", "fvals")
##View(liver_SexAge_fvalues)
##double check to see that corresponding fvalues are with their right gene.




liver.tmp <- left_join(liver_sex_fvalues, liver_age_fvalues, by = "gene_names")
liver_fvalues <- left_join(liver.tmp, liver_SexAge_fvalues, by = "gene_names")

colnames(liver_fvalues) <- c("gene_names", "sex_fvals", "age_fvals", "SexAge_fvals")

##Brain F's
## extract sex effect for each gene into big list

brain_sex_fvals <- c()
for (i in 1:(NCOL(brain_cpms_w_groups.noX)-3)){
  brain_sex_fvals[i] <- brain_anova_summaries[[i]]$`F value`[2]
}

brain_gene_names.noX <-  as.data.frame(colnames(brain_cpms_w_groups.noX[,1:(NCOL(brain_cpms_w_groups.noX)-3)]))

brain_sex_fvalues <- as.data.frame(cbind(brain_gene_names.noX, brain_sex_fvals))
colnames(brain_sex_fvalues) <- c("gene_names", "fvals")
##View(brain_sex_fvalues)


## extract age effect for each gene into big list

brain_age_fvals <- c()
for (i in 1:(NCOL(brain_cpms_w_groups.noX)-3)){
  brain_age_fvals[i] <- brain_anova_summaries[[i]]$`F value`[3]
}

brain_age_fvalues <- as.data.frame(cbind(brain_gene_names.noX, brain_age_fvals))
colnames(brain_age_fvalues) <- c("gene_names", "fvals")
##View(brain_age_fvalues)


## extract Sex by Age effect for each gene into big list

brain_SexAge_fvals <- c()
for (i in 1:(NCOL(brain_cpms_w_groups.noX)-3)){
  brain_SexAge_fvals[i] <- brain_anova_summaries[[i]]$`F value`[4]
}

brain_SexAge_fvalues <- as.data.frame(cbind(brain_gene_names.noX, brain_SexAge_fvals))
colnames(brain_SexAge_fvalues) <- c("gene_names", "fvals")
##View(brain_SexAge_fvalues)
##double check to see that corresponding fvalues are with their right gene.


brain.tmp <- left_join(brain_sex_fvalues, brain_age_fvalues, by = "gene_names")
brain_fvalues <- left_join(brain.tmp, brain_SexAge_fvalues, by = "gene_names")

colnames(brain_fvalues) <- c("gene_names", "sex_fvals", "age_fvals", "SexAge_fvals")



##Muscle F's
## extract sex effect for each gene into big list

muscle_sex_fvals <- c()
for (i in 1:(NCOL(muscle_cpms_w_groups.noX)-3)){
  muscle_sex_fvals[i] <- muscle_anova_summaries[[i]]$`F value`[2]
}

muscle_gene_names.noX <-  as.data.frame(colnames(muscle_cpms_w_groups.noX[,1:(NCOL(muscle_cpms_w_groups.noX)-3)]))

muscle_sex_fvalues <- as.data.frame(cbind(muscle_gene_names.noX, muscle_sex_fvals))
colnames(muscle_sex_fvalues) <- c("gene_names", "fvals")
##View(muscle_sex_fvalues)


## extract age effect for each gene into big list

muscle_age_fvals <- c()
for (i in 1:(NCOL(muscle_cpms_w_groups.noX)-3)){
  muscle_age_fvals[i] <- muscle_anova_summaries[[i]]$`F value`[3]
}

muscle_age_fvalues <- as.data.frame(cbind(muscle_gene_names.noX, muscle_age_fvals))
colnames(muscle_age_fvalues) <- c("gene_names", "fvals")
##View(muscle_age_fvalues)


## extract Sex by Age effect for each gene into big list

muscle_SexAge_fvals <- c()
for (i in 1:(NCOL(muscle_cpms_w_groups.noX)-3)){
  muscle_SexAge_fvals[i] <- muscle_anova_summaries[[i]]$`F value`[4]
}

muscle_SexAge_fvalues <- as.data.frame(cbind(muscle_gene_names.noX, muscle_SexAge_fvals))
colnames(muscle_SexAge_fvalues) <- c("gene_names", "fvals")
##View(muscle_SexAge_fvalues)
##double check to see that corresponding fvalues are with their right gene.

muscle.tmp <- left_join(muscle_sex_fvalues, muscle_age_fvalues, by = "gene_names")
muscle_fvalues <- left_join(muscle.tmp, muscle_SexAge_fvalues, by = "gene_names")

colnames(muscle_fvalues) <- c("gene_names", "sex_fvals", "age_fvals", "SexAge_fvals")




liver_fvalues %<>%
  mutate(tissue = "Liver")
brain_fvalues %<>%
  mutate(tissue = "Brain")
muscle_fvalues %<>%
  mutate(tissue = "Muscle")


all_fvalues <- rbind(liver_fvalues, brain_fvalues, muscle_fvalues)



##maybe add tau in this data frame. Try calculating tau by age group and merge into fvalue dataset and then eventually add popgen stats to this dataset

liver_cpms_w_groups.noX
brain_cpms_w_groups.noX
muscle_cpms_w_groups.noX




##create average logCPMs for each group

liver_cpms_w_groups.noX %>%
  group_by(liver_groups_cpm_scun) %>%
  summarise_at(vars(c(1:(NCOL(liver_cpms_w_groups.noX)-3))),funs(mean)) -> liver_cpm_means
liver_cpm_means %<>%
  mutate(Sex = case_when(liver_groups_cpm_scun %in% c("A","C","E") ~ "Female",
                         liver_groups_cpm_scun %in% c("B","D","F") ~ "Male"),
         Age = case_when(liver_groups_cpm_scun %in% c("A","B") ~ "Adult",
                         liver_groups_cpm_scun %in% c("C","D") ~ "Neonate",
                         liver_groups_cpm_scun %in% c("E","F") ~ "Maturing")
         )



brain_cpms_w_groups.noX %>%
  group_by(brain_groups_cpm_scun) %>%
  summarise_at(vars(c(1:(NCOL(brain_cpms_w_groups.noX)-3))),funs(mean)) -> brain_cpm_means
brain_cpm_means %<>%
  mutate(Sex = case_when(brain_groups_cpm_scun %in% c("M","O","Q") ~ "Female",
                         brain_groups_cpm_scun %in% c("N","P","R") ~ "Male"),
         Age = case_when(brain_groups_cpm_scun %in% c("M","N") ~ "Adult",
                         brain_groups_cpm_scun %in% c("O","P") ~ "Neonate",
                         brain_groups_cpm_scun %in% c("Q","R") ~ "Maturing")
  )



muscle_cpms_w_groups.noX %>%
  group_by(muscle_groups_cpm_scun) %>%
  summarise_at(vars(c(1:(NCOL(muscle_cpms_w_groups.noX)-3))),funs(mean)) -> muscle_cpm_means
muscle_cpm_means %<>%
  mutate(Sex = case_when(muscle_groups_cpm_scun %in% c("S","U","W") ~ "Female",
                         muscle_groups_cpm_scun %in% c("T","V","X") ~ "Male"),
         Age = case_when(muscle_groups_cpm_scun %in% c("S","T") ~ "Adult",
                         muscle_groups_cpm_scun %in% c("U","V") ~ "Neonate",
                         muscle_groups_cpm_scun %in% c("W","X") ~ "Maturing")
  )





## Load in pi and TajD

window_pi_chrom1 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/windows/NC_056522.1_10kb_pi.windowed.pi", sep='\t', header=TRUE)
window_pi_chrom2 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/windows/NC_056523.1_10kb_pi.windowed.pi", sep='\t', header=TRUE)
window_pi_chrom3 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/windows/NC_056524.1_10kb_pi.windowed.pi", sep='\t', header=TRUE)
window_pi_chrom4 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/windows/NC_056525.1_10kb_pi.windowed.pi", sep='\t', header=TRUE)
window_pi_chrom5 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/windows/NC_056526.1_10kb_pi.windowed.pi", sep='\t', header=TRUE)
window_pi_chrom6 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/windows/NC_056527.1_10kb_pi.windowed.pi", sep='\t', header=TRUE)
window_pi_chrom7 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/windows/NC_056528.1_10kb_pi.windowed.pi", sep='\t', header=TRUE)
window_pi_chrom8 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/windows/NC_056529.1_10kb_pi.windowed.pi", sep='\t', header=TRUE)
window_pi_chrom9 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/windows/NC_056530.1_10kb_pi.windowed.pi", sep='\t', header=TRUE)
window_pi_chromX <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/windows/NC_056531.1_10kb_pi.windowed.pi", sep='\t', header=TRUE)
window_pi_chrom11 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/windows/NC_056532.1_10kb_pi.windowed.pi", sep='\t', header=TRUE)

allChromWindowPi <- rbind(window_pi_chrom1, window_pi_chrom2, window_pi_chrom3, window_pi_chrom4, window_pi_chrom5, window_pi_chrom6, window_pi_chrom7, window_pi_chrom8, window_pi_chrom9, window_pi_chromX, window_pi_chrom11)

site_pi_chrom1 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/sites/NC_056522.1.sites.pi", sep='\t', header=TRUE)
site_pi_chrom2 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/sites/NC_056523.1.sites.pi", sep='\t', header=TRUE)
site_pi_chrom3 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/sites/NC_056524.1.sites.pi", sep='\t', header=TRUE)
site_pi_chrom4 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/sites/NC_056525.1.sites.pi", sep='\t', header=TRUE)
site_pi_chrom5 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/sites/NC_056526.1.sites.pi", sep='\t', header=TRUE)
site_pi_chrom6 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/sites/NC_056527.1.sites.pi", sep='\t', header=TRUE)
site_pi_chrom7 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/sites/NC_056528.1.sites.pi", sep='\t', header=TRUE)
site_pi_chrom8 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/sites/NC_056529.1.sites.pi", sep='\t', header=TRUE)
site_pi_chrom9 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/sites/NC_056530.1.sites.pi", sep='\t', header=TRUE)
site_pi_chromX <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/sites/NC_056531.1.sites.pi", sep='\t', header=TRUE)
site_pi_chrom11 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/pi/sites/NC_056532.1.sites.pi", sep='\t', header=TRUE)

allChromSitePi <- rbind(site_pi_chrom1, site_pi_chrom2, site_pi_chrom3, site_pi_chrom4, site_pi_chrom5, site_pi_chrom6, site_pi_chrom7, site_pi_chrom8, site_pi_chrom9, site_pi_chromX, site_pi_chrom11)

TajD_chrom1 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/TajD/windows/NC_056522.1_10kb_TajD.Tajima.D", sep='\t', header=TRUE)
TajD_chrom2 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/TajD/windows/NC_056523.1_10kb_TajD.Tajima.D", sep='\t', header=TRUE)
TajD_chrom3 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/TajD/windows/NC_056524.1_10kb_TajD.Tajima.D", sep='\t', header=TRUE)
TajD_chrom4 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/TajD/windows/NC_056525.1_10kb_TajD.Tajima.D", sep='\t', header=TRUE)
TajD_chrom5 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/TajD/windows/NC_056526.1_10kb_TajD.Tajima.D", sep='\t', header=TRUE)
TajD_chrom6 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/TajD/windows/NC_056527.1_10kb_TajD.Tajima.D", sep='\t', header=TRUE)
TajD_chrom7 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/TajD/windows/NC_056528.1_10kb_TajD.Tajima.D", sep='\t', header=TRUE)
TajD_chrom8 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/TajD/windows/NC_056529.1_10kb_TajD.Tajima.D", sep='\t', header=TRUE)
TajD_chrom9 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/TajD/windows/NC_056530.1_10kb_TajD.Tajima.D", sep='\t', header=TRUE)
TajD_chromX <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/TajD/windows/NC_056531.1_10kb_TajD.Tajima.D", sep='\t', header=TRUE)
TajD_chrom11 <- read.delim(file="/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/TajD/windows/NC_056532.1_10kb_TajD.Tajima.D", sep='\t', header=TRUE)

allChromTajD <- rbind(TajD_chrom1, TajD_chrom2, TajD_chrom3, TajD_chrom4, TajD_chrom5, TajD_chrom6, TajD_chrom7, TajD_chrom8, TajD_chrom9, TajD_chromX, TajD_chrom11)


chromosome_names <- read.csv("/scratch/dtn2an/workspace/refseq_aligned/SceUnd1.1_gene_Chr_locations.csv", header=TRUE)

chromosome_names$window_pi <- NA
for (i in 1:NROW(chromosome_names)){
  for(j in 1:NROW(allChromWindowPi)){
   if ((chromosome_names$chromosome[i] == allChromWindowPi$CHROM[j]) & (chromosome_names$start[i] > allChromWindowPi$BIN_START[j]) & (allChromWindowPi$BIN_END[j] > chromosome_names$stop[i])){
    chromosome_names$window_pi[i] <- allChromWindowPi$PI[j]
    break
    }else{
    }
  }
}

chromosome_names$site_pi <- NA
for (i in 1:NROW(chromosome_names)){
  for(j in 1:NROW(allChromSitePi)){
    if ((chromosome_names$chromosome[i] == allChromSitePi$CHROM[j]) & (chromosome_names$start[i] < allChromSitePi$POS[j]) & (allChromSitePi$POS[j] < chromosome_names$stop[i])){
      chromosome_names$site_pi[i] <- allChromSitePi$PI[j]
      break
    }else{
    }
  }
}

chromosome_names$TajD <- NA
for (i in 1:NROW(chromosome_names)){
  for(j in 1:NROW(allChromTajD)){
    if ((chromosome_names$chromosome[i] == allChromTajD$CHROM[j]) & (chromosome_names$start[i] > allChromTajD$BIN_START[j]) & (allChromTajD$BIN_START[j]+9999 > chromosome_names$stop[i])){
      chromosome_names$TajD[i] <- allChromTajD$TajimaD[j]
      break
    }else{
    }
  }
}

colnames(chromosome_names) <- c("gene_names", "chromosome", "start", "stop", "window_pi", "site_pi", "TajD")

write.csv(chromosome_names, "/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/Pi.TajD.gene.names.csv", sep=',')



pi.TajD.F <- left_join(chromosome_names, all_fvalues, by="gene_names")

write.csv(pi.TajD.F, "/scratch/dtn2an/workspace/refseq_aligned/liver_snp_filtered/main_chroms/Pi.TajD.F.csv", sep=',')






