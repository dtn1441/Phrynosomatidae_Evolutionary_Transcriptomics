#!/bin/bash
#
#
#SBATCH -J merge_tissues
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH -t 24:00:00
#SBATCH --mem 15G
#SBATCH -o /scratch/dtn2an/workspace/refseq_aligned/stderr/merge_tissues.out
#SBATCH -e /scratch/dtn2an/workspace/refseq_aligned/stderr/merge_tissues.err
#SBATCH --account=coxlab
#SBATCH --partition=standard

cd /scratch/dtn2an/workspace/refseq_aligned/merge_runs/



module load picard
module load samtools


samtools_n=/scratch/dtn2an/tools/samtools-1.16.1/samtools



java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_12_merged_sorted.bam \
O=Scun_12_merged_sorted.fix.bam \
RGLB="AF3" \
RGPL="AF3" \
RGPU="AF3" \
RGSM="AF3"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_14_merged_sorted.bam \
O=Scun_14_merged_sorted.fix.bam \
RGLB="AF3" \
RGPL="AF3" \
RGPU="AF3" \
RGSM="AF3"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_16_merged_sorted.bam \
O=Scun_16_merged_sorted.fix.bam \
RGLB="AF3" \
RGPL="AF3" \
RGPU="AF3" \
RGSM="AF3"

samtools merge \
AF3.bam \
Scun_12_merged_sorted.fix.bam  \
Scun_14_merged_sorted.fix.bam \
Scun_16_merged_sorted.fix.bam

#############
java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_18_merged_sorted.bam \
O=Scun_18_merged_sorted.fix.bam \
RGLB="AF4" \
RGPL="AF4" \
RGPU="AF4" \
RGSM="AF4"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_20_merged_sorted.bam \
O=Scun_20_merged_sorted.fix.bam \
RGLB="AF4" \
RGPL="AF4" \
RGPU="AF4" \
RGSM="AF4"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_22_merged_sorted.bam \
O=Scun_22_merged_sorted.fix.bam \
RGLB="AF4" \
RGPL="AF4" \
RGPU="AF4" \
RGSM="AF4"

samtools merge \
AF4.bam \
Scun_18_merged_sorted.fix.bam  \
Scun_20_merged_sorted.fix.bam \
Scun_22_merged_sorted.fix.bam

################

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_23_merged_sorted.bam \
O=Scun_23_merged_sorted.fix.bam \
RGLB="AM5" \
RGPL="AM5" \
RGPU="AM5" \
RGSM="AM5"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_25_merged_sorted.bam \
O=Scun_25_merged_sorted.fix.bam \
RGLB="AM5" \
RGPL="AM5" \
RGPU="AM5" \
RGSM="AM5"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_27_merged_sorted.bam \
O=Scun_27_merged_sorted.fix.bam \
RGLB="AM5" \
RGPL="AM5" \
RGPU="AM5" \
RGSM="AM5"

samtools merge \
AM5.bam \
Scun_23_merged_sorted.fix.bam  \
Scun_25_merged_sorted.fix.bam \
Scun_27_merged_sorted.fix.bam


####################

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_28_merged_sorted.bam \
O=Scun_28_merged_sorted.fix.bam \
RGLB="AM6" \
RGPL="AM6" \
RGPU="AM6" \
RGSM="AM6"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_30_merged_sorted.bam \
O=Scun_30_merged_sorted.fix.bam \
RGLB="AM6" \
RGPL="AM6" \
RGPU="AM6" \
RGSM="AM6"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_32_merged_sorted.bam \
O=Scun_32_merged_sorted.fix.bam \
RGLB="AM6" \
RGPL="AM6" \
RGPU="AM6" \
RGSM="AM6"

samtools merge \
AM6.bam \
Scun_28_merged_sorted.fix.bam  \
Scun_30_merged_sorted.fix.bam \
Scun_32_merged_sorted.fix.bam



###############

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_33_merged_sorted.bam \
O=Scun_33_merged_sorted.fix.bam \
RGLB="AF7" \
RGPL="AF7" \
RGPU="AF7" \
RGSM="AF7"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_35_merged_sorted.bam \
O=Scun_35_merged_sorted.fix.bam \
RGLB="AF7" \
RGPL="AF7" \
RGPU="AF7" \
RGSM="AF7"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_37_merged_sorted.bam \
O=Scun_37_merged_sorted.fix.bam \
RGLB="AF7" \
RGPL="AF7" \
RGPU="AF7" \
RGSM="AF7"

samtools merge \
AF7.bam \
Scun_33_merged_sorted.fix.bam  \
Scun_35_merged_sorted.fix.bam \
Scun_37_merged_sorted.fix.bam

#####################

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_38_merged_sorted.bam \
O=Scun_38_merged_sorted.fix.bam \
RGLB="AM8" \
RGPL="AM8" \
RGPU="AM8" \
RGSM="AM8"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_40_merged_sorted.bam \
O=Scun_40_merged_sorted.fix.bam \
RGLB="AM8" \
RGPL="AM8" \
RGPU="AM8" \
RGSM="AM8"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_42_merged_sorted.bam \
O=Scun_42_merged_sorted.fix.bam \
RGLB="AM8" \
RGPL="AM8" \
RGPU="AM8" \
RGSM="AM8"

samtools merge \
AM8.bam \
Scun_38_merged_sorted.fix.bam  \
Scun_40_merged_sorted.fix.bam \
Scun_42_merged_sorted.fix.bam

#############

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_43_merged_sorted.bam \
O=Scun_43_merged_sorted.fix.bam \
RGLB="AM9" \
RGPL="AM9" \
RGPU="AM9" \
RGSM="AM9"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_45_merged_sorted.bam \
O=Scun_45_merged_sorted.fix.bam \
RGLB="AM9" \
RGPL="AM9" \
RGPU="AM9" \
RGSM="AM9"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_47_merged_sorted.bam \
O=Scun_47_merged_sorted.fix.bam \
RGLB="AM9" \
RGPL="AM9" \
RGPU="AM9" \
RGSM="AM9"

samtools merge \
AM9.bam \
Scun_43_merged_sorted.fix.bam  \
Scun_45_merged_sorted.fix.bam \
Scun_47_merged_sorted.fix.bam

####################

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_48_merged_sorted.bam \
O=Scun_48_merged_sorted.fix.bam \
RGLB="YM21" \
RGPL="YM21" \
RGPU="YM21" \
RGSM="YM21"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_50_merged_sorted.bam \
O=Scun_50_merged_sorted.fix.bam \
RGLB="YM21" \
RGPL="YM21" \
RGPU="YM21" \
RGSM="YM21"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_52_merged_sorted.bam \
O=Scun_52_merged_sorted.fix.bam \
RGLB="YM21" \
RGPL="YM21" \
RGPU="YM21" \
RGSM="YM21"

samtools merge \
YM21.bam \
Scun_48_merged_sorted.fix.bam  \
Scun_50_merged_sorted.fix.bam \
Scun_52_merged_sorted.fix.bam


#############

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_54_merged_sorted.bam \
O=Scun_54_merged_sorted.fix.bam \
RGLB="AF11" \
RGPL="AF11" \
RGPU="AF11" \
RGSM="AF11"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_56_merged_sorted.bam \
O=Scun_56_merged_sorted.fix.bam \
RGLB="AF11" \
RGPL="AF11" \
RGPU="AF11" \
RGSM="AF11"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_58_merged_sorted.bam \
O=Scun_58_merged_sorted.fix.bam \
RGLB="AF11" \
RGPL="AF11" \
RGPU="AF11" \
RGSM="AF11"

samtools merge \
AF11.bam \
Scun_54_merged_sorted.fix.bam  \
Scun_56_merged_sorted.fix.bam \
Scun_58_merged_sorted.fix.bam

##################

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_60_merged_sorted.bam \
O=Scun_60_merged_sorted.fix.bam \
RGLB="AM12" \
RGPL="AM12" \
RGPU="AM12" \
RGSM="AM12"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_62_merged_sorted.bam \
O=Scun_62_merged_sorted.fix.bam \
RGLB="AM12" \
RGPL="AM12" \
RGPU="AM12" \
RGSM="AM12"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_64_merged_sorted.bam \
O=Scun_64_merged_sorted.fix.bam \
RGLB="AM12" \
RGPL="AM12" \
RGPU="AM12" \
RGSM="AM12"

samtools merge \
AM12.bam \
Scun_60_merged_sorted.fix.bam  \
Scun_62_merged_sorted.fix.bam \
Scun_64_merged_sorted.fix.bam

############################



java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_66_merged_sorted.bam \
O=Scun_66_merged_sorted.fix.bam \
RGLB="AF13" \
RGPL="AF13" \
RGPU="AF13" \
RGSM="AF13"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_68_merged_sorted.bam \
O=Scun_68_merged_sorted.fix.bam \
RGLB="AF13" \
RGPL="AF13" \
RGPU="AF13" \
RGSM="AF13"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_70_merged_sorted.bam \
O=Scun_70_merged_sorted.fix.bam \
RGLB="AF13" \
RGPL="AF13" \
RGPU="AF13" \
RGSM="AF13"

samtools merge \
AF13.bam \
Scun_66_merged_sorted.fix.bam  \
Scun_68_merged_sorted.fix.bam \
Scun_70_merged_sorted.fix.bam

##################

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_72_merged_sorted.bam \
O=Scun_72_merged_sorted.fix.bam \
RGLB="AM14" \
RGPL="AM14" \
RGPU="AM14" \
RGSM="AM14"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_74_merged_sorted.bam \
O=Scun_74_merged_sorted.fix.bam \
RGLB="AM14" \
RGPL="AM14" \
RGPU="AM14" \
RGSM="AM14"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_76_merged_sorted.bam \
O=Scun_76_merged_sorted.fix.bam \
RGLB="AM14" \
RGPL="AM14" \
RGPU="AM14" \
RGSM="AM14"

samtools merge \
AM14.bam \
Scun_72_merged_sorted.fix.bam  \
Scun_74_merged_sorted.fix.bam \
Scun_76_merged_sorted.fix.bam

###################



java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_78_merged_sorted.bam \
O=Scun_78_merged_sorted.fix.bam \
RGLB="AM15" \
RGPL="AM15" \
RGPU="AM15" \
RGSM="AM15"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_80_merged_sorted.bam \
O=Scun_80_merged_sorted.fix.bam \
RGLB="AM15" \
RGPL="AM15" \
RGPU="AM15" \
RGSM="AM15"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_82_merged_sorted.bam \
O=Scun_82_merged_sorted.fix.bam \
RGLB="AM15" \
RGPL="AM15" \
RGPU="AM15" \
RGSM="AM15"

samtools merge \
AM15.bam \
Scun_78_merged_sorted.fix.bam  \
Scun_80_merged_sorted.fix.bam \
Scun_82_merged_sorted.fix.bam



##################



java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_101_merged_sorted.bam \
O=Scun_101_merged_sorted.fix.bam \
RGLB="YM21" \
RGPL="YM21" \
RGPU="YM21" \
RGSM="YM21"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_102_merged_sorted.bam \
O=Scun_102_merged_sorted.fix.bam \
RGLB="YM21" \
RGPL="YM21" \
RGPU="YM21" \
RGSM="YM21"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_104_merged_sorted.bam \
O=Scun_104_merged_sorted.fix.bam \
RGLB="YM21" \
RGPL="YM21" \
RGPU="YM21" \
RGSM="YM21"

samtools merge \
YM21.bam \
Scun_101_merged_sorted.fix.bam  \
Scun_102_merged_sorted.fix.bam \
Scun_104_merged_sorted.fix.bam




########

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_105_merged_sorted.bam \
O=Scun_105_merged_sorted.fix.bam \
RGLB="YM22" \
RGPL="YM22" \
RGPU="YM22" \
RGSM="YM22"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_106_merged_sorted.bam \
O=Scun_106_merged_sorted.fix.bam \
RGLB="YM22" \
RGPL="YM22" \
RGPU="YM22" \
RGSM="YM22"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_108_merged_sorted.bam \
O=Scun_108_merged_sorted.fix.bam \
RGLB="YM22" \
RGPL="YM22" \
RGPU="YM22" \
RGSM="YM22"

samtools merge \
YM22.bam \
Scun_105_merged_sorted.fix.bam  \
Scun_106_merged_sorted.fix.bam \
Scun_108_merged_sorted.fix.bam


#########

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_109_merged_sorted.bam \
O=Scun_109_merged_sorted.fix.bam \
RGLB="YM23" \
RGPL="YM23" \
RGPU="YM23" \
RGSM="YM23"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_110_merged_sorted.bam \
O=Scun_110_merged_sorted.fix.bam \
RGLB="YM23" \
RGPL="YM23" \
RGPU="YM23" \
RGSM="YM23"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_112_merged_sorted.bam \
O=Scun_112_merged_sorted.fix.bam \
RGLB="YM23" \
RGPL="YM23" \
RGPU="YM23" \
RGSM="YM23"

samtools merge \
YM23.bam \
Scun_109_merged_sorted.fix.bam  \
Scun_110_merged_sorted.fix.bam \
Scun_112_merged_sorted.fix.bam




#############



java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_113_merged_sorted.bam \
O=Scun_113_merged_sorted.fix.bam \
RGLB="YM24" \
RGPL="YM24" \
RGPU="YM24" \
RGSM="YM24"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_114_merged_sorted.bam \
O=Scun_114_merged_sorted.fix.bam \
RGLB="YM24" \
RGPL="YM24" \
RGPU="YM24" \
RGSM="YM24"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_116_merged_sorted.bam \
O=Scun_116_merged_sorted.fix.bam \
RGLB="YM24" \
RGPL="YM24" \
RGPU="YM24" \
RGSM="YM24"

samtools merge \
YM24.bam \
Scun_113_merged_sorted.fix.bam  \
Scun_114_merged_sorted.fix.bam \
Scun_116_merged_sorted.fix.bam







#####################


java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_117_merged_sorted.bam \
O=Scun_117_merged_sorted.fix.bam \
RGLB="YM25" \
RGPL="YM25" \
RGPU="YM25" \
RGSM="YM25"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_118_merged_sorted.bam \
O=Scun_118_merged_sorted.fix.bam \
RGLB="YM25" \
RGPL="YM25" \
RGPU="YM25" \
RGSM="YM25"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_120_merged_sorted.bam \
O=Scun_120_merged_sorted.fix.bam \
RGLB="YM25" \
RGPL="YM25" \
RGPU="YM25" \
RGSM="YM25"

samtools merge \
YM25.bam \
Scun_117_merged_sorted.fix.bam  \
Scun_118_merged_sorted.fix.bam \
Scun_120_merged_sorted.fix.bam





################################


java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_121_merged_sorted.bam \
O=Scun_121_merged_sorted.fix.bam \
RGLB="YF26" \
RGPL="YF26" \
RGPU="YF26" \
RGSM="YF26"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_122_merged_sorted.bam \
O=Scun_122_merged_sorted.fix.bam \
RGLB="YF26" \
RGPL="YF26" \
RGPU="YF26" \
RGSM="YF26"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_124_merged_sorted.bam \
O=Scun_124_merged_sorted.fix.bam \
RGLB="YF26" \
RGPL="YF26" \
RGPU="YF26" \
RGSM="YF26"

samtools merge \
YF26.bam \
Scun_121_merged_sorted.fix.bam  \
Scun_122_merged_sorted.fix.bam \
Scun_124_merged_sorted.fix.bam




################


java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_125_merged_sorted.bam \
O=Scun_125_merged_sorted.fix.bam \
RGLB="YF27" \
RGPL="YF27" \
RGPU="YF27" \
RGSM="YF27"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_126_merged_sorted.bam \
O=Scun_126_merged_sorted.fix.bam \
RGLB="YF27" \
RGPL="YF27" \
RGPU="YF27" \
RGSM="YF27"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_128_merged_sorted.bam \
O=Scun_128_merged_sorted.fix.bam \
RGLB="YF27" \
RGPL="YF27" \
RGPU="YF27" \
RGSM="YF27"

samtools merge \
YF27.bam \
Scun_125_merged_sorted.fix.bam  \
Scun_126_merged_sorted.fix.bam \
Scun_128_merged_sorted.fix.bam





###############



java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_129_merged_sorted.bam \
O=Scun_129_merged_sorted.fix.bam \
RGLB="YM28" \
RGPL="YM28" \
RGPU="YM28" \
RGSM="YM28"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_130_merged_sorted.bam \
O=Scun_130_merged_sorted.fix.bam \
RGLB="YM28" \
RGPL="YM28" \
RGPU="YM28" \
RGSM="YM28"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_132_merged_sorted.bam \
O=Scun_132_merged_sorted.fix.bam \
RGLB="YM28" \
RGPL="YM28" \
RGPU="YM28" \
RGSM="YM28"

samtools merge \
YM28.bam \
Scun_129_merged_sorted.fix.bam  \
Scun_130_merged_sorted.fix.bam \
Scun_132_merged_sorted.fix.bam




###################


java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_133_merged_sorted.bam \
O=Scun_133_merged_sorted.fix.bam \
RGLB="YF29" \
RGPL="YF29" \
RGPU="YF29" \
RGSM="YF29"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_134_merged_sorted.bam \
O=Scun_134_merged_sorted.fix.bam \
RGLB="YF29" \
RGPL="YF29" \
RGPU="YF29" \
RGSM="YF29"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_136_merged_sorted.bam \
O=Scun_136_merged_sorted.fix.bam \
RGLB="YF29" \
RGPL="YF29" \
RGPU="YF29" \
RGSM="YF29"

samtools merge \
YF29.bam \
Scun_133_merged_sorted.fix.bam  \
Scun_134_merged_sorted.fix.bam \
Scun_136_merged_sorted.fix.bam





#############



java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_137_merged_sorted.bam \
O=Scun_137_merged_sorted.fix.bam \
RGLB="YF31" \
RGPL="YF31" \
RGPU="YF31" \
RGSM="YF31"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_138_merged_sorted.bam \
O=Scun_138_merged_sorted.fix.bam \
RGLB="YF31" \
RGPL="YF31" \
RGPU="YF31" \
RGSM="YF31"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_140_merged_sorted.bam \
O=Scun_140_merged_sorted.fix.bam \
RGLB="YF31" \
RGPL="YF31" \
RGPU="YF31" \
RGSM="YF31"

samtools merge \
YF31.bam \
Scun_137_merged_sorted.fix.bam  \
Scun_138_merged_sorted.fix.bam \
Scun_140_merged_sorted.fix.bam





############


java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_141_merged_sorted.bam \
O=Scun_141_merged_sorted.fix.bam \
RGLB="YF32" \
RGPL="YF32" \
RGPU="YF32" \
RGSM="YF32"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_142_merged_sorted.bam \
O=Scun_142_merged_sorted.fix.bam \
RGLB="YF32" \
RGPL="YF32" \
RGPU="YF32" \
RGSM="YF32"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_144_merged_sorted.bam \
O=Scun_144_merged_sorted.fix.bam \
RGLB="YF32" \
RGPL="YF32" \
RGPU="YF32" \
RGSM="YF32"

samtools merge \
YF32.bam \
Scun_141_merged_sorted.fix.bam  \
Scun_142_merged_sorted.fix.bam \
Scun_144_merged_sorted.fix.bam



###########

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_145_merged_sorted.bam \
O=Scun_145_merged_sorted.fix.bam \
RGLB="YF33" \
RGPL="YF33" \
RGPU="YF33" \
RGSM="YF33"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_146_merged_sorted.bam \
O=Scun_146_merged_sorted.fix.bam \
RGLB="YF33" \
RGPL="YF33" \
RGPU="YF33" \
RGSM="YF33"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_148_merged_sorted.bam \
O=Scun_148_merged_sorted.fix.bam \
RGLB="YF33" \
RGPL="YF33" \
RGPU="YF33" \
RGSM="YF33"

samtools merge \
YF33.bam \
Scun_145_merged_sorted.fix.bam  \
Scun_146_merged_sorted.fix.bam \
Scun_148_merged_sorted.fix.bam


#########


java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_151_merged_sorted.bam \
O=Scun_151_merged_sorted.fix.bam \
RGLB="HF40" \
RGPL="HF40" \
RGPU="HF40" \
RGSM="HF40"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_152_merged_sorted.bam \
O=Scun_152_merged_sorted.fix.bam \
RGLB="HF40" \
RGPL="HF40" \
RGPU="HF40" \
RGSM="HF40"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_153_merged_sorted.bam \
O=Scun_153_merged_sorted.fix.bam \
RGLB="HF40" \
RGPL="HF40" \
RGPU="HF40" \
RGSM="HF40"

samtools merge \
HF40.bam \
Scun_151_merged_sorted.fix.bam  \
Scun_152_merged_sorted.fix.bam \
Scun_153_merged_sorted.fix.bam



#############

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_154_merged_sorted.bam \
O=Scun_154_merged_sorted.fix.bam \
RGLB="HF41" \
RGPL="HF41" \
RGPU="HF41" \
RGSM="HF41"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_155_merged_sorted.bam \
O=Scun_155_merged_sorted.fix.bam \
RGLB="HF41" \
RGPL="HF41" \
RGPU="HF41" \
RGSM="HF41"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_156_merged_sorted.bam \
O=Scun_156_merged_sorted.fix.bam \
RGLB="HF41" \
RGPL="HF41" \
RGPU="HF41" \
RGSM="HF41"

samtools merge \
HF41.bam \
Scun_154_merged_sorted.fix.bam  \
Scun_155_merged_sorted.fix.bam \
Scun_156_merged_sorted.fix.bam


##########

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_157_merged_sorted.bam \
O=Scun_157_merged_sorted.fix.bam \
RGLB="HF42" \
RGPL="HF42" \
RGPU="HF42" \
RGSM="HF42"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_158_merged_sorted.bam \
O=Scun_158_merged_sorted.fix.bam \
RGLB="HF42" \
RGPL="HF42" \
RGPU="HF42" \
RGSM="HF42"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_159_merged_sorted.bam \
O=Scun_159_merged_sorted.fix.bam \
RGLB="HF42" \
RGPL="HF42" \
RGPU="HF42" \
RGSM="HF42"

samtools merge \
HF42.bam \
Scun_157_merged_sorted.fix.bam  \
Scun_158_merged_sorted.fix.bam \
Scun_159_merged_sorted.fix.bam




###############

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_160_merged_sorted.bam \
O=Scun_160_merged_sorted.fix.bam \
RGLB="HM43" \
RGPL="HM43" \
RGPU="HM43" \
RGSM="HM43"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_161_merged_sorted.bam \
O=Scun_161_merged_sorted.fix.bam \
RGLB="HM43" \
RGPL="HM43" \
RGPU="HM43" \
RGSM="HM43"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_162_merged_sorted.bam \
O=Scun_162_merged_sorted.fix.bam \
RGLB="HM43" \
RGPL="HM43" \
RGPU="HM43" \
RGSM="HM43"

samtools merge \
HM43.bam \
Scun_160_merged_sorted.fix.bam  \
Scun_161_merged_sorted.fix.bam \
Scun_162_merged_sorted.fix.bam


#######
java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_163_merged_sorted.bam \
O=Scun_163_merged_sorted.fix.bam \
RGLB="HF44" \
RGPL="HF44" \
RGPU="HF44" \
RGSM="HF44"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_164_merged_sorted.bam \
O=Scun_164_merged_sorted.fix.bam \
RGLB="HF44" \
RGPL="HF44" \
RGPU="HF44" \
RGSM="HF44"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_165_merged_sorted.bam \
O=Scun_165_merged_sorted.fix.bam \
RGLB="HF44" \
RGPL="HF44" \
RGPU="HF44" \
RGSM="HF44"

samtools merge \
HF44.bam \
Scun_163_merged_sorted.fix.bam  \
Scun_164_merged_sorted.fix.bam \
Scun_165_merged_sorted.fix.bam



######
java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_166_merged_sorted.bam \
O=Scun_166_merged_sorted.fix.bam \
RGLB="HM45" \
RGPL="HM45" \
RGPU="HM45" \
RGSM="HM45"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_167_merged_sorted.bam \
O=Scun_167_merged_sorted.fix.bam \
RGLB="HM45" \
RGPL="HM45" \
RGPU="HM45" \
RGSM="HM45"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_168_merged_sorted.bam \
O=Scun_168_merged_sorted.fix.bam \
RGLB="HM45" \
RGPL="HM45" \
RGPU="HM45" \
RGSM="HM45"

samtools merge \
HM45.bam \
Scun_166_merged_sorted.fix.bam  \
Scun_167_merged_sorted.fix.bam \
Scun_168_merged_sorted.fix.bam



##############
java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_169_merged_sorted.bam \
O=Scun_169_merged_sorted.fix.bam \
RGLB="HM46" \
RGPL="HM46" \
RGPU="HM46" \
RGSM="HM46"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_170_merged_sorted.bam \
O=Scun_170_merged_sorted.fix.bam \
RGLB="HM46" \
RGPL="HM46" \
RGPU="HM46" \
RGSM="HM46"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_171_merged_sorted.bam \
O=Scun_171_merged_sorted.fix.bam \
RGLB="HM46" \
RGPL="HM46" \
RGPU="HM46" \
RGSM="HM46"

samtools merge \
HM46.bam \
Scun_169_merged_sorted.fix.bam  \
Scun_170_merged_sorted.fix.bam \
Scun_171_merged_sorted.fix.bam


###########

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_172_merged_sorted.bam \
O=Scun_172_merged_sorted.fix.bam \
RGLB="HM47" \
RGPL="HM47" \
RGPU="HM47" \
RGSM="HM47"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_173_merged_sorted.bam \
O=Scun_173_merged_sorted.fix.bam \
RGLB="HM47" \
RGPL="HM47" \
RGPU="HM47" \
RGSM="HM47"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_174_merged_sorted.bam \
O=Scun_174_merged_sorted.fix.bam \
RGLB="HM47" \
RGPL="HM47" \
RGPU="HM47" \
RGSM="HM47"

samtools merge \
HM47.bam \
Scun_172_merged_sorted.fix.bam  \
Scun_173_merged_sorted.fix.bam \
Scun_174_merged_sorted.fix.bam





#########

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_175_merged_sorted.bam \
O=Scun_175_merged_sorted.fix.bam \
RGLB="HM48" \
RGPL="HM48" \
RGPU="HM48" \
RGSM="HM48"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_176_merged_sorted.bam \
O=Scun_176_merged_sorted.fix.bam \
RGLB="HM48" \
RGPL="HM48" \
RGPU="HM48" \
RGSM="HM48"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_177_merged_sorted.bam \
O=Scun_177_merged_sorted.fix.bam \
RGLB="HM48" \
RGPL="HM48" \
RGPU="HM48" \
RGSM="HM48"

samtools merge \
HM48.bam \
Scun_175_merged_sorted.fix.bam  \
Scun_176_merged_sorted.fix.bam \
Scun_177_merged_sorted.fix.bam



##############
java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_178_merged_sorted.bam \
O=Scun_178_merged_sorted.fix.bam \
RGLB="HF49" \
RGPL="HF49" \
RGPU="HF49" \
RGSM="HF49"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_179_merged_sorted.bam \
O=Scun_179_merged_sorted.fix.bam \
RGLB="HF49" \
RGPL="HF49" \
RGPU="HF49" \
RGSM="HF49"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_180_merged_sorted.bam \
O=Scun_180_merged_sorted.fix.bam \
RGLB="HF49" \
RGPL="HF49" \
RGPU="HF49" \
RGSM="HF49"

samtools merge \
HF49.bam \
Scun_178_merged_sorted.fix.bam  \
Scun_179_merged_sorted.fix.bam \
Scun_180_merged_sorted.fix.bam


##############
java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_181_merged_sorted.bam \
O=Scun_181_merged_sorted.fix.bam \
RGLB="HM50" \
RGPL="HM50" \
RGPU="HM50" \
RGSM="HM50"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_182_merged_sorted.bam \
O=Scun_182_merged_sorted.fix.bam \
RGLB="HM50" \
RGPL="HM50" \
RGPU="HM50" \
RGSM="HM50"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_183_merged_sorted.bam \
O=Scun_183_merged_sorted.fix.bam \
RGLB="HM50" \
RGPL="HM50" \
RGPU="HM50" \
RGSM="HM50"

samtools merge \
HM50.bam \
Scun_181_merged_sorted.fix.bam  \
Scun_182_merged_sorted.fix.bam \
Scun_183_merged_sorted.fix.bam


##########
java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_184_merged_sorted.bam \
O=Scun_184_merged_sorted.fix.bam \
RGLB="HF51" \
RGPL="HF51" \
RGPU="HF51" \
RGSM="HF51"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_185_merged_sorted.bam \
O=Scun_185_merged_sorted.fix.bam \
RGLB="HF51" \
RGPL="HF51" \
RGPU="HF51" \
RGSM="HF51"

java -jar $PICARD AddOrReplaceReadGroups \
I=Scun_186_merged_sorted.bam \
O=Scun_186_merged_sorted.fix.bam \
RGLB="HF51" \
RGPL="HF51" \
RGPU="HF51" \
RGSM="HF51"

samtools merge \
HF51.bam \
Scun_184_merged_sorted.fix.bam  \
Scun_185_merged_sorted.fix.bam \
Scun_186_merged_sorted.fix.bam









