######
### leave-one-out cross-validation
######

#!/bin/bash
module load bcftools/1.11
module load tabix/1.2.1 shapeit4/4.1.3
module load java/1.8.0 beagle/4.1
module load beagle/5.1
module load perl/5.24.0 vcftools/0.1.16

# defining files
CHROM='chromosome_name'

# delete 1 of the validation samples from the reference panel
for SAMPLE in CA19_13 CA19_18 CA24_01 CC19_18 CA03_14 CA06_05 CA13_03 CA21_09 CB20_03 CC09_13 CB07_06 CC20_14
do
 bcftools view -s ^${SAMPLE}C1a -Oz -o ${SAMPLE}_miss_${CHROM}.vcf.gz HD_${CHROM}.vcf.gz
 bcftools index ${SAMPLE}_miss_${CHROM}.vcf.gz
done

# phase the reference panels
for SAMPLE in CA19_13 CA19_18 CA24_01 CC19_18 CA03_14 CA06_05 CA13_03 CA21_09 CB20_03 CC09_13 CB07_06 CC20_14
do
shapeit4 --input ${SAMPLE}_miss.vcf.gz --region ${CHROM} --output ${SAMPLE}_filt_phased.vcf.gz
tabix ${SAMPLE}_filt_phased.vcf.gz
done


# likelihoods update
for SAMPLE in CA19_13 CA19_18 CA24_01 CC19_18 CA03_14 CA06_05 CA13_03 CA21_09 CB20_03 CC09_13 CB07_06 CC20_14
do
java -Xmx180g -jar /services/tools/beagle/4.1/beagle.27Jul16.86a.jar gl=LD_${CHROM}_SNPs.vcf.gz map=imputing.map ref=${SAMPLE}_filt_phased.vcf.gz gprobs=true chrom=${CHROM} out=${SAMPLE}_prob_${CHROM}
bcftools index ${SAMPLE}_prob_${CHROM}.vcf.gz
bcftools +setGT ${SAMPLE}_prob_${CHROM}.vcf.gz -- -t q -n . -e'FORMAT/GP>=0.99' > ${SAMPLE}_prob_${CHROM}_filt.vcf
module unload bcftools/1.11
module load anaconda3/4.4.0
bgzip ${SAMPLE}_prob_${CHROM}_filt.vcf
module unload anaconda3/4.4.0
module load bcftools/1.11
done

# impute missing variants
for SAMPLE in CA19_13 CA19_18 CA24_01 CC19_18 CA03_14 CA06_05 CA13_03 CA21_09 CB20_03 CC09_13 CB07_06 CC20_14
do
  tabix ${SAMPLE}_prob_${CHROM}_filt.vcf.gz
  java -Xmx180g -jar /services/tools/beagle/5.1/beagle-5.1.jar gt=${SAMPLE}_prob_${CHROM}_filt.vcf.gz ref=${SAMPLE}_filt_phased.vcf.gz gp=true chrom=${CHROM} out=${SAMPLE}_imputed_${CHROM}
  tabix ${SAMPLE}_imputed_${CHROM}.vcf.gz
  module load bcftools/1.11
  bcftools +setGT ${SAMPLE}_imputed_${CHROM}.vcf.gz -- -t q -n . -e'FORMAT/GP>=0.99' >  ${SAMPLE}_imputed__${CHROM}_filt.vcf
  module unload bcftools/1.11
  module load anaconda3/4.4.0
  bgzip ${SAMPLE}_imputed__${CHROM}_filt.vcf
  module unload anaconda3/4.4.0
  module load bcftools/1.11
done

# calculate concordance tables
for SAMPLE in CA19_13 CA19_18 CA24_01 CC19_18 CA03_14 CA06_05 CA13_03 CA21_09 CB20_03 CC09_13 CB07_06 CC20_14
do
  bcftools view -s ${SAMPLE}C1a -Oz -o concor_probs/${SAMPLE}C1a.vcf.gz originals/Combined_filt.vcf.gz
  bcftools view -s ${SAMPLE}F1a -Oz -o concor_probs/${SAMPLE}F1a.vcf.gz ${SAMPLE}_prob_clean.vcf.gz
  bcftools view -h concor_probs/${SAMPLE}C1a.vcf.gz > concor_probs/${SAMPLE}.header
  bcftools reheader -h concor_probs/${SAMPLE}.header -o concor_probs/${SAMPLE}F1a_mod.vcf.gz concor_probs/${SAMPLE}F1a.vcf.gz
  vcftools --gzvcf concor_probs/${SAMPLE}C1a.vcf.gz --gzdiff concor_probs/${SAMPLE}F1a_mod.vcf.gz --diff-discordance-matrix --out concor_probs/${SAMPLE}_probs
done

# continue in R
