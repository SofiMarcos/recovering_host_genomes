#!/bin/bash
####
# Custom reference panel design
####

# paths
workdir='path/to/directory'

# directories
mkdir "${workdir}"/RefPanel/{06_Variants,07_PanelDesign,08_Phased}
mkdir "${workdir}"/RefPanel/06_Variants/"${dataset}"
mkdir "${workdir}"/RefPanel/07_PanelDesign/"${panel}"
mkdir "${workdir}"/RefPanel/08_Phased/"${panel}"

# parameters
threads='40'

# 06 - variant calling
### A - Internal reference samples
##### paths
dataset1='type_of_dataset' # internal reference samples
ref='path/to/reference'
sample='sample_codes'
chrom='chromosome_name'

##### divide internal reference samples if needed
samtools view -h -b ${workdir}/PrePro/${dataset}/03_MapToRef/${sample}_map2host.bam ${chrom} > ${workdir}/PrePro/${dataset}/03_MapToRef/${sample}_${chrom}.bam
samtools index ${workdir}/PrePro/${dataset}/03_MapToRef/${sample}_${chrom}.bam

##### variant calling
find ${workdir}/PrePro/${dataset}/03_MapToRef/*${chrom}.bam > ${workdir}/PrePro/${dataset}/03_MapToRef/${chrom}.txt
sed -i -e "s#^#${workdir}/PrePro/${dataset}/03_MapToRef/#" ${chrom}.txt

bcftools mpileup -b ${workdir}/PrePro/${dataset}/03_MapToRef/${chrom}.txt -C 50 -q 30 -Q 20 -Ou -f ${ref} -r ${chrom} | bcftools call -m -v -Ou | bcftools filter -o ${workdir}/RefPanel/06_Variants/${dataset1}/${chrom}_all.vcf.gz -s LowQual -e '%QUAL<30 || DP<(AVG(DP)*3)' -Oz --threads ${threads}
bcftools view ${workdir}/RefPanel/06_Variants/${dataset1}/${chrom}_all.vcf.gz -o ${workdir}/RefPanel/06_Variants/${dataset1}/${chrom}_snps.vcf.gz -m2 -M2 -v snps -Oz --threads ${threads}
rm *${chrom}.bam

### B - External reference samples
##### paths
dataset2='type_of_dataset' # external reference samples
ref='path/to/reference'
chrom='chromosome_name'

##### download from database
module load anaconda3/4.4.0 enabrowsertools/1.6
enaDataGet -f fastq -d ${workdir}/PrePro/00_RawData/${sample}
rm ${workdir}/PrePro/00_RawData/${sample}.fastq.gz

#### alignment
module load samtools/1.11 bwa/0.7.16a
bwa mem -t 8 -R "@RG\\tID:tID\\tCN:tCN\\tDS:tDS\\tPL:tPL\\tSM:${sample}" ${ref} ${workdir}/PrePro/00_RawData/${sample}_1.fastq.gz ${workdir}/PrePro/00_RawData/${sample}_2.fastq.gz | samtools sort -o ${workdir}/PrePro/${dataset}/03_MapToRef/${SAMPLE}.bam --threads ${threads}

##### genotyping
module load java/1.8.0 gatk/4.0.8.1
find ${workdir}/PrePro/${dataset}/03_MapToRef/*.bam > ${workdir}/PrePro/${dataset}/03_MapToRef/samplefile.txt
sed -i -e "s#^#${workdir}/PrePro/${dataset}/03_MapToRef/#" samplefile.txt

gatk IndexFeaturefile -F ${workdir}/RefPanel/06_Variants/${dataset2}/${chrom}_snps.vcf.gz
gatk HaplotypeCaller --java-options "-Xmx180g" -I ${workdir}/PrePro/${dataset}/03_MapToRef/samplefile.txt --output ${workdir}/RefPanel/06_Variants/${dataset2}/${chrom}_all.vcf.gz --reference ${ref} --alleles ${workdir}/RefPanel/int/06_Variants/${chrom}_snps.vcf.gz --intervals ${workdir}/RefPanel/int/06_Variants/${chrom}_snps.vcf.gz --sample-ploidy 2 --min-base-quality-score 20 --standard-min-confidence-threshold-for-calling 30.0
gatk SelectVariants -R ${ref} -V ${workdir}/RefPanel/06_Variants/${dataset2}/${chrom}_all.vcf.gz  -O ${workdir}/RefPanel/06_Variants/${dataset2}/${chrom}_snps.vcf.gz --select-type-to-include SNP


# 07 - combining datasets to design different reference panels (diverse panel as an example)
module load plink2/1.90beta6.17 bcftools/1.11
bcftools merge ${workdir}/RefPanel/06_Variants/${dataset1}/${chrom}_snps.vcf.gz ${workdir}/RefPanel/06_Variants/${dataset2}/${chrom}_snps.vcf.gz -o ${workdir}/RefPanel/07_PanelDesign/${panel}/${chrom}_snps.vcf.gz -Oz --threads ${threads}
bcftools index ${workdir}/RefPanel/07_PanelDesign/${panel}/${chrom}_snps.vcf.gz
plink --vcf ${workdir}/RefPanel/07_PanelDesign/${panel}/${chrom}_snps.vcf.gz --out ${workdir}/RefPanel/07_PanelDesign/${panel}/${chrom}_snps --double-id --make-bed --allow-extra-chr --keep-allele-order  --real-ref-alleles --set-missing-var-ids '@:#\$1,\$2'
plink --bfile ${workdir}/RefPanel/07_PanelDesign/${panel}/${chrom}_snps --out ${workdir}/RefPanel/07_PanelDesign/${panel}/${chrom}_snps_filt --double-id --allow-extra-chr --keep-allele-order  --real-ref-alleles --geno 0 --recode vcf-iid bgz
find ${workdir}/RefPanel/07_PanelDesign/${panel}/ -name "*nosex" -delete

# 08 - phasing
## paths
panel='some' # internal, external, combined or diverse

module load tabix/1.2.1 shapeit4/4.1.3
tabix ${workdir}/RefPanel/07_PanelDesign/${panel}/${chrom}_snps_filt.vcf.gz
shapeit4 --input ${workdir}/RefPanel/07_PanelDesign/${panel}/${chrom}_snps_filt.vcf.gz --output ${workdir}/RefPanel/08_Phased/${panel}/${chrom}_phased.vcf.gz --region ${chrom}
tabix ${workdir}/RefPanel/08_Phased/${panel}/${chrom}_phased.vcf.gz
