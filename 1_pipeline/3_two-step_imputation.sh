#!/bin/bash
####
# Two-step imputation
####

# directories
mkdir "${workdir}"/Imp/"${panel}"/{09_Variants,10_LikeliUpdated,11_Imputed,12_TargetPop}

# paths
workdir='path/to/directory'
panel='internal' # internal, external, combined or diverse
chrom='chromosome_name'
ref='path/to/reference'


# parameters

# 06 - variant calling
module load bcftools/1.11
find ${workdir}/PrePro/${dataset}/03_MapToRef/*map2host.bam > ${workdir}/PrePro/${dataset}/03_MapToRef/samplefile.txt
sed -i -e "s#^#${workdir}/PrePro/${dataset}/03_MapToRef/#" samplefile.txt.txt

bcftools mpileup -C 50 -q 30 -Q 20 -Ou -f ${ref} -r ${chrom} -b samplefile.txt | bcftools call -o ${workdir}/Imp/${panel}/09_Variants/${chrom}_all.vcf.gz -m -v -Oz
bcftools view ${workdir}/Imp/${panel}/09_Variants/${chrom}_all.vcf.gz -o ${workdir}/Imp/${panel}/09_Variants/${chrom}_snps.vcf.gz -m2 -M2 -v snps -Oz

# 07 - likelihood update
module load java/1.8.0 tabix/1.2.1 bcftools/1.11 anaconda3/4.4.0
java -Xss5m -Xmx180g -jar /services/tools/beagle/4.1/beagle.27Jul16.86a.jar gl=${workdir}/Imp/${panel}/09_Variants/${chrom}_snps.vcf.gz out=${workdir}/Imp/${panel}/09_Variants/${chrom}_prob ref=${workdir}/RefPanel/${panel}/08_Phased/${panel}_phased.vcf.gz gprobs=true chrom=${chrom}
bcftools index ${workdir}/04_Probs/${panel}/${chrom}_prob.vcf.gz
bcftools +setGT ${workdir}/04_Probs/${panel}/${chrom}_prob.vcf.gz -- -t q -n . -e'FORMAT/GP>=0.99' > ${workdir}/04_Probs/${panel}/${chrom}_prob_filt.vcf
bgzip -f ${workdir}/04_Probs/${panel}/${chrom}_prob_filt.vcf



# 08 - imputation
module load java/1.8.0 tabix/1.2.1 bcftools/1.11 anaconda3/4.4.0
tabix ${workdir}/04_Probs/${panel}/${chrom}_prob_filt.vcf.gz
java -Xmx180g -jar /services/tools/beagle/5.1/beagle-5.1.jar gt=${workdir}/04_Probs/${panel}/${chrom}_prob_filt.vcf.gz ref=${panel}/${panel}_phased.vcf.gz gp=true chrom=${chrom} out=${workdir}/05_Imp/${panel}/${chrom}_imputed
bcftools index ${workdir}/05_Imp/${panel}/${chrom}_imputed.vcf.gz
bcftools +setGT ${workdir}/05_Imp/${panel}/${chrom}_imputed.vcf.gz -- -t q -n . -e'FORMAT/GP>=0.99' >  ${workdir}/05_Imp/${panel}/${chrom}_imputed_filt.vcf
bgzip -f ${workdir}/05_Imp/${panel}/${chrom}_imputed_filt.vcf
