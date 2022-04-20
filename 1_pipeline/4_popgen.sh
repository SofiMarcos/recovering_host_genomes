module load plink2/1.90beta6.17
WORKDIR='/home/projects/ku-cbd/data/HoloFood/Sofi/LD_imp/all-chroms'

mkdir -p new-PopGen/{Het/NucDiv/Pw/Kinship/Fst}

## checking
awk '{print $1,$1}' allRossCobb.txt > allRossCobb.file

# sending jobs:
pop='/home/projects/ku-cbd/data/HoloFood/Sofi/LD_imp/all-chroms'
qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $pop/ld_decay_20.err -o $pop/ld_decay_20.out -l nodes=1:ppn=40,mem=180gb,walltime=0:03:00:00 -N ld_20 $pop/ld_decay.sh

# filtering VCF file:
plink --vcf ${WORKDIR}/05_PopGen/ALL_SNPs.vcf.gz --double-id --make-bed --allow-extra-chr --keep-allele-order  --real-ref-alleles --set-missing-var-ids  '@:#\$1,\$2' --out ${WORKDIR}/05_PopGen/ALL_SNPs
plink --bfile ${WORKDIR}/05_PopGen/ALL_SNPs --double-id --allow-extra-chr --keep-allele-order  --real-ref-alleles --geno 0 --make-bed --out ${WORKDIR}/05_PopGen/ALL_SNPs_filt
plink --bfile ALL_SNPs_filt --double-id --allow-extra-chr --keep-allele-order  --real-ref-alleles --recode vcf-iid bgz --out ${WORKDIR}/05_PopGen/ALL_SNPs_filt

##  0. Nucleotide diversity
#####################
module load perl/5.24.0 vcftools/0.1.16
vcftools --gzvcf ${WORKDIR}/05_PopGen/Fst/ALL_SNPs_filt.vcf.gz --window-pi 40000 --window-pi-step 20000 --keep ${WORKDIR}/05_PopGen/files/LD_r.txt --out Nuc_div/ross
vcftools --gzvcf ${WORKDIR}/05_PopGen/Fst/ALL_SNPs_filt.vcf.gz --window-pi 40000 --window-pi-step 20000 --keep ${WORKDIR}/05_PopGen/files/LD_c.txt --out Nuc_div/cobb
vcftools --gzvcf ${WORKDIR}/05_PopGen/Fst/ALL_SNPs_filt.vcf.gz --window-pi 40000 --window-pi-step 20000 --keep ${WORKDIR}/05_PopGen/files/Q_1.txt --out Nuc_div/Q_1
vcftools --gzvcf ${WORKDIR}/05_PopGen/Fst/ALL_SNPs_filt.vcf.gz --window-pi 40000 --window-pi-step 20000 --keep ${WORKDIR}/05_PopGen/files/Q_2.txt --out Nuc_div/Q_2

Average ± standard deviation of nucleotide diversity estimated in 40kb windows.
##  1. Linkage decay for each breed (couple of chromosomes, all markers)
#####################
# LD with plink
#!/bin/bash
module load plink2/1.90beta6.17
WORKDIR='/home/projects/ku-cbd/data/HoloFood/Sofi/LD_imp/all-chroms'
# --stdout
CHROM='NC_006107.5'
START='829'
STOP='3000011'
plink --bfile ${WORKDIR}/05_PopGen/ALL_SNPs_filt --maf 0.1 --keep ${WORKDIR}/05_PopGen/files/LD_c.file --r2 dprime --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --chr ${CHROM} --from-bp ${START} --to-bp ${STOP} --allow-extra-chr --keep-allele-order --real-ref-alleles --out ${WORKDIR}/05_PopGen/LD_decay/cobb_${CHROM}
awk '{ print $2,$5,$7}' ${WORKDIR}/05_PopGen/LD_decay/cobb_${CHROM}.ld > ${WORKDIR}/05_PopGen/LD_decay/cobb_${CHROM}.edit.ld
awk 'BEGIN { OFS = "\t" } NR == 1 { $4 = "dist" } NR >= 2 { $4 = $2 - $1 } 1' ${WORKDIR}/05_PopGen/LD_decay/cobb_${CHROM}.edit.ld | awk '{ OFS = "\t"} { print $3, $4}' > ${WORKDIR}/05_PopGen/LD_decay/cobb_${CHROM}.OK.ld
rm ${WORKDIR}/05_PopGen/LD_decay/cobb_${CHROM}.edit.ld ${WORKDIR}/05_PopGen/LD_decay/cobb_${CHROM}.ld

plink --bfile ${WORKDIR}/05_PopGen/ALL_SNPs_filt --maf 0.1 --keep ${WORKDIR}/05_PopGen/files/LD_r.file --r2 dprime --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --chr ${CHROM} --from-bp ${START} --to-bp ${STOP} --allow-extra-chr --keep-allele-order --real-ref-alleles --out ${WORKDIR}/05_PopGen/LD_decay/ross_${CHROM}
awk '{ print $2,$5,$7}' ${WORKDIR}/05_PopGen/LD_decay/ross_${CHROM}.ld > ${WORKDIR}/05_PopGen/LD_decay/ross_${CHROM}.edit.ld
awk 'BEGIN { OFS = "\t" } NR == 1 { $4 = "dist" } NR >= 2 { $4 = $2 - $1 } 1' ${WORKDIR}/05_PopGen/LD_decay/ross_${CHROM}.edit.ld | awk '{ OFS = "\t"} { print $3, $4}' > ${WORKDIR}/05_PopGen/LD_decay/ross_${CHROM}.OK.ld
rm ${WORKDIR}/05_PopGen/LD_decay/ross_${CHROM}.edit.ld ${WORKDIR}/05_PopGen/LD_decay/ross_${CHROM}.ld

plink --bfile ${WORKDIR}/05_PopGen/ALL_SNPs_filt --maf 0.1 --keep ${WORKDIR}/05_PopGen/files/Q_1.file --r2 dprime --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --chr ${CHROM} --from-bp ${START} --to-bp ${STOP} --allow-extra-chr --keep-allele-order --real-ref-alleles --out ${WORKDIR}/05_PopGen/LD_decay/q_1_${CHROM}
awk '{ print $2,$5,$7}' ${WORKDIR}/05_PopGen/LD_decay/q_1_${CHROM}.ld > ${WORKDIR}/05_PopGen/LD_decay/q_1_${CHROM}.edit.ld
awk 'BEGIN { OFS = "\t" } NR == 1 { $4 = "dist" } NR >= 2 { $4 = $2 - $1 } 1' ${WORKDIR}/05_PopGen/LD_decay/q_1_${CHROM}.edit.ld | awk '{ OFS = "\t"} { print $3, $4}' > ${WORKDIR}/05_PopGen/LD_decay/q_1_${CHROM}.OK.ld
rm ${WORKDIR}/05_PopGen/LD_decay/q_1_${CHROM}.edit.ld ${WORKDIR}/05_PopGen/LD_decay/q_1_${CHROM}.ld

plink --bfile ${WORKDIR}/05_PopGen/ALL_SNPs_filt --maf 0.1 --keep ${WORKDIR}/05_PopGen/files/Q_2.file --r2 dprime --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --chr ${CHROM} --from-bp ${START} --to-bp ${STOP} --allow-extra-chr --keep-allele-order --real-ref-alleles --out ${WORKDIR}/05_PopGen/LD_decay/q_2_${CHROM}
awk '{ print $2,$5,$7}' ${WORKDIR}/05_PopGen/LD_decay/q_2_${CHROM}.ld > ${WORKDIR}/05_PopGen/LD_decay/q_2_${CHROM}.edit.ld
awk 'BEGIN { OFS = "\t" } NR == 1 { $4 = "dist" } NR >= 2 { $4 = $2 - $1 } 1' ${WORKDIR}/05_PopGen/LD_decay/q_2_${CHROM}.edit.ld | awk '{ OFS = "\t"} { print $3, $4}' > ${WORKDIR}/05_PopGen/LD_decay/q_2_${CHROM}.OK.ld
rm ${WORKDIR}/05_PopGen/LD_decay/q_2_${CHROM}.edit.ld ${WORKDIR}/05_PopGen/LD_decay/q_2_${CHROM}.ld

rm ${WORKDIR}/05_PopGen/LD_decay/*nosex

## LD with vcftools
#!/bin/bash
# module load perl/5.24.0 vcftools/0.1.16
# CHROM='NC_006107.5'
# vcftools --gzvcf ${WORKDIR}/05_PopGen/ALL_SNPs.vcf.gz --maf 0.1 --max-missing 1 --keep ${WORKDIR}/05_PopGen/files/LD_c.txt --chr ${CHROM} --ld-window-bp 1000000 --ld-window 99999  --geno-r2 --min-r2 0 --phased --out ${WORKDIR}/05_PopGen/LD_decay/${CHROM}_c


##  2. PCA with PLINK with ALL SAMPLES (including reference panel breeds)
#####################
# https://speciationgenomics.github.io/pca/  --> check pruning details
# linkage pruning
    # --indep-pairwise
    # The first argument, 50 denotes we have set a window of 50 Kb.
    # The second argument, 10 is our window step size - meaning we move 10 bp each time we calculate linkage.
    # Finally, we set an r2 threshold - i.e. the threshold of linkage we are willing to tolerate.
plink --bfile ${WORKDIR}/05_PopGen/ALL_SNPs_filt --double-id --allow-extra-chr --keep-allele-order  --real-ref-alleles --indep-pairwise 100 50 0.1 --out ${WORKDIR}/05_PopGen/PCA/ALL_SNPs_prunned_test4
# pca
plink --bfile ${WORKDIR}/05_PopGen/ALL_SNPs_filt --double-id --allow-extra-chr --keep-allele-order  --real-ref-alleles  --extract ${WORKDIR}/05_PopGen/PCA/ALL_SNPs_prunned_test4.prune.in --make-bed --pca --out ${WORKDIR}/05_PopGen/PCA/prune_pca_test4

## HeatMap - Distance square
#####################
plink --bfile ${WORKDIR}/05_PopGen/ALL_SNPs_filt --maf 0.1 --keep ${WORKDIR}/05_PopGen/files/LD_r.file --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --distance square --out ${WORKDIR}/05_PopGen/HeatMap/Ross
plink --bfile ${WORKDIR}/05_PopGen/ALL_SNPs_filt --maf 0.1 --keep ${WORKDIR}/05_PopGen/files/LD_c.file --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --distance square --out ${WORKDIR}/05_PopGen/HeatMap/Cobb

plink --bfile ${WORKDIR}/05_PopGen/ALL_SNPs_filt --maf 0.1 --keep ${WORKDIR}/05_PopGen/files/allRossCobb.file --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --distance square ibs --out ${WORKDIR}/05_PopGen/HeatMap/LD_HD


plink --bfile ${WORKDIR}/05_PopGen/ALL_SNPs_filt --maf 0.1 --double-id --allow-extra-chr --keep-allele-order  --extract ${WORKDIR}/05_PopGen/PCA/ALL_SNPs_prunned_test4.prune.in --real-ref-alleles --distance square ibs --out ${WORKDIR}/05_PopGen/HeatMap/all
##  3. Heterocigosity per individual with VCFtools. Each breed individually.
#####################
# perform Het for all breeds
plink --bfile ${WORKDIR}/05_PopGen/ALL_SNPs_filt --maf 0.1 --keep ${WORKDIR}/05_PopGen/files/LD_r.file --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --het --out ${WORKDIR}/05_PopGen/Het/Ross
plink --bfile ${WORKDIR}/05_PopGen/ALL_SNPs_filt --maf 0.1 --keep ${WORKDIR}/05_PopGen/files/LD_c.file --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --het --out ${WORKDIR}/05_PopGen/Het/Cobb
plink --bfile ${WORKDIR}/05_PopGen/ALL_SNPs_filt --maf 0.1 --keep ${WORKDIR}/05_PopGen/files/Q_1.file --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --het --out ${WORKDIR}/05_PopGen/Het/Q_1
plink --bfile ${WORKDIR}/05_PopGen/ALL_SNPs_filt --maf 0.1 --keep ${WORKDIR}/05_PopGen/files/Q_2.file --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --het --out ${WORKDIR}/05_PopGen/Het/Q_2


plink --bfile ALL_SNPs_diverse_filt --keep ${WORKDIR}/05_PopGen/files/LD_r.file --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --extract pruned.prune.in --het --out het/Ross
plink --bfile ALL_SNPs_diverse_filt --keep ${WORKDIR}/05_PopGen/files/LD_c.file --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --extract pruned.prune.in --het --out het/Cobb
plink --bfile ALL_SNPs_diverse_filt --keep ${WORKDIR}/05_PopGen/files/Q_1.file --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --extract pruned.prune.in --het --out het/Br1
plink --bfile ALL_SNPs_diverse_filt --keep ${WORKDIR}/05_PopGen/files/Q_2.file --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --extract pruned.prune.in --het --out het/Br2
plink --bfile ALL_SNPs_diverse_filt --keep ${WORKDIR}/05_PopGen/files/L1.file --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --extract pruned.prune.in --het --out het/L1
plink --bfile ALL_SNPs_diverse_filt --keep ${WORKDIR}/05_PopGen/files/L2.file --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --extract pruned.prune.in --het --out het/L2
plink --bfile ALL_SNPs_diverse_filt --keep ${WORKDIR}/05_PopGen/files/RJF.file --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --extract pruned.prune.in --het --out het/RJF

##  4. ROH
#####################
# Check out ROH, runs of homozygosity, in plink. I am not familiar with this as in honeybees it doesn’t work well. But it is much used for farm animals.
# https://www.g3journal.org/content/10/12/4615.abstract
# https://www.cambridge.org/core/journals/animal/article/abs/relationship-of-runs-of-homozygosity-with-adaptive-and-production-traits-in-a-paternal-broiler-line/EC2D3150F3D4C7FBC087991EB2FC32AA
# https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6463-x#Sec17
# https://onlinelibrary.wiley.com/doi/pdf/10.1111/jbg.12508?casa_token=pYKa8TvvmTYAAAAA:3ZqueSCkh6Yu8HBWdQzrI2OHKf6ePjbHfLF0onamz9HhqHq62ZhGDt6kwYI7yH8cQQf0ev-lghmvM3ko

# PLINK 1.9
# https://www.cog-genomics.org/plink/1.9/ibd#homozyg
# some studies perform maf and LD prunning (optional).LD prunning for inbreeding populations leads to a biased ROH analysis.
# steps
  # 1.- scanning window defined:
        ## --homozyg-window-snp > scanning window --> increased window size leads to a decrease in estimated Froh
        ## --homozyg-window-het > number of het SNPs
        ## --homozyg-window-missing > a maximal number of missing SNPs

  # 2.- segments of homozygous SNPs are identified
        ## --homozyg--window-threshold > scanning window hit rate

  # 3.- extra constraints set to these homozygous segments
        ## --homozyg-gap > maximal interval between a segment 500kb
        ## --homozyg-het > maximal amount of het allowed in a segment

  # 4.- not filter passed ROH are re-evaluated
        ## --homozyg-density > minimal density requireement
        ## --homozyg-kb
        ## --homozyg-snp
plink --bfile ${WORKDIR}/05_PopGen/ALL_SNPs_filt --maf 0.1 --keep ${WORKDIR}/05_PopGen/files/Q_1.file --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --homozyg --out ${WORKDIR}/05_PopGen/ROH/test
# found 730 ROH
plink --bfile ${WORKDIR}/05_PopGen/ALL_SNPs_filt --maf 0.1 --keep ${WORKDIR}/05_PopGen/files/Q_1.file --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --homozyg --homozyg-window-snp 30 --homozyg-density 1000 --homozyg-kb 10 --homozyg-gap 1000 --out ${WORKDIR}/05_PopGen/ROH/test2
# found 46809 ROH
# BCFTOOLS 1.9
bcftools roh [OPTIONS] file.vcf.gz




##  5. Fst
#####################
# 2 broiler breeds against each other
# 1 breed against all others combined
module load perl/5.24.0 vcftools/0.1.16
vcftools --gzvcf ALL_SNPs_filt.vcf.gz --mac 3 --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop ${WORKDIR}/05_PopGen/files/LD_r.txt --weir-fst-pop ${WORKDIR}/05_PopGen/files/LD_c.txt --out ${WORKDIR}/05_PopGen/Fst/Fst_R_vs_C_40kb_winds
vcftools --gzvcf ${WORKDIR}/05_PopGen/Fst/ALL_SNPs_filt.vcf.gz --mac 3 --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop ${WORKDIR}/05_PopGen/files/LD_r.txt --weir-fst-pop ${WORKDIR}/05_PopGen/files/allQanbari.txt --out ${WORKDIR}/05_PopGen/Fst/Fst_R_vs_Qanbari_40kb_winds
vcftools --gzvcf ${WORKDIR}/05_PopGen/Fst/ALL_SNPs_filt.vcf.gz --mac 3 --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop ${WORKDIR}/05_PopGen/files/LD_c.txt --weir-fst-pop ${WORKDIR}/05_PopGen/files/allQanbari.txt --out ${WORKDIR}/05_PopGen/Fst/Fst_R_vs_Qanbari_40kb_winds

qan='/home/projects/ku-cbd/data/HoloFood/Sofi/LD_imp/all-chroms/05_PopGen/Fst'
qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $qan/fstwindow.err -o $qan/fstwindow.out -l nodes=1:ppn=40,mem=180gb,walltime=0:02:00:00 -N fstwindow $qan/fstwindow.sh


# to write allele frequencies:
less ALL_SNPs_diverse.bim | grep NC_006089.5:149116928
plink --bfile ALL_SNPs_diverse  --keep ${WORKDIR}/05_PopGen/files/LD.file --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --freq --snp 'NC_006089.5:149116928[GRCg6a]\G,\T' --out L
plink --bfile ALL_SNPs_diverse  --keep ${WORKDIR}/05_PopGen/files/RJF.file --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --freq --snp 'NC_006089.5:149116928[GRCg6a]\G,\T' --out R
more L.frq
more R.frq


module load perl/5.24.0 vcftools/0.1.16
vcftools --gzvcf ALL_SNPs_diverse_filt.vcf.gz --mac 3 --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop ${WORKDIR}/05_PopGen/files/LD.txt --weir-fst-pop ${WORKDIR}/05_PopGen/files/RJF.txt --out fst/Fst_R-C_vs_RJF_40kb_winds
vcftools --gzvcf ALL_SNPs_filt.vcf.gz --mac 3 --fst-window-size 40000 --fst-window-step 20000 --weir-fst-pop ${WORKDIR}/05_PopGen/files/LD.txt --weir-fst-pop ${WORKDIR}/05_PopGen/files/allQanbari.txt --out Fst_R-C_vs_broiler_40kb_winds

#Fst analysis
vcftools --gzvcf ${VCF} --weir-fst-pop Cobbs.txt --weir-fst-pop Rosses.txt --fst-window-size interger --out cobb_vs_ross
# filtering with VCFtools
vcftools --gzvcf ALL_SNPs_diverse.vcf.gz --maf 0.1 --max-missing 1 --recode --recode-INFO-all --out ALL_SNPs_filt
module load tabix/1.2.1
gzip ALL_SNPs_filt.vcf

bcftools view -r NC_006100.5:18040001-18120000 -S ${WORKDIR}/05_PopGen/files/LD.txt -Oz -o sweepchr13_third_q99.99.vcf.gz ALL_SNPs_filt.vcf.gz

bcftools view -r NC_006100.5:18840001-18960000 -S ${WORKDIR}/05_PopGen/files/LD.txt -Oz -o sweepchr13_q99.99.vcf.gz ALL_SNPs_filt.vcf.gz
bcftools view -r NC_006100.5:18040001-19000000 -S ${WORKDIR}/05_PopGen/files/LD.txt -Oz -o sweepchr13_q99.9.vcf.gz ALL_SNPs_filt.vcf.gz
bcftools view -r NC_006100.5:15460001-15680000 -S ${WORKDIR}/05_PopGen/files/LD.txt -Oz -o sweepchr13_second_q99.9.vcf.gz ALL_SNPs_filt.vcf.gz

bcftools view -r NC_006090.5:16940001-16980000 -S ${WORKDIR}/05_PopGen/files/LD.txt -Oz -o sweepchr3_q99.9.vcf.gz ALL_SNPs_filt.vcf.gz
bcftools view -r NC_006091.5:1520001-1560000 -S ${WORKDIR}/05_PopGen/files/LD.txt -Oz -o sweepchr4_q99.9.vcf.gz ALL_SNPs_filt.vcf.gz
bcftools view -r NC_006095.5:10580001-10740000 -S ${WORKDIR}/05_PopGen/files/LD.txt -Oz -o sweepchr8_q99.9.vcf.gz ALL_SNPs_filt.vcf.gz
bcftools view -r NC_006092.5:31300001-31340000 -S ${WORKDIR}/05_PopGen/files/LD.txt -Oz -o Annotations/sweepchr5_q99.9.vcf.gz ALL_SNPs_filt.vcf.gz

 bcftools view -H -r NC_028740.2:1920001-1960000 ALL_SNPs_diverse_filt.vcf.gz | awk '{print $1,$2,".",$4,$5,". . ."}' | less

### kinship
plink2/2.00alpha20210203
plink2 --bfile erbel_all --make-king square --out erbel_allKing --allow-extra-chr --keep filt/fathers+workers.keep





##  6. GWAS
#####################
# https://www.youtube.com/watch?v=7QMSZx3io-Q&amp;ab_channel=mathetal
# Quantitative traits (start with a simple indicator on growth or body weight) with population stratification (breed structure)

##  7. hapFLK
#####################
# Input plink binary format (bed, fam, bim). When transforming from vcf to plink, to keep the correct phase you need to include –keep-allele-order.
# You need to manually transform the fam files so that the first column states the breed, and the second the individual ID. https://github.com/bcm-uga/SSMPG2017/blob/master/Presentations/hapflk/hapflk.org
# The version that I used still used phyton2, so on the server I needed to activate phyton2.
# I think the current version is still on phyton2, I cannot find the information in the website.
# But this you can easily test with “hapflk –help”.
# If not even help works, you need to change the pyhton version.
# https://pypi.org/project/hapflk/
          # if not installed in CPH
              # requisites: python >= 2.7, the numpy and scipy packages and the C compiler installed.
              # Unpack the archive (tar -xvzf hapflk-version.tar.gz)
              # cd in the directory created and issue the command:
              # sudo python setup.py install
python -c "import hapflk"
pip install hapflk

# 3 files needed to run hapflk: 1) .ped 2) .map and 3). the kinship matrix estimated from Reynolds distances
# hapmap
hapflk --file hapmap3-lct --kinship kinship.txt -K 15 --nfit=1 --ncpu=2 -p hapmap_tutorial # -K number of clusters to use. --nfit specifies number of EM runs

# building local trees
python local_reynolds.py -p hapmap_tutorial
python local_reynolds.py -p hapmap_tutorial -o hapmap_tutorial -l 135.3e6 -r 136.7e6

# local_trees.R
######################## USER INPUT #########################
### Pre filled with files from the tutorial
## tree file : the whole genome tree as obtained
## from the hapflk analysis
tree_file='hapmap_tree.txt'
## wg_dist_file : the Global Reynolds distances between population,
## as obtained from the hapflk analysis
wg_dist_file='Reynolds.txt'
## Name of the outgroup in the previous file (if any)
outgroup='YRI'
## local reynolds distances for SNPs and Haplotype clusters
## To be computed using the local_reynolds.py script
loc_dist_file_snp='hapflk_snp_reynolds.txt'
loc_dist_file_hap='hapflk_hap_reynolds.txt'
## prefix for outputfiles
prefix='hapmap_tutorial'
##################### END USER INPUT #########################

########## comand for xephh
MGAP=200000
SGAP=20000
WIN=100000
MAX=1000000

# run for each chr
cluster=1
selscan --xpehh --trunc-ok --ehh-win ${WIN} --pmap --max-gap ${MGAP}
--max-extend ${MAX} --gap-scale ${SGAP} --vcf
Cluster1_${CHR}.recode.vcf.gz --vcf-ref
Restofsamples(exceptcluster1).recode.vcf.gz --out
outCluster_${cluster}/chr${CHR} --threads 12;


#############
#     R     #
#############
# modules need to be loaded
module load intel/perflibs/
module load gcc/9.3.0
module load R/4.0.0
# change script mode (only if this file will be ejecuted)
chmod +x {file.r}

allc='/home/projects/ku-cbd/data/HoloFood/Sofi/LD_imp/all-chroms/05_PopGen/LD_decay/chr20'
qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $allc/chr20.err -o $allc/chr20.out -l nodes=1:ppn=40,mem=180gb,walltime=0:04:00:00 -N test4 $allc/chr20.R

##  6. Haplotype blocks
#####################
# checking haplotypes
module load bcftools/1.11
# take the region out
bcftools view -r NC_006100.5:18400000-19000000 -S ${WORKDIR}/05_PopGen/files/allRossCobbordered.txt -Oz -o test ALL_SNPs_filt.vcf.gz --force-samples
# convert to hap file
bcftools convert --hapsample --vcf-ids Hapotype_blocks/sweep13_HD_LD.vcf.gz -o Hapotype_blocks/sweep13_HD_LD



## haplotype diversity
para haplotype diversity usa --hapcount en vcftools. hay que usar bed para las regiones. que es en cada linia chr start end



# Population Genetics Inference
# PLINK
module load plink2/1.90beta6.17
VCF=
BFILE=

## Transforming and filtering ##
plink --vcf ${VCF} --double-id --make-bed --geno 0.1 --mind 0.5 --maf 0.1  --allow-extra-chr --keep-allele-order  --real-ref-alleles --set-missing-var-ids  '@:#[GRCg6a]\$1,\$2' --out ${BFILE}_filt

plink --vcf ${VCF} --double-id --make-bed --allow-extra-chr --keep-allele-order --real-ref-alleles --set-missing-var-ids  '@:#[GRCg6a]\$1,\$2' --out ${BFILE}
plink --bfile ${BFILE} --double-id --make-bed --geno 0.1 --mind 0.5 --maf 0.1  --allow-extra-chr --keep-allele-order --real-ref-alleles --out ${BFILE}_filt
# check missing values
plink --bfile ${BFILE} --double-id --keep-allele-order --allow-extra-chr --keep-allele-order --missing --out ${BFILE}
# eso te da estadisticas de cuantos SNPs missing por muestra y tambien de cuantos genotipos missing por SNPs de esto ultimo tienes que hacer un historgram en R para ver la distribucion
# dependiendo de como de mal o bien esté ponemos un filtro por el missing value antes de hacer la PCA

## PCA graph ##
# perform linkage pruning - i.e. identify prune sites (not necessary for already created bed files)
plink --bfile ${BFILE}_filt --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --pca --out ${BFILE}

# Distance square for a heatmap
plink --bfile ${BFILE}_filt --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --distance square --out ${BFILE}

# Excluding sex chromosomes
plink --bfile ${BFILE}_filt --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --not-chr NC_006127.5 NC_006126.5 --distance square --out ${BFILE}

# allele frequencies
plink --bfile ${BFILE}_filt --double-id --keep-allele-order --allow-extra-chr --real-ref-alleles --freq --out SNPs_cobb_filt

# heterozygosity
 plink --bfile ${BFILE}_filt --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --het --out SNPs_cobb_filt # add keep if we need to add samples

 # non binary file > PED file
 plink --bfile ${BFILE}_filt --recode --out ${BFILE}

# Convert again to vcf
plink --bfile ${BFILE}_filt --double-id --allow-extra-chr --keep-allele-order  --real-ref-alleles --recode vcf-iid bgz --out ${BFILE}_filt

# QWAS analysis
plink --bfile ${BFILE}_filt --double-id --allow-extra-chr --keep-allele-order --assoc --out chicken_SNPs_PCA_out_prun
plink --bfile IC_all_SNPs --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --extract 17K.snps --recode vcf --out IC_17K_SNPs_plink

# Keep a number of SNPs
plink --bfile ${BFILE}_filt --thin-count 1393289 --recode vcf-iid bgz --out HD_4M
################################################################################

# Fst
module load perl/5.24.0 vcftools/0.1.16
#Fst analysis
vcftools --gzvcf ${VCF} --weir-fst-pop Cobbs.txt --weir-fst-pop Rosses.txt --fst-window-size interger --out cobb_vs_ross
# filtering with VCFtools
vcftools --gzvcf ${VCF} --maf 0.1 --max-missing 1 --recode recode-INFO-all --out HD_filt_snps.vcf
## some quality check tests:
# --maf --max-maf include only sites with minor allele frequency greated than or equal to the maf value and less than or equal to the max-maf.
# allele freq: the number of times an allele appears over all individuals at that site, divided by the total number of non-missing alleles at that time.
# --max-missing-count: exclude sites with more than this number of missing genotypes over all individuals.

# extracting SNPs
vcftools --gzvcf ${VCF} --positions 17K_vcf.snps --recode --out IC_17K_SNPs.vcf
cut -f2 bim
plink --bfile IC_all_SNPs --double-id --allow-extra-chr --keep-allele-order --real-ref-alleles --extract 17K.snps --recode vcf --out IC_17K_SNPs_plink

# extract sample from the whole file
bcftools view -s CA19_18C1a -Oz -o CA19_18C1a_IC__10K_SNPs.vcf.gz IC_17K_SNPs.vcf.recode.vcf




# Linkage disequilibirum decay
#genotype correlation with plink
#(there is a more accurate way to do it with phased haplotypes in vcftools, but we can do it later)
#needs to be calculated by breed. You can calculate it for the largest chromosome only. --ld-window-kb need to be adjusted according to your species (for honeybee 10kb is way sufficient). Roughly I think it could be an order of magnitude larger, so start with 100 for chicken.
plink --bfile chickens --keep cobbindividuals.list --chr [name of chr] --maf 0.1 --r2 dprime --ld-window-kb 10 --ld-window 99999 --ld-window-r2 0 --out LD_breed1
plink --bfile all_SNPs_filt --keep ross_filt.txt --chr NC_006088.5,NC_006089.5,NC_006090.5,NC_006091.5,NC_006092.5,NC_006093.5,NC_006094.5,NC_006095.5,NC_006096.5,NC_006097.5 --r2 dprime --ld-window-kb 100 --ld-window 99999 --ld-window-r2 0 --allow-extra-chr --keep-allele-order --real-ref-alleles --out macro_LD_ross


# SNP annotation
SNPeff
SNPeff='/home/people/sofbas/CustomSoftware/snpEff'
cd $SNPeff
snpEff build -gff3 -v Amel_HAv3.1
/home/people/sofbas/CustomSoftware/snpEff/data

# Mouse genome, version mm37.61
mm37.61.genome : Mouse
# Chicken genome, version 6
GRCg6a.genome : Chicken GRCg6a (Galgal6)

/home/people/sofbas/CustomSoftware/snpEff/snpEff.jar
module load java/1.8.0
workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/LD_imp/all-chroms/05_PopGen'
java -Xmx8g -jar /home/people/sofbas/CustomSoftware/snpEff/snpEff.jar build -gtf22 -v GRCg6a

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz
mv GCF_000002315.6_GRCg6a_genomic.gff.gz genes.gff.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.fna.gz
mv GCF_000002315.6_GRCg6a_genomic.fna.gz GRCg6a.fa.gz

cd ${workdir}/SnpEff/
java -Xmx8g -jar /home/people/sofbas/CustomSoftware/snpEff/snpEff.jar -v GRCg6a.99 ${workdir}/Fst/ALL_SNPs_filt.vcf > all_annotated.vcf
gzip ${workdir}/SnpEff/all_annotated.vcf


# pFst
module load vcflib/1.0.0-rc2
# https://github.com/vcflib/vcflib
### https://github.com/vcflib/vcflib/blob/master/doc/pFst.md
# Esto genera un archivo mucho más grande. Primero hacer por windows, después pensamos si hacer esto.
# fist we needed to filter with vcftools instead of using plink because we need gp regions.
module load perl/5.24.0 vcftools/0.1.16
vcftools --gzvcf ALL_SNPs.vcf.gz --maf 0.1 --max-missing 1 --recode --recode-INFO-all --out all_filt_withvcf.vcf
pFst --target 53,55,56,57,60,64,65,69,71,72,79,82,83,84,85,87,89,90,91,95,101,102,103,104,105,106,108,117,118,119,120,121,122,124,125,127 --background 54,58,59,61,62,63,66,67,68,70,72,73,74,76,77,78,80,81,86,88,92,93,94,96,97,98,99,100,107,109,110,111,112,113,114,115,116,123,126,128,129,130 --file all_filt_withvcf.vcf --deltaaf 0.1 --type GP > test.out


## neighbor-joining (NJ) trees
module load rapidnj/2.3.2
rapidnj -i
(pd = distance matrix in phylip format)
sth = multiple alignment in (single line) stockholm format.
