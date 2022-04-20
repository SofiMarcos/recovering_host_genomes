# Recovering high-quality host genomes from gut metagenomic data through genotype imputation

This repository contains the bioinfomatic resources related to the genomic analyses for retrieving host genomes from low‐coverage host data generated in metagenomic analysis.

## Raw data
Raw data arewill be available from European Nucleotide Archive (ENA), with BioProject accession no. PRJEB43192 (https://www.ebi.ac.uk/ena/browser/view/PRJEB43192?show = component‐projects). Until the release date, data will be made available upon request. ENA_sample codes of the used samples can be found in the ena_codes_host_bams.tsv.

## 1_pipeline - bioinformatic pipeline
All files for the R scripts can be found in 2_loocv directory.

**1_pre-processing** - divide host DNA sequences from metagenomic data.

**2_variant_calling-reference_samples** - design custom reference panels.

**3_two-step_imputation** - perform two-step imputation.

**4_popgen** - code for popgen analysis. 

**s1_damage_profiler** - code to see if read ends had undergone DNA damage.

**s2_table_host_dna** - download and calculate host DNA from previous works. 

## 2_loocv - Imputation accuracy using 12 validation samples
**1_loocv** - leave‐one‐out cross‐validation bash script.

**2_concordance** - calculate concordance for the 12 validation samples.

**3_t-test** - plots and t-tests for the 12 validation samples.

**4_missing_percentage** - plot of the missing percentage.

  

## 3_population_genetics_inference - R codes for ploting
All files for the R scripts can be found in the data directory.

**1_allele_freq** - allele frequencies distribution.

**2_heterozygosity** - heterozygosity of the imputed population and the 12 validation samples.

**3_nucleotide_diversity** - nucleotide diversity of the imputed population and the 12 validation samples.

**4_ibs_heatmap** - identity by state of the imputed population and the 12 validation samples.

**5_kinship** - kinship of the imputed population and the 12 validation samples.

**6_correlogram** - correlation plots of ibs and kinship.

**7_fst-comparison** - Venn diagrams for fst.


