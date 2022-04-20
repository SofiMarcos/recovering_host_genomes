#!/bin/bash
####
# Preprocessing steps for metagenomic data
####

# directories
mkdir ${workdir}/PrePro/${dataset}
mkdir ${workdir}/PrePro/${dataset}/{00_RawData,01_QualityFiltered,02_DupRemoved,03_MapToRef,04_BamQC,05_MetaGfqs}

# paths
workdir='path/to/directory'
dataset='type_of_dataset' # Target population or Internal reference samples
ref='path/to/reference'
sample='CA00_00' # sample codes

# parameters
adapter1="adapter1_sequence"
adapter2="adapter2_sequence"
threads="40"

# 01 - quality filtering
module load gcc/8.2.0 AdapterRemoval/2.2.4
AdapterRemoval --file1 ${workdir}/PrePro/00_RawData/${sample}_1.fq.gz --file2 ${workdir}/PrePro/00_RawData/${sample}_2.fq.gz --settings ${workdir}/PrePro/01_QualityFiltered/${sample}.settings --singleton ${workdir}/PrePro/01_QualityFiltered/${sample}.singletons --discarded ${workdir}/PrePro/01_QualityFiltered/${sample}.discarded --output1 ${workdir}/PrePro/01_QualityFiltered/${sample}_filtered_1.fq --output2 ${workdir}/PrePro/01_QualityFiltered/${sample}_filtered_2.fq --trimqualities --trimns --maxns 5 --minquality 30 --minlength 100 --threads ${threads} --adapter1 ${adapter1} --adapter2 ${adapter2}

# 02 - duplicates removal
rm ${workdir}/PrePro/01_QualityFiltered/${sample}.settings ${workdir}/PrePro/01_QualityFiltered/${sample}.singletons ${workdir}/PrePro/01_QualityFiltered/${sample}.discarded
module load pigz/2.3.4 seqkit/0.7.1
paste -d '^' ${workdir}/PrePro/01_QualityFiltered/${sample}_filtered_1.fq ${workdir}/PrePro/01_QualityFiltered/${sample}_filtered_2.fq | seqkit rmdup -s -j ${threads} -o ${workdir}/PrePro/02_DupRemoved/${sample}.fq.tmp
< ${workdir}/PrePro/02_DupRemoved/${sample}.fq.tmp tee >(cut -d '^' -f 2 > ${workdir}/PrePro/02_DupRemoved/${sample}_filtered_rmdup_2.fq) | cut -d '^' -f 1 > ${workdir}/PrePro/02_DupRemoved/${sample}_filtered_rmdup_1.fq

# 03 - alignment
rm ${workdir}/PrePro/02_DupRemoved/${sample}.fq.tmp
module load samtools/1.9 bwa/0.7.16a
bwa mem -t ${threads} -R "@RG\\tID:tID\\tCN:tCN\\tDS:tDS\\tPL:tPL\\tSM:${sample}" ${ref} ${workdir}/PrePro/02_DupRemoved/${sample}_filtered_rmdup_1.fq ${workdir}/PrePro/02_DupRemoved/${sample}_filtered_rmdup_2.fq | samtools sort --threads ${threads} -o ${workdir}/PrePro/03_MapToRef/${sample}_all.bam
samtools view -b -F12 ${workdir}/PrePro/03_MapToRef/${sample}_all.bam > ${workdir}/PrePro/03_MapToRef/${sample}_map2host.bam
samtools view -b -f12 ${workdir}/PrePro/03_MapToRef/${sample}_all.bam > ${workdir}/PrePro/03_MapToRef/${sample}_metaG.bam
samtools index ${workdir}/PrePro/03_MapToRef/${sample}_map2host.bam

# 04 - quiality control
module load intel/perflibs/2019_update5 gcc/8.2.0 java/1.8.0-openjdk R/3.6.1 qualimap/2.2.1
qualimap bamqc -bam ${workdir}/PrePro/03_MapToRef/${sample}_all.bam -outdir ${workdir}/PrePro/04_BamQC/${sample}_all -outformat html --java-mem-size=16G
qualimap bamqc -bam ${workdir}/PrePro/03_MapToRef/${sample}_map2host.bam -outdir ${workdir}/PrePro/04_BamQC/${sample}_map2host -outformat html --java-mem-size=16G
mv ${workdir}/PrePro/03_MapToRef/${sample}_metaG.bam ${workdir}/PrePro/04_BamQC/${sample}_metaG.bam

# 05 - exporting metaG
module load samtools/1.11
samtools fastq --threads ${threads} -1 ${workdir}/PrePro/05_MetaGfqs/${sample}_metaG_1.fq.gz -2 ${workdir}/PrePro/05_MetaGfqs/${sample}_metaG_2.fq.gz ${workdir}/PrePro/04_BamQC/${sample}_metaG.bam
