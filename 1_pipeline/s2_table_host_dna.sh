# Examples of previous datasets

# Cattle
# 1-
#!/bin/bash
module load anaconda3/4.4.0 enabrowsertools/1.6
Cattle='/home/projects/ku-cbd/data/HoloFood/Sofi/Host/cattle'
cd ${Cattle}/
for SAMPLE in ERR3201375 ERR3201376	ERR3201377 ERR3201378 ERR3201379 ERR3201380	ERR3201381 ERR3201382 ERR3201383 ERR3201384
do
  enaDataGet -f fastq -d ${Cattle}/ ${SAMPLE}
  rm ${Cattle}/${SAMPLE}/${SAMPLE}.fastq.gz
done
#ena='/home/projects/ku-cbd/data/HoloFood/Sofi/Host'
#qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $ena/ena.err -o $ena/ena.out -l nodes=1:ppn=40,mem=180gb,walltime=0:10:00:00 -N ena_cattle $ena/ena.sh #3 hours for 10 samples
mv ERR*/*fastq.gz 00_RawData/
rm -r ERR*

# 2-
#!/bin/bash
module load bwa/0.7.16a
cd /home/projects/ku-cbd/data/HoloFood/Sofi/Host/ref_genomes/bos_taurus
bwa index GCF_002263795.1_ARS-UCD1.2_genomic.fna
#workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/Host/ref_genomes/bos_taurus'
#qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $workdir/bos.err -o $workdir/bos.out -l nodes=1:ppn=40,mem=180gb,walltime=0:50:00:00 -N ref_cattle $workdir/bos.sh #3 hours for 10 samples

# 3-
#!/bin/bash
module load samtools/1.11
module load bwa/0.7.16a
workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/Host/cattle'
ref='/home/projects/ku-cbd/data/HoloFood/Sofi/Host/ref_genomes/bos_taurus/GCF_002263795.1_ARS-UCD1.2_genomic.fna'
for sample in ERR3201375 ERR3201376	ERR3201377 ERR3201378 ERR3201379 ERR3201380	ERR3201381 ERR3201382 ERR3201383 ERR3201384
do
  bwa mem -t 40 -R "@RG\\tID:HoloFood\\tCN:SenLi\\tDS:Mappingt\\tPL:Illumina1.9\\tSM:SM" ${ref} ${workdir}/00_RawData/${sample}_1.fastq.gz ${workdir}/00_RawData/${sample}_2.fastq.gz | samtools sort --threads 40 -o ${workdir}/01_Mapped/${sample}.bam
done
#workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/Host/cattle/01_Mapped'
#qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $workdir/bos_map.err -o $workdir/bos_map.out -l nodes=1:ppn=40,mem=180gb,walltime=0:50:00:00 -N map_cattle $workdir/map.sh #3 hours for 10 samples

# 4-
#!/bin/bash
module load samtools/1.11
workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/Host/cattle'
for sample in ERR3201375 ERR3201376	ERR3201377 ERR3201378 ERR3201379 ERR3201380	ERR3201381 ERR3201382 ERR3201383 ERR3201384
do
  samtools flagstat ${workdir}/01_Mapped/${sample}.bam > ${workdir}/02_Stats/${sample}.stats
  samtools depth ${workdir}/01_Mapped/${sample}.bam > ${workdir}/02_Stats/${sample}.depth
done
#workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/Host'
#qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $workdir/stats.err -o $workdir/stats.out -l nodes=1:ppn=40,mem=180gb,walltime=0:50:00:00 -N stats $workdir/stats.sh #3 hours for 10 samples

# 5- Depth
for sample in ERR3201375 ERR3201376	ERR3201377 ERR3201378 ERR3201379 ERR3201380	ERR3201381 ERR3201382 ERR3201383 ERR3201384
do
  less ${sample}.depth | awk '{c++;s+=$3}END{print s/c}'
done

# Pig
# 0-
#!/bin/bash
workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/Host'
 tar -zcvf ${workdir}/cattle_10_samples.tar.gz ${workdir}/cattle/
#workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/Host'
#qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $workdir/com.err -o $workdir/com.out -l nodes=1:ppn=40,mem=180gb,walltime=0:05:00:00 -N compress $workdir/compreess.sh

# 0.1-
#!/bin/bash
cd /home/projects/ku-cbd/data/HoloFood/Sofi/Host/pig/00_RawData/
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178547/CNX0130695/CNR0162624/FCR.HF1_300_1.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178547/CNX0130695/CNR0162624/FCR.HF1_300_2.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178548/CNX0130696/CNR0162625/FCR.HF2_300_1.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178548/CNX0130696/CNR0162625/FCR.HF2_300_2.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178549/CNX0130697/CNR0162626/FCR.HF3_300_1.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178549/CNX0130697/CNR0162626/FCR.HF3_300_2.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178550/CNX0130698/CNR0162627/FCR.HF4_300_1.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178550/CNX0130698/CNR0162627/FCR.HF4_300_2.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178551/CNX0130699/CNR0162628/FCR.HF5_300_1.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178551/CNX0130699/CNR0162628/FCR.HF5_300_2.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178552/CNX0130700/CNR0162629/FCR.HM1_300_1.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178552/CNX0130700/CNR0162629/FCR.HM1_300_2.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178553/CNX0130701/CNR0162630/FCR.HM2_300_1.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178553/CNX0130701/CNR0162630/FCR.HM2_300_2.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178554/CNX0130702/CNR0162631/FCR.HM3_300_1.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178554/CNX0130702/CNR0162631/FCR.HM3_300_2.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178555/CNX0130703/CNR0162632/FCR.HM4_300_1.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178555/CNX0130703/CNR0162632/FCR.HM4_300_2.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178556/CNX0130704/CNR0162633/FCR.HM5_300_1.fq.gz
wget https://ftp.cngb.org/pub/CNSA/data2/CNP0000824/CNS0178556/CNX0130704/CNR0162633/FCR.HM5_300_2.fq.gz
#qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $workdir/pig.err -o $workdir/pig.out -l nodes=1:ppn=40,mem=180gb,walltime=0:05:00:00 -N pig_data2 $workdir/download_pig.sh

# 1-
#!/bin/bash
module load bwa/0.7.16a
cd /home/projects/ku-cbd/data/HoloFood/Sofi/Host/ref_genomes/sus_scrofa
bwa index GCF_000003025.6_Sscrofa11.1_genomic.fna
#workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/Host'
#qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $workdir/sus.err -o $workdir/sus.out -l nodes=1:ppn=40,mem=180gb,walltime=0:10:00:00 -N ref_pig $workdir/sus.sh

# 2-
#!/bin/bash
module load samtools/1.11 bwa/0.7.16a
workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/Host/pig'
ref='/home/projects/ku-cbd/data/HoloFood/Sofi/Host/ref_genomes/sus_scrofa/GCF_000003025.6_Sscrofa11.1_genomic.fna'
for sample in FCR.HF1_300 FCR.HF2_300 FCR.HF3_300 FCR.HF4_300 FCR.HF5_300 FCR.HM1_300 FCR.HM2_300 FCR.HM3_300 FCR.HM4_300 FCR.HM5_300
do
  bwa mem -t 40 -R "@RG\\tID:HoloFood\\tCN:SenLi\\tDS:Mappingt\\tPL:Illumina1.9\\tSM:SM" ${ref} ${workdir}/00_RawData/${sample}_1.fastq.gz ${workdir}/00_RawData/${sample}_2.fastq.gz | samtools sort --threads 40 -o ${workdir}/01_Mapped/${sample}.bam
  samtools flagstat ${workdir}/01_Mapped/${sample}.bam > ${workdir}/02_Stats/${sample}.stats
  samtools depth -a ${workdir}/01_Mapped/${sample}.bam > ${workdir}/02_Stats/${sample}.depth
done
#workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/Host'
#qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $workdir/sus_map.err -o $workdir/sus_map.out -l nodes=1:ppn=40,mem=180gb,walltime=0:50:00:00 -N map_pig $workdir/sus_map.sh

# 5- Depth
#!/bin/bash
cd /home/projects/ku-cbd/data/HoloFood/Sofi/Host/pig/02_Stats/
for sample in FCR.HF1_300 FCR.HF2_300 FCR.HF3_300 FCR.HF4_300 FCR.HF5_300 FCR.HM1_300 FCR.HM2_300 FCR.HM3_300 FCR.HM4_300 FCR.HM5_300
do
  less ${sample}.depth | awk '{c++;s+=$3}END{print s/c}'
done
#workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/Host'
#qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $workdir/stats_pig.err -o $workdir/stats_pig.out -l nodes=1:ppn=40,mem=180gb,walltime=0:05:00:00 -N stats_pig $workdir/stats.sh


# Human
# 0-
#!/bin/bash
module load anaconda3/4.4.0 enabrowsertools/1.6
Human='/home/projects/ku-cbd/data/HoloFood/Sofi/Host/human_1'
cd ${Human}/
for SAMPLE in ERR4327984 ERR4327985 ERR4327986 ERR4327987 ERR4327988 ERR4327989 ERR4327990 ERR4327991	ERR4327992 ERR4327993
do
  enaDataGet -f fastq -d ${Human}/ ${SAMPLE}
  rm ${Human}/${SAMPLE}/${SAMPLE}.fastq.gz
done
ena='/home/projects/ku-cbd/data/HoloFood/Sofi/Host'
qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $ena/down_human.err -o $ena/down_human.out -l nodes=1:ppn=40,mem=180gb,walltime=0:10:00:00 -N down_human $ena/down_human.sh #3 hours for 10 samples
mv ERR*/*fastq.gz 00_RawData/
rm -r ERR*

# 1-
#!/bin/bash
module load bwa/0.7.16a
cd /home/projects/ku-cbd/data/HoloFood/Sofi/Host/ref_genomes/homo_sapiens
bwa index GCF_000001405.39_GRCh38.p13_genomic.fna
workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/Host'
qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $workdir/homo.err -o $workdir/homo.out -l nodes=1:ppn=40,mem=180gb,walltime=0:10:00:00 -N ref_human $workdir/homo_ref.sh


# 2-
#!/bin/bash
module load samtools/1.11 bwa/0.7.16a
workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/Host/human_1'
ref='/home/projects/ku-cbd/data/HoloFood/Sofi/Host/ref_genomes/homo_sapiens/GCF_000001405.39_GRCh38.p13_genomic.fna'
for sample in ERR4327984 ERR4327985 ERR4327986 ERR4327987 ERR4327988 ERR4327989 ERR4327990 ERR4327991	ERR4327992 ERR4327993
do
  bwa mem -t 40 -R "@RG\\tID:HoloFood\\tCN:SenLi\\tDS:Mappingt\\tPL:Illumina1.9\\tSM:SM" ${ref} ${workdir}/00_RawData/${sample}_1.fastq.gz ${workdir}/00_RawData/${sample}_2.fastq.gz | samtools sort --threads 40 -o ${workdir}/01_Mapped/${sample}.bam
  samtools flagstat ${workdir}/01_Mapped/${sample}.bam > ${workdir}/02_Stats/${sample}.stats
  samtools depth -a ${workdir}/01_Mapped/${sample}.bam > ${workdir}/02_Stats/${sample}.depth
done
workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/Host'
qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $workdir/map_human.err -o $workdir/map_human.out -l nodes=1:ppn=40,mem=180gb,walltime=0:50:00:00 -N map_human $workdir/map_human.sh

# 3-
#!/bin/bash
cd /home/projects/ku-cbd/data/HoloFood/Sofi/Host/human_1/02_Stats/
for sample in ERR4327984 ERR4327985 ERR4327986 ERR4327987 ERR4327988 ERR4327989 ERR4327990 ERR4327991	ERR4327992 ERR4327993
do
  less ${sample}.depth | awk '{c++;s+=$3}END{print s/c}'
done
workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/Host'
qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $workdir/stats_human.err -o $workdir/stats_human.out -l nodes=1:ppn=40,mem=180gb,walltime=0:05:00:00 -N stats_human $workdir/stats_human.sh


qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $workdir/stats_2.err -o $workdir/stats_2.out -l nodes=1:ppn=40,mem=180gb,walltime=0:15:00:00 -N stats_2 $workdir/stats_2.sh


# Rat
# 1-
#!/bin/bash
module load anaconda3/4.4.0 enabrowsertools/1.6
Rat='/home/projects/ku-cbd/data/HoloFood/Sofi/Host/rat'
cd ${Rat}/
for SAMPLE in SRR11223902	SRR11223903	SRR11223904	SRR11223905	SRR11223906	SRR11223907	SRR11223908	SRR11223909	SRR11223910	SRR11223911
do
  enaDataGet -f fastq -d ${Rat}/ ${SAMPLE}
  rm ${Rat}/${SAMPLE}/${SAMPLE}.fastq.gz
done

Fish='/home/projects/ku-cbd/data/HoloFood/Sofi/Host/fish'
cd ${Fish}/
for SAMPLE in SRR10827562 SRR10827563 SRR10827564	SRR10827565 SRR10827566	SRR10827567	SRR10827568	SRR10827569	SRR10827570	SRR10827571
do
  enaDataGet -f fastq -d ${Fish}/ ${SAMPLE}
  rm ${Fish}/${SAMPLE}/${SAMPLE}.fastq.gz
done

module load bwa/0.7.16a
cd /home/projects/ku-cbd/data/HoloFood/Sofi/Host/ref_genomes/danio_rerio
bwa index GCF_000002035.6_GRCz11_genomic.fna

cd /home/projects/ku-cbd/data/HoloFood/Sofi/Host/ref_genomes/rattus_rattus
bwa index GCF_011064425.1_Rrattus_CSIRO_v1_genomic.fna

ena='/home/projects/ku-cbd/data/HoloFood/Sofi/Host'
qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $ena/prepare.err -o $ena/prepare.out -l nodes=1:ppn=40,mem=180gb,walltime=0:10:00:00 -N prepare $ena/prepare.sh #3 hours for 10 samples
mv ERR*/*fastq.gz 00_RawData/
rm -r ERR*


# Zebrafish
# 1- 
#!/bin/bash
module load samtools/1.11 bwa/0.7.16a
workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/Host/rat/'
ref='/home/projects/ku-cbd/data/HoloFood/Sofi/Host/ref_genomes/rattus_rattus/GCF_011064425.1_Rrattus_CSIRO_v1_genomic.fna'
for sample in SRR11223902	SRR11223903	SRR11223904	SRR11223905	SRR11223906	SRR11223907	SRR11223908	SRR11223909	SRR11223910	SRR11223911
do
  bwa mem -t 40 -R "@RG\\tID:HoloFood\\tCN:SenLi\\tDS:Mappingt\\tPL:Illumina1.9\\tSM:SM" ${ref} ${workdir}/00_RawData/${sample}_1.fastq.gz ${workdir}/00_RawData/${sample}_2.fastq.gz | samtools sort --threads 40 -o ${workdir}/01_Mapped/${sample}.bam
  samtools flagstat ${workdir}/01_Mapped/${sample}.bam > ${workdir}/02_Stats/${sample}.stats
done

workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/Host/fish/'
ref='/home/projects/ku-cbd/data/HoloFood/Sofi/Host/ref_genomes/danio_rerio/GCF_000002035.6_GRCz11_genomic.fna'
for sample in SRR10827562 SRR10827563 SRR10827564	SRR10827565 SRR10827566	SRR10827567	SRR10827568	SRR10827569	SRR10827570	SRR10827571
do
  bwa mem -t 40 -R "@RG\\tID:HoloFood\\tCN:SenLi\\tDS:Mappingt\\tPL:Illumina1.9\\tSM:SM" ${ref} ${workdir}/00_RawData/${sample}_1.fastq.gz ${workdir}/00_RawData/${sample}_2.fastq.gz | samtools sort --threads 40 -o ${workdir}/01_Mapped/${sample}.bam
  samtools flagstat ${workdir}/01_Mapped/${sample}.bam > ${workdir}/02_Stats/${sample}.stats
done

workdir='/home/projects/ku-cbd/data/HoloFood/Sofi/Host'
qsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e $workdir/map.err -o $workdir/map.out -l nodes=1:ppn=40,mem=180gb,walltime=0:50:00:00 -N map_ratfish $workdir/map.sh
