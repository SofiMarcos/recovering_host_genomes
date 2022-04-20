# Damage profiler
# https://damageprofiler.readthedocs.io/en/latest/

bam="/home/projects/ku-cbd/data/HoloFood/ChickenGutContent/MAG_caecum/03_MapToRef/"
ref="/home/projects/ku-cbd/data/HoloFood/Chicken_REF_Genome/GCF_000002315.6_GRCg6a_genomic.fna"

java -jar /home/people/sofbas/CustomSoftware/DamageProfiler/DamageProfiler-1.1-java11.jar -i ${bam}CB15_17F1a_map2host.bam -o CB15_17 -r ${ref} -title CB15_17
java -jar /home/people/sofbas/CustomSoftware/DamageProfiler/DamageProfiler-1.1-java11.jar -i ${bam}CB10_16F1a_map2host.bam -o CB10_16 -r ${ref} -title CB10_16
java -jar /home/people/sofbas/CustomSoftware/DamageProfiler/DamageProfiler-1.1-java11.jar -i ${bam}CC14_08F1a_map2host.bam -o CC14_08 -r ${ref} -title CC14_08
java -jar /home/people/sofbas/CustomSoftware/DamageProfiler/DamageProfiler-1.1-java11.jar -i ${bam}CA22_12F1a_map2host.bam -o CA22_12 -r ${ref} -title CA22_12
java -jar /home/people/sofbas/CustomSoftware/DamageProfiler/DamageProfiler-1.1-java11.jar -i ${bam}CA08_08F1a_map2host.bam -o CA08_08 -r ${ref} -title CA08_08
java -jar /home/people/sofbas/CustomSoftware/DamageProfiler/DamageProfiler-1.1-java11.jar -i ${bam}CA06_05F1a_map2host.bam -o CA06_05 -r ${ref} -title CA06_05
java -jar /home/people/sofbas/CustomSoftware/DamageProfiler/DamageProfiler-1.1-java11.jar -i ${bam}CC19_18F1a_map2host.bam -o CC19_18 -r ${ref} -title CC19_18
java -jar /home/people/sofbas/CustomSoftware/DamageProfiler/DamageProfiler-1.1-java11.jar -i ${bam}CC11_07F1a_map2host.bam -o CC11_07 -r ${ref} -title CC11_07
java -jar /home/people/sofbas/CustomSoftware/DamageProfiler/DamageProfiler-1.1-java11.jar -i ${bam}CB08_08F1a_map2host.bam -o CB08_08 -r ${ref} -title CB08_08
java -jar /home/people/sofbas/CustomSoftware/DamageProfiler/DamageProfiler-1.1-java11.jar -i ${bam}CB01_15F1a_map2host.bam -o CB01_15 -r ${ref} -title CB01_15

#samples:
CB10_16F1a_map2host.bam
CC14_08F1a_map2host.bam
CA22_12F1a_map2host.bam
CA08_08F1a_map2host.bam
CC20_14F1a_map2host.bam
CA17_02F1a_map2host.bam
CA19_18F1a_map2host.bam
CB22_10F1a_map2host.bam

CC19_18F1a_map2host.bam
CC11_07F1a_map2host.bam
CB08_08F1a_map2host.bam
CB01_15F1a_map2host.bam
CA03_14F1a_map2host.bam
CA08_16F1a_map2host.bam
CA06_05F1a_map2host.bam
