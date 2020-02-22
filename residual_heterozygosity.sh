#!/bin/bash

##############################################################################
#######Counts the residual heterozygosity in the C. remaeni PX506 line #######
##############################################################################

module load java bedtools bwa samtools



ref="/path-to/ref/CR.ncbi.fasta"
REPEATS="/path-to/MASKS/CR.repeats.regions.bed" #Supplementary file S3
tmp="/path-to/tmp"

picard="/path-to/scripts/picard/build/libs/picard.jar"



cd reads

module load easybuild  icc/2017.1.132-GCC-6.3.0-2.27  impi/2017.1.132 skewer


#filter nextera adapters

skewer -x CTGTCTCTTATA -t 12 -l 36 -r 0.01 -d 0.01 -q 20 -o PX506.oct.filt PX506_1.oct.fq PX506_2.oct.fq;
skewer -x CTGTCTCTTATA -t 12 -l 36 -r 0.01 -d 0.01 -q 20 -o PX506.filt PX506_1.fq PX506_2.fq;


mv PX506.filt-trimmed-pair1.fastq to PX506.fp1.fastq
mv PX506.filt-trimmed-pair2.fastq to PX506.fp2.fastq
mv PX506.oct.filt-trimmed-pair1.fastq to PX506.oct.fp1.fastq
mv PX506.oct.filt-trimmed-pair2.fastq to PX506.oct.fp2.fastq


#map reads
for file in *fp1.fastq;do

    name=${file/.fp1.fastq/}

    #map
    bwa mem -M -t 16 -R "@RG\tID:$name\tSM:$name\tPL:Illumina\tPI:330" $ref $file ${file/fp1./fp2.} > ${file/fp1.fastq/CR.sam} 2>${file/fp1.fastq/CR.bwa.log}

    #filter mapped reads
    samtools view -@ 15 -F 4 -bS -q 20 ${file/fp1.fastq/CR.sam} | samtools sort -@ 16 -o ${file/fp1.fastq/CR.s.bam} -

    samtools index ${file/fp1.fastq/CR.s.bam};

    #estimate the mean coverage
    samtools depth -a ${file/fp1.fastq/CR.s.bam} | awk '{sum+=$3} END { print "Average = ",sum/NR; }' > ${file/fp1.fastq/cov_mean};

    # remove duplicates
    java -Xmx4g -jar $picard MarkDuplicates INPUT=${file/fp1.fastq/CR.s.bam} OUTPUT=${file/fp1.fastq/ded.bam} METRICS_FILE=${file/.fp1.fastq/.ded.metrics.txt} REMOVE_DUPLICATES=true;
    samtools index ${file/fp1.fastq/ded.bam};

done

##rm *sam *CR.ba*

#merge

java -Xmx4g -jar $picard MergeSamFiles SORT_ORDER=coordinate ASSUME_SORTED=true \
    I=PX506.oct.ded.bam \
    I=PX506.ded.bam \
    O=PX506.M.bam

##rm *ded.ba*

#change the reading group

java -Xmx4g -jar $picard AddOrReplaceReadGroups I=PX506.M.bam O=PX506.re.bam RGPL=illumina RGLB=PX506 RGPU=NONE RGSM=PX506;

samtools index PX506.re.bam;

#realing with GATK 3.8
file=PX506.re.bam

java -Xmx5g -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -I $file -o ${file/.re.bam/.intervals} ;
java -Xmx5g -Djava.io.tmpdir=$tmp  -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -I $file -R $ref -targetIntervals ${file/.re.bam/.intervals} -o ${file/.re.bam/.ind.bam};


#call variants


java -Xmx6g -jar $EBROOTGATK/GenomeAnalysisTK.jar -R $ref -T HaplotypeCaller -I ${file/.re.bam/.ind.bam}  -o ${file/.re.bam/.raw.g.vcf} --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000;

java -Xmx30g -jar $EBROOTGATK/GenomeAnalysisTK.jar -R $ref -T GenotypeGVCFs -o PX506_resid_het.raw.vcf  -V PX506.M.raw.g.vcf;


#filter SNPs (biallelic, high quality, not near indels, not in repeats)
java -Xmx30g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T SelectVariants \
   -R $ref \
   -V PX506_resid_het.raw.vcf \
   -selectType SNP \
     -o PX506_resid_het.snp.vcf \
    --excludeIntervals $REPEATS \
      --excludeFiltered


java -Xmx30g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
    -T VariantFiltration \
   -R $ref \
   -V PX506_resid_het.snp.vcf \
                --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
                --filterName "hard_filter" \
                -o PX506_resid_het.filtsnp.vcf


grep -P "#|PASS" PX506_resid_het.filtsnp.vcf > PX506_resid_het.filtsnp.short.vcf


#count the number of polymorphic siter per 100 Kb windows

bedtools makewindows -g ${ref}.fai -w 100000 > winds.100kb.bed
bedtools map -b PX506_resid_het.filtsnp.short.vcf -a winds.100kb.bed -c 4 -o count > PX506_resid_het.counts.txt
