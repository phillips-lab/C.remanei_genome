#!/bin/bash

##################################################################
#########Generates input files for the futher Hi-C analysis#######
##################################################################


#module load python2/2.7.15 samtools bwa bedtools java

ref="/path-to/ref/CR.ncbi.fasta"
#samtools="/path-to/samtools"
#bedtools="/path-to/bedtools"
#bwa="/path-to/bwa"
#python2.7
picard="/path-to/scripts/picard/build/libs/picard.jar"
FILTER="/path-to/ARIMA/mapping_pipeline/filter_five_end.pl"
COMBINER="/path-to/ARIMA/mapping_pipeline/two_read_bam_combiner.pl"
STATS="/path-to/ARIMA/mapping_pipeline/get_stats.pl"
tmp="/path-to/tmp"


#C. elegans data from SRR5341677
#our C. remanei Hi-C data
LISTFILES=(*1.fastq)

R1=${LISTFILES[$SLURM_ARRAY_TASK_ID]}

#name=${R1/.filt-trimmed-pair1.fastq/}
name=${R1/_1.fastq/}


#######################################################################################################
############### next part is ARIMA mapping pipeline from https://github.com/ArimaGenomics/mapping_pipeline

#### Step 1.A: FASTQ to BAM (1st) #https://github.com/ArimaGenomics/mapping_pipeline
bwa index $ref
bwa mem -M -t 20 $ref $R1 | samtools view -@ 20 -Sb - > ${name}_1.bam
bwa mem -M -t 20 $ref ${R1/_1/_2} | samtools view -@ 20 -Sb - > ${name}_2.bam

#### Step 2.A: Filter 5' end (1st) #https://github.com/ArimaGenomics/mapping_pipeline
samtools view -h ${name}_1.bam | perl $FILTER | samtools view -Sb - > ${name}_1.filt5.bam
samtools view -h ${name}_2.bam | perl $FILTER | samtools view -Sb - > ${name}_2.filt5.bam

#### Step 3A: Pair reads & mapping quality filter #https://github.com/ArimaGenomics/mapping_pipeline
perl $COMBINER ${name}_1.filt5.bam ${name}_2.filt5.bam samtools 30 | samtools view -bS -t $ref.fai - | samtools sort -@ 20 -o ${name}.TEMP.bam -

#### Step 4: Mark duplicates #https://github.com/ArimaGenomics/mapping_pipeline

java -Xmx95g -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$tmp -jar $picard MarkDuplicates INPUT=${name}.TEMP.bam OUTPUT=${name}.rep.bam METRICS_FILE=${name}.rep.metrics.txt REMOVE_DUPLICATES=true ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT TMP_DIR=$tmp;

samtools index ${name}.rep.bam
#perl $STATS ${name}.rep.bam > ${name}.rep.bam.stats #https://github.com/ArimaGenomics/mapping_pipeline


bamToBed -i ${name}.rep.bam > ${name}.rep.bed
sort -k 4 ${name}.rep.bed >  ${name}.rep.s.bed

################## the end of ARIMA mapping pipeline ########
#############################################################




#get odd and even lines, then join the files

awk 'NR % 2 == 1'  ${name}.rep.s.bed >${name}.ODD
awk 'NR % 2 == 0'  ${name}.rep.s.bed >${name}.EVEN
paste ${name}.ODD ${name}.EVEN > ${name}.PAIRED

awk '{if($1==$7){print $1"\t"$2-1"\t"$2"\t"$8-$2}}' ${name}.PAIRED > ${name}.same


#replace the chromosome names with numbers
sed -e 's/MtDNA/7/g' ${name}.PAIRED> ${name}.PAIRED.names
sed -i 's/III/3/g' ${name}.PAIRED.names
sed -i 's/II/2/g' ${name}.PAIRED.names
sed -i 's/IV/4/g' ${name}.PAIRED.names
sed -i 's/V/5/g' ${name}.PAIRED.names
sed -i 's/I/1/g' ${name}.PAIRED.names
sed -i 's/X/6/g' ${name}.PAIRED.names

#links on different chromosomes
#get count per 100Kb windows
awk '{if($1!=$7){print $1"\t"int($2/100000)"\t"$7"\t"int($8/100000)}}' ${name}.PAIRED.names | awk '{if($3<$1){print $3"\t"$4"\t"$1"\t"$2;} else {print;}}' - |sort -k1 -k2n -| uniq -c > ${name}.differ



###########################################################################
### the median distance between paired reads mapped to the same chromosome

sed -i 's/-//g' ${name}.same
grep -v "MtDNA" ${name}.same | sort -k 1,1 -k2,2n - | awk '{if($3>$2)print;}' - > ${name}.s.same

#clean up
rm ${name}.same ${name}.PAIRED ${name}.PAIRED.names ${name}.ODD ${name}.EVEN

#100 kb windows
samtools faidx ref.fasta
bedtools makewindows -g ref.fasta.fai -w 100000 > ref_100kb_winds.bed
bedtools map -b ${name}.s.same -a ref_100kb_winds.bed -c 4 -o median > ${name}_100.bed
