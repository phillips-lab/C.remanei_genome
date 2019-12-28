#!/bin/bash

###Script to generate fractions in C.remanei
###

ref="/path-to/ref/CR.ncbi.fasta"
refmask="/path-to/ref/CR.ncbi.hardmasked.fasta"
ann="/path-to/ref/CR.ncbi.gff"
sp="CR"

### change these aprameters for C. elegans
####C.elegans data:
##C.elegans reference genome https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB28388
##C.elegans annotations from Supplemental_Data_File_S6 https://genome.cshlp.org/content/29/6/1009/suppl/DC1
##bed files with masked intervals are in Supplimentary files S3 and S4
##(https://figshare.com/projects/Supplementary_Tables_and_Figures_of_Chromosome-level_assembly_of_the_Caenorhabditis_remanei_genome_reveals_conserved_patterns_of_nematode_genome_organization_/73518)


####estimate fractions of exonic, intronic DNA and other features.

samtools faidx $ref
bedtools makewindows -g ${ref}.fai -w 100000 > ${sp}_100kb_winds.bed

##estimate exon/intron/gene coverage
grep "exon" $ann |awk '{print$1"\t"$4-1"\t"$5-1}' -  > ${sp}.exon.BED
grep "intron" $ann |awk '{print$1"\t"$4-1"\t"$5-1}' -  > ${sp}.intron.BED
grep "gene" $ann |awk '{print$1"\t"$4-1"\t"$5-1}' -  > ${sp}.gene.BED

bedtools coverage -a ${sp}_100kb_winds.bed -b ${sp}.exon.BED >${sp}_exon_cov.bed
bedtools coverage -a ${sp}_100kb_winds.bed -b ${sp}.intron.BED >${sp}_intron_cov.bed
bedtools coverage -a ${sp}_100kb_winds.bed -b ${sp}.gene.BED >${sp}_gene_cov.bed


#estimate repatitive fractions
bedtools nuc -bed {sp}_100kb_winds.bed -fi $refmask > ${sp}_N_content.bed

#estimate GC content
bedtools nuc -bed {sp}_100kb_winds.bed -fi $ref > ${sp}_GC_content.bed

#overlap between introns and reperirive, script generate_masked_ranges.py is available on https://gist.github.com/danielecook/cfaa5c359d99bcad3200
python2 generate_masked_ranges.py $refmaks > ${sp}.masked.bed
bedtools intersect -a ${sp}.masked.bed -b ${sp}.intron.BED > ${sp}.intron_repeat_overlap.BED

cat ${sp}.intron.BED | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' > ${sp}_INTRONS_TOTAL.txt
cat ${sp}.intron_repeat_overlap.BED | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' > ${sp}_REP_OVER.txt
