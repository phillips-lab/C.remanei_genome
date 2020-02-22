
######################################################################
############### RNA-seq data, RPKM calculation########################
######################################################################

#' Estimates RPKM values and produces a plot for RNA_seq data (Supplementary Figure S6)
setwd("~/Documents/Phillips_lab/drafts/CR_genome/RNA_L1")
library(lsr)
library(dplyr)
library(ggplot2)

###Analyze L1 RNA for C. remanei and C. elegans

CEC<-read.csv("CE_L1_counts.txt",header=F,sep="\t") #raw RNA-seq counts
CRC<-read.csv("CR_L1_counts.txt",header=F,sep="\t") #raw RNA-seq counts

#CRPC00001.1	445
#CRPC00002.1	1203
#CRPC00003.1	2030
#CRPC00004.1	854
#CRPC00005.1	456
#CRPC00006.1	1075

CEP<-read.csv("CE.gene.BED",header=F,sep="\t") #Gene coordinates
CRP<-read.csv("CR.gene.BED",header=F,sep="\t") #Gene coordinates

#I	0	1094
#I	1407	4076
#I	6997	18219
#I	19098	21265
#I	22678	32516
#I	46683	47742

#'To calculate gene statistics from gff3 files I used genestats script from https://gist.github.com/darencard/fcb32168c243b92734e85c5f8b59a1c3

CEG<-read.csv("CE_genesize_stats.txt",header=F,sep="\t")
CRG<-read.csv("CR_geneStats.txt",header=F,sep="\t")

#CRPC00001.1	1094	2	209	1	884	2	209	0	0
#CRPC00002.1	2669	6	978	5	1686	6	978	0	0
#CRPC00003.1	11222	8	2038	7	9177	8	2038	0	0
#CRPC00004.1	2167	5	1501	4	662	5	1501	0	0
#CRPC00005.1	9838	8	949	7	8882	8	949	0	0
#CRPC00006.1	1059	2	478	1	580	2	478	0	0


#combine tables
CE<-cbind(CEG,CEP)
colnames(CE)<-c("ID","GENE_SIZE","N_EXONS","EXONS_TOTAL","N_INTRONS","INTRONS_TOTAL","N_CDS","CDS_TOTAL","5UTR","5UTR_TOTAL","3UTR","3UTR_TOTAL","CHR","START","END")
colnames(CEC)<-c("ID","COUNT")
CE<-merge(CE,CEC,by="ID",sort = FALSE,all.x =T)
CE$RPKM<-CE$COUNT/((CE$EXONS_TOTAL/1000)*(sum(CE$COUNT)/1000000)) #estimate RPKMs
CE$CHR<-gsub("_pilon","",CE$CHR)
CE$CHR<-gsub("chr","",CE$CHR)
CE<-CE[CE$CHR!="M",] #remove mitochondrial data

CRP<-CRP[CRP$V1!="MtDNA",] #remove mitochondrial data
CR<-cbind(CRG,CRP)
colnames(CR)<-c("ID","GENE_SIZE","N_EXONS","EXONS_TOTAL","N_INTRONS","INTRONS_TOTAL","N_CDS","CDS_TOTAL","5UTR","5UTR_TOTAL","3UTR","3UTR_TOTAL","CHR","START","END")
colnames(CRC)<-c("ID","COUNT")
CR<-merge(CR,CRC,by="ID",sort = FALSE, all.x=T)
CR$RPKM<-CR$COUNT/((CR$EXONS_TOTAL/1000)*(sum(CR$COUNT)/1000000)) #estimate RPKMs




#Sizes of C. remanei chromosomes
cr1=17247545
cr2=19935723
cr3=17877849
cr4=25790997
cr5=22502457
cr6=21501900

#Sizes of C. elegans chromosomes
ce1=15331301
ce2=15525148
ce3=14108536
ce4=17759200
ce5=21243235
ce6=18110855



CE$CHRSIZE<-0
CR$CHRSIZE<-0


CE[CE$CHR=="I" ,]$CHRSIZE<-ce1
CE[CE$CHR=="II" ,]$CHRSIZE<-ce2
CE[CE$CHR=="III" ,]$CHRSIZE<-ce3
CE[CE$CHR=="IV" ,]$CHRSIZE<-ce4
CE[CE$CHR=="V" ,]$CHRSIZE<-ce5
CE[CE$CHR=="X" ,]$CHRSIZE<-ce6

CR[CR$CHR=="I" ,]$CHRSIZE<-cr1
CR[CR$CHR=="II" ,]$CHRSIZE<-cr2
CR[CR$CHR=="III" ,]$CHRSIZE<-cr3
CR[CR$CHR=="IV" ,]$CHRSIZE<-cr4
CR[CR$CHR=="V" ,]$CHRSIZE<-cr5
CR[CR$CHR=="X" ,]$CHRSIZE<-cr6


CE$POSNORM<-abs((CE$START-CE$CHRSIZE/2)/CE$CHRSIZE)

CE$POSCLASS<-"CENTER"
CE[CE$POSNORM>=0.25,]$POSCLASS<-"ARM"

CR$POSNORM<-abs((CR$START-CR$CHRSIZE/2)/CR$CHRSIZE)

CR$POSCLASS<-"CENTER"
CR[CR$POSNORM>=0.25,]$POSCLASS<-"ARM"

#fix C. elegans counts with multiple isoforms
CEGENENAME<-read.csv("GENE_NAMES_TO_ISOFORMS.txt",sep="\t",header=F) #gene and transcript names
colnames(CEGENENAME)<-c("ID","GENE")
CE3<-merge(CE,CEGENENAME,by="ID",all.x=T)


 CE3<-CE3 %>%
   group_by(GENE,POSCLASS,CHR) %>%
   summarize(START=min(START), RPKM = sum(RPKM),COUNTS=sum(COUNT),intron=max(INTRONS_TOTAL),exons=max(EXONS_TOTAL),Species="C. elegans")
 CE3<-as.data.frame(CE3)
 CE3 %>%  filter_all(all_vars(!is.infinite(.))) ->CE3

 CE3$RPKM<-CE3$COUNTS/((CE3$exons/1000)*(sum(CE3$COUNTS)/1000000))
CR3<-CR[,c("ID","POSCLASS","CHR", "START","RPKM","COUNT","INTRONS_TOTAL","EXONS_TOTAL")]

CE3$Species<-"C.elegans"
CR3$Species<-"C.remanei"
colnames(CR3)<-colnames(CE3)

COMBO<-rbind(CE3,CR3)

######################################################################
############ Statistics from Supplementary Table S4 ##################
######################################################################

cohensD(RPKM ~POSCLASS, data=COMBO[COMBO$COUNT>10 & COMBO$Species=="C.elegans",], method="corrected")
cohensD(RPKM ~POSCLASS, data=COMBO[COMBO$COUNT>10 & COMBO$Species=="C.remanei",], method="corrected")

wilcox.test(
  COMBO[COMBO$COUNT>10 & COMBO$Species=="C.elegans" & COMBO$POSCLASS == "ARM",]$RPKM,
  COMBO[COMBO$COUNT>10 & COMBO$Species=="C.elegans" & COMBO$POSCLASS == "CENTER",]$RPKM, conf.int = T
)
wilcox.test(
  COMBO[COMBO$COUNT>10 & COMBO$Species=="C.remanei" & COMBO$POSCLASS == "ARM",]$RPKM,
  COMBO[COMBO$COUNT>10 & COMBO$Species=="C.remanei" & COMBO$POSCLASS == "CENTER",]$RPKM, conf.int = T
)

######################################################################
################# Supplementary figure S6 ############################
######################################################################

ggplot(COMBO[COMBO$COUNT>10,], aes_string(x = "POS", y = log(COMBO[COMBO$COUNT>10,]$RPKM), ordered = FALSE)) + facet_grid(Species~CHR, scale="free") +theme_bw() + labs(title ="", x = "Genome Position (Mb)", y = "RPKM") + theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5, size=8), plot.title = element_text(hjust = 0.5)) + theme(legend.position = "none")  + scale_x_continuous(labels = function(x) x/1000000) +theme(strip.background = element_rect(colour="white", fill="white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(strip.text= element_text(size = 10))
+geom_hex( bins = 35) + scale_fill_gradient(low = "grey40", high = "#FF7C01") ->rpkm

ggsave("CE_CR_RPKM.pdf",rpkm, width=7, height=3, units="in", scale=1.15)

#Additional figure to compare summary stats
#ggplot(COMBO[COMBO$COUNT>10,], aes_string(x = "POSCLASS", y = log(COMBO[COMBO$COUNT>10,]$RPKM), ordered = FALSE,color="POSCLASS")) + facet_grid(Species~CHR, space="free") +theme_bw() + labs(title ="", x = "Genome Position (Mb)", y = "Number of Exons") + theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5, size=8), plot.title = element_text(hjust = 0.5)) + theme(legend.position = "none")   +theme(strip.background = element_rect(colour="white", fill="white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + theme(strip.text= element_text(size = 10))  + scale_colour_manual(values = c(rep(c("#FF7C01","grey40"),3),"#4daf4a","#e31a1c", "#ff7f00","#984ea3","#a65628")) + scale_fill_manual(values = c(rep(c("#FF7C01","grey40"),3),"#377eb8","#4daf4a","#e31a1c", "#ff7f00","#984ea3","#a65628")) +geom_jitter(size=0.5, alpha=0.1,col="#FF7C01")+ geom_boxplot(size=0.5,outlier.shape = NA,col="grey40",fill=NA)




######################################################################
####### Stats for total introns lengths from Supplementary Table S4
######################################################################
 wilcox.test(
   CE[CE$POSCLASS == "ARM" & CE$intron>10 , ]$intron,
   CE[CE$POSCLASS == "CENTER" & CE$intron>10, ]$intron, conf.int = T
 )

 wilcox.test(
   CR[CR$POSCLASS == "ARM" & CE$intron>10 , ]$intron,
   CR[CR$POSCLASS == "CENTER" & CE$intron>10, ]$intron, conf.int = T
 )


cohensD(intron ~POSCLASS, data=CE[CE$intron>10,], method="corrected")
cohensD(intron ~POSCLASS, data=CR[CR$intron>10,], method="corrected")
