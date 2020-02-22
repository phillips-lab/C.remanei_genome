################################################################
###### Generates Figure 2, Figure S5, a Part of Table S4 #######
################################################################


setwd("~/Documents/Phillips_lab/drafts/CR_genome/Fractions")

library(lsr)
library(ggplot2)

###################
###Load the data###
###################

#exon/intron/repeats
CEE<-read.csv("CE_exon.bed",header=F,sep ="\t")
CEE<-CEE[,c("V1","V2","V3","V7")]
colnames(CEE)<-c("CHR","START","END","Value")
CEE$Type<-"Exon"
CEE$Species<-"C.elegans"

CEI<-read.csv("CE_intron_cov.bed",header=F,sep ="\t")
CEI<-CEI[,c("V1","V2","V3","V7")]
colnames(CEI)<-c("CHR","START","END","Value")
CEI$Type<-"Intron"
CEI$Species<-"C.elegans"
CEN<-read.csv("CE_N_content.bed",header=F,sep ="\t",skip=1)
CEN$NEW<-CEN$V10/100000
CEN<-CEN[,c("V1","V2","V3","NEW")]
colnames(CEN)<-c("CHR","START","END","Value")
CEN$Type<-"Repeatitive DNA"
CEN$Species<-"C.elegans"

CRE<-read.csv("CR_exon_cov.bed",header=F,sep ="\t")
CRE<-CRE[,c("V1","V2","V3","V7")]
colnames(CRE)<-c("CHR","START","END","Value")
CRE$Type<-"Exon"
CRE$Species<-"C.remanei"

CRI<-read.csv("CR_intron_cov.bed",header=F,sep ="\t")
CRI<-CRI[,c("V1","V2","V3","V7")]
colnames(CRI)<-c("CHR","START","END","Value")
CRI$Type<-"Intron"
CRI$Species<-"C.remanei"


CRN<-read.csv("CR_N_content.bed",header=F,sep ="\t",skip=1)
CRN$NEW<-CRN$V10/100000
CRN<-CRN[,c("V1","V2","V3","NEW")]
colnames(CRN)<-c("CHR","START","END","Value")
CRN$Type<-"Repeatitive DNA"
CRN$Species<-"C.remanei"

#gene/GC
CRG<-read.csv("CR_gene_cov.bed",header=F,sep ="\t")
CRG<-CRG[,c("V1","V2","V3","V7")]
colnames(CRG)<-c("CHR","START","END","Value")
CRG$Type<-"Gene"
CRG$Species<-"C.remanei"

CRGco<-read.csv("CR_gene_cov.bed",header=F,sep ="\t")
CRGco<-CRGco[,c("V1","V2","V3","V4")]
colnames(CRGco)<-c("CHR","START","END","Value")
CRGco$Type<-"Gene_count"
CRGco$Species<-"C.remanei"


CEGco<-read.csv("CE_gene.bed",header=F,sep ="\t")
CEGco<-CEGco[,c("V1","V2","V3","V4")]
colnames(CEGco)<-c("CHR","START","END","Value")
CEGco$Type<-"Gene_count"
CEGco$Species<-"C.elegans"


CRGC<-read.csv("CR_GC_content.bed",header=F,sep ="\t",skip=1)
CRGC<-CRGC[,c("V1","V2","V3","V5")]
colnames(CRGC)<-c("CHR","START","END","Value")
CRGC$Type<-"GC-content"
CRGC$Species<-"C.remanei"

CEG<-read.csv("CE_gene.bed",header=F,sep ="\t")
CEG<-CEG[,c("V1","V2","V3","V7")]
colnames(CEG)<-c("CHR","START","END","Value")
CEG$Type<-"Gene"
CEG$Species<-"C.elegans"

CEGC<-read.csv("CE_GC_content.bed",header=F,sep ="\t",skip=1)
CEGC<-CEGC[,c("V1","V2","V3","V5")]
colnames(CEGC)<-c("CHR","START","END","Value")
CEGC$Type<-"GC-content"
CEGC$Species<-"C.elegans"



######combine data


EIN<-rbind(CEE,CEI)
EIN<-rbind(EIN,CEN)
EIN<-rbind(EIN,CRE)
EIN<-rbind(EIN,CRI)
EIN<-rbind(EIN,CRN)
EIN$CHR<-gsub("_pilon","",EIN$CHR)
EIN$CHR<-gsub("chr","",EIN$CHR)
EIN$CHR<-gsub("MtDNA","M",EIN$CHR)
EIN$POS<-(EIN$START+EIN$END)/2



GGC<-rbind(CEG,CEGC)
GGC<-rbind(GGC,CRG)
GGC<-rbind(GGC,CRGC)
GGC$POS<-(GGC$START+GGC$END)/2

GGC$CHR<-gsub("_pilon","",GGC$CHR)
GGC$CHR<-gsub("chr","",GGC$CHR)
GGC$CHR<-gsub("MtDNA","M",GGC$CHR)


GENECOUNT<-rbind(CEGco,CRGco)
GENECOUNT$POS<-(GENECOUNT$START+GENECOUNT$END)/2
GENECOUNT$CHR<-gsub("_pilon","",GENECOUNT$CHR)
GENECOUNT$CHR<-gsub("chr","",GENECOUNT$CHR)
GENECOUNT$CHR<-gsub("MtDNA","M",GENECOUNT$CHR)



####################################
############# FIGURES ##############
####################################

CROMOSOME <-
  data.frame(CHR = rep(c("I","I","II", "II", "III", "III","IV","IV", "V","V","X", "X"), 2),
             Boundary = c(cr1 / 4, cr1 * 3/ 4, cr2 / 4, cr2 * 3 / 4, cr3 / 4, cr3 * 3 / 4, cr4 / 4, cr4 * 3 / 4, cr5 / 4, cr5 * 3 / 4, cr6 / 4, cr6 * 3 / 4, ce1 / 4, ce1 * 3 / 4, ce2 / 4, ce2 * 3/ 4, ce3 / 4, ce3 * 3 / 4, ce4 / 4, ce4 * 3 / 4, ce5 / 4, ce5 * 3 / 4, ce6 / 4, ce6 * 3 / 4), Species=c(rep("C.remanei",12),rep("C.elegans",12)))




#################################################
########exon/intron/repeat Figure 2 #############
#################################################


ggplot(EIN[EIN$CHR!="M" ,], aes_string(x = "POS", y = EIN[EIN$CHR!="M" ,]$Value*100, color= "Type", fill="Type", ordered = FALSE)) + facet_grid(Species~CHR,scales  ="free") +theme_bw() + labs(title ="", x = "Genome Position (Mb)", y = "Fraction (%)") + theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0.5)) + theme(legend.position = "none") + scale_colour_manual(values = c(rep(c("#FF7C01","#622278","grey70"),3),"#4daf4a","#e31a1c", "#ff7f00","#984ea3","#a65628"),name="") + scale_fill_manual(values = c(rep(c("#FF7C01","#622278","grey70"),3),"#377eb8","#4daf4a","#e31a1c", "#ff7f00","#984ea3","#a65628"),name="") + scale_x_continuous(labels = function(x) x/1000000) +theme(strip.background = element_rect(colour="white", fill="white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + geom_smooth() + theme(strip.text= element_text(size = 10)) + geom_vline(data=CROMOSOME,aes(xintercept =Boundary), colour="grey50",size=0.5,linetype="dashed") ->eir

ggsave("CE_CR_exon_intron_repeats_content_NOlegend.pdf",eir, width=7, height=3, units="in", scale=1.15)


############################################
######### GC-content, Figure S5 A ##########
############################################

ggplot(GGC[GGC$CHR!="M" & GGC$Type=="GC-content",], aes_string(x = "POS", y = GGC[GGC$CHR!="M" & GGC$Type=="GC-content",]$Value*100, color= "Type", fill="Type", ordered = FALSE)) + facet_grid(Species~CHR, scale="free") +theme_bw() + labs(title ="", x = "Genome Position (Mb)", y = "GC-content (%)") + theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0.5)) + theme(legend.position = "none") + scale_colour_manual(values = c(rep(c("#A49C74","grey40"),3),"#4daf4a","#e31a1c", "#ff7f00","#984ea3","#a65628")) + scale_fill_manual(values = c(rep(c("#A49C74","grey40"),3),"#377eb8","#4daf4a","#e31a1c", "#ff7f00","#984ea3","#a65628")) + scale_x_continuous(labels = function(x) x/1000000) +theme(strip.background = element_rect(colour="white", fill="white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + geom_point(size=0.75,alpha=0.75) + geom_smooth(color="grey30",linetype="solid",se=F,size=0.6) + theme(strip.text= element_text(size = 10))  + theme(strip.text= element_text(size = 10)) + geom_vline(data=CROMOSOME,aes(xintercept =Boundary), colour="grey50",size=0.5,linetype="dashed") ->gc

ggsave("CE_CR_GC_content.pdf",gc, width=7, height=3, units="in", scale=1.15)



#################################################
########## Gene fraction, Figure S5 B ###########
#################################################

ggplot(GGC[GGC$CHR!="M" & GGC$Type=="Gene",], aes_string(x = "POS", y = GGC[GGC$CHR!="M" & GGC$Type=="Gene",]$Value*100, color= "Type", fill="Type", ordered = FALSE)) + facet_grid(Species~CHR, scale="free") +theme_bw() + labs(title ="", x = "Genome Position (Mb)", y = "Gene Fraction (%)") + theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0.5)) + theme(legend.position = "none") + scale_colour_manual(values = c(rep(c("#FF7C01","grey40"),3),"#4daf4a","#e31a1c", "#ff7f00","#984ea3","#a65628")) + scale_fill_manual(values = c(rep(c("#FF7C01","grey40"),3),"#377eb8","#4daf4a","#e31a1c", "#ff7f00","#984ea3","#a65628")) + scale_x_continuous(labels = function(x) x/1000000) +theme(strip.background = element_rect(colour="white", fill="white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + geom_point(size=0.75,alpha=0.75) + geom_smooth(color="grey30",linetype="solid",se=F,size=0.6) + theme(strip.text= element_text(size = 10))  + theme(strip.text= element_text(size = 10)) + geom_vline(data=CROMOSOME,aes(xintercept =Boundary), colour="grey50",size=0.5,linetype="dashed") ->gf

ggsave("CE_CR_gene_fraction.pdf",gf, width=7, height=3, units="in", scale=1.15)



####################################################
######### The number of genes, Figure S5 C##########
####################################################


ggplot(GENECOUNT[GENECOUNT$CHR!="M",], aes_string(x = "POS", y = GENECOUNT[GENECOUNT$CHR!="M" ,]$Value, color= "Type", fill="Type", ordered = FALSE)) + facet_grid(Species~CHR, scale="free") +theme_bw() + labs(title ="", x = "Genome Position (Mb)", y = "Number of genes") + theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0.5)) + theme(legend.position = "none") + scale_colour_manual(values = c(rep(c("#A0750C","grey40"),3),"#4daf4a","#e31a1c", "#ff7f00","#984ea3","#a65628")) + scale_fill_manual(values = c(rep(c("#A0750C","grey40"),3),"#377eb8","#4daf4a","#e31a1c", "#ff7f00","#984ea3","#a65628")) + scale_x_continuous(labels = function(x) x/1000000) +theme(strip.background = element_rect(colour="white", fill="white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + geom_point(size=0.75,alpha=0.75) + geom_smooth(color="grey30",linetype="solid",se=F,size=0.6) + theme(strip.text= element_text(size = 10))  + theme(strip.text= element_text(size = 10)) + geom_vline(data=CROMOSOME,aes(xintercept =Boundary), colour="grey50",size=0.5,linetype="dashed") ->gn
ggsave("CE_CR_gene_number.pdf",gn, width=7, height=3, units="in", scale=1.15)



########################################################
############# Arm/Center domain statistics #############
########################################################

#normalize positions

cr1=17247545
cr2=19935723
cr3=17877849
cr4=25790997
cr5=22502457
cr6=21501900

ce1=15331301
ce2=15525148
ce3=14108536
ce4=17759200
ce5=21243235
ce6=18110855

COMBO<-rbind(EIN,GGC)
COMBO<-rbind(COMBO,GENECOUNT)
COMBO<-COMBO[COMBO$CHR!="M",]


COMBO$CHRSIZE<-0

COMBO[COMBO$Species=="C.elegans" & COMBO$CHR=="I",]$CHRSIZE<-ce1
COMBO[COMBO$Species=="C.elegans" & COMBO$CHR=="II",]$CHRSIZE<-ce2
COMBO[COMBO$Species=="C.elegans" & COMBO$CHR=="III",]$CHRSIZE<-ce3
COMBO[COMBO$Species=="C.elegans" & COMBO$CHR=="IV",]$CHRSIZE<-ce4
COMBO[COMBO$Species=="C.elegans" & COMBO$CHR=="V",]$CHRSIZE<-ce5
COMBO[COMBO$Species=="C.elegans" & COMBO$CHR=="X",]$CHRSIZE<-ce6

COMBO[COMBO$Species=="C.remanei" & COMBO$CHR=="I",]$CHRSIZE<-cr1
COMBO[COMBO$Species=="C.remanei" & COMBO$CHR=="II",]$CHRSIZE<-cr2
COMBO[COMBO$Species=="C.remanei" & COMBO$CHR=="III",]$CHRSIZE<-cr3
COMBO[COMBO$Species=="C.remanei" & COMBO$CHR=="IV",]$CHRSIZE<-cr4
COMBO[COMBO$Species=="C.remanei" & COMBO$CHR=="V",]$CHRSIZE<-cr5
COMBO[COMBO$Species=="C.remanei" & COMBO$CHR=="X",]$CHRSIZE<-cr6


COMBO$POSNORM<-abs((COMBO$POS-COMBO$CHRSIZE/2)/COMBO$CHRSIZE)

COMBO$POSCLASS<-"CENTER"
COMBO[COMBO$POSNORM>=0.25,]$POSCLASS<-"ARM"


######################################
######## effect sizes ################
######################################

cohensD(Value ~POSCLASS, data=COMBO[COMBO$Species=="C.elegans" & COMBO$Type=="Gene_count",], method="corrected")
cohensD(Value ~POSCLASS, data=COMBO[COMBO$Species=="C.remanei" & COMBO$Type=="Gene_count",], method="corrected")
cohensD(Value ~POSCLASS, data=COMBO[COMBO$Species=="C.elegans" & COMBO$Type=="GC-content",], method="corrected")
cohensD(Value ~POSCLASS, data=COMBO[COMBO$Species=="C.remanei" & COMBO$Type=="GC-content",], method="corrected")
cohensD(Value ~POSCLASS, data=COMBO[COMBO$Species=="C.elegans" & COMBO$Type=="Repetitive DNA",], method="corrected")
cohensD(Value ~POSCLASS, data=COMBO[COMBO$Species=="C.remanei" & COMBO$Type=="Repetitive DNA",], method="corrected")
cohensD(Value ~POSCLASS, data=COMBO[COMBO$Species=="C.elegans" & COMBO$Type=="Exon",], method="corrected")
cohensD(Value ~POSCLASS, data=COMBO[COMBO$Species=="C.remanei" & COMBO$Type=="Exon",], method="corrected")
cohensD(Value ~POSCLASS, data=COMBO[COMBO$Species=="C.elegans" & COMBO$Type=="Intron",], method="corrected")
cohensD(Value ~POSCLASS, data=COMBO[COMBO$Species=="C.remanei" & COMBO$Type=="Intron",], method="corrected")
cohensD(Value ~POSCLASS, data=COMBO[COMBO$Species=="C.elegans" & COMBO$Type=="Gene",], method="corrected")
cohensD(Value ~POSCLASS, data=COMBO[COMBO$Species=="C.remanei" & COMBO$Type=="Gene",], method="corrected")


################################
##### Mann_Whitney test ########
################################

wilcox.test(
  COMBO[COMBO$Species == "C.elegans" &
          COMBO$Type == "Gene" &
          COMBO$POSCLASS == "ARM" , ]$Value,
  COMBO[COMBO$Species == "C.elegans" &
          COMBO$Type == "Gene" &
          COMBO$POSCLASS == "CENTER" , ]$Value, conf.int = T
)

wilcox.test(
  COMBO[COMBO$Species == "C.elegans" &
          COMBO$Type == "Gene_count" &
          COMBO$POSCLASS == "ARM" , ]$Value,
  COMBO[COMBO$Species == "C.elegans" &
          COMBO$Type == "Gene_count" &
          COMBO$POSCLASS == "CENTER" , ]$Value, conf.int = T
)

wilcox.test(
  COMBO[COMBO$Species == "C.elegans" &
          COMBO$Type == "Exon" &
          COMBO$POSCLASS == "ARM" , ]$Value,
  COMBO[COMBO$Species == "C.elegans" &
          COMBO$Type == "Exon" &
          COMBO$POSCLASS == "CENTER" , ]$Value, conf.int = T
)
wilcox.test(
  COMBO[COMBO$Species == "C.elegans" &
          COMBO$Type == "Intron" &
          COMBO$POSCLASS == "ARM" , ]$Value,
  COMBO[COMBO$Species == "C.elegans" &
          COMBO$Type == "Intron" &
          COMBO$POSCLASS == "CENTER" , ]$Value, conf.int = T
)
wilcox.test(
  COMBO[COMBO$Species == "C.elegans" &
          COMBO$Type == "GC-content" &
          COMBO$POSCLASS == "ARM" , ]$Value,
  COMBO[COMBO$Species == "C.elegans" &
          COMBO$Type == "GC-content" &
          COMBO$POSCLASS == "CENTER" , ]$Value, conf.int = T
)
wilcox.test(
  COMBO[COMBO$Species == "C.elegans" &
          COMBO$Type == "Repetitive DNA" &
          COMBO$POSCLASS == "ARM" , ]$Value,
  COMBO[COMBO$Species == "C.elegans" &
          COMBO$Type == "Repetitive DNA" &
          COMBO$POSCLASS == "CENTER" , ]$Value, conf.int = T
)


#C.remanei

wilcox.test(
  COMBO[COMBO$Species == "C.remanei" &
          COMBO$Type == "Gene" &
          COMBO$POSCLASS == "ARM" , ]$Value,
  COMBO[COMBO$Species == "C.remanei" &
          COMBO$Type == "Gene" &
          COMBO$POSCLASS == "CENTER" , ]$Value, conf.int = T
)

wilcox.test(
  COMBO[COMBO$Species == "C.remanei" &
          COMBO$Type == "Gene_count" &
          COMBO$POSCLASS == "ARM" , ]$Value,
  COMBO[COMBO$Species == "C.remanei" &
          COMBO$Type == "Gene_count" &
          COMBO$POSCLASS == "CENTER" , ]$Value, conf.int = T
)

wilcox.test(
  COMBO[COMBO$Species == "C.remanei" &
          COMBO$Type == "Exon" &
          COMBO$POSCLASS == "ARM" , ]$Value,
  COMBO[COMBO$Species == "C.remanei" &
          COMBO$Type == "Exon" &
          COMBO$POSCLASS == "CENTER" , ]$Value, conf.int = T
)
wilcox.test(
  COMBO[COMBO$Species == "C.remanei" &
          COMBO$Type == "Intron" &
          COMBO$POSCLASS == "ARM" , ]$Value,
  COMBO[COMBO$Species == "C.remanei" &
          COMBO$Type == "Intron" &
          COMBO$POSCLASS == "CENTER" , ]$Value, conf.int = T
)
wilcox.test(
  COMBO[COMBO$Species == "C.remanei" &
          COMBO$Type == "GC-content" &
          COMBO$POSCLASS == "ARM" , ]$Value,
  COMBO[COMBO$Species == "C.remanei" &
          COMBO$Type == "GC-content" &
          COMBO$POSCLASS == "CENTER" , ]$Value, conf.int = T
)
wilcox.test(
  COMBO[COMBO$Species == "C.remanei" &
          COMBO$Type == "Repetitive DNA" &
          COMBO$POSCLASS == "ARM" , ]$Value,
  COMBO[COMBO$Species == "C.remanei" &
          COMBO$Type == "Repetitive DNA" &
          COMBO$POSCLASS == "CENTER" , ]$Value, conf.int = T
)
