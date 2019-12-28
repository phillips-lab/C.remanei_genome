######################################################################
################ Produces Figure 3 and Hi-C stats ####################
######################################################################

setwd("~/Documents/Phillips_lab/drafts/CR_genome/HIC/")
library(dplyr)
library(lsr)
library(ggplot2)
library(circlize)
library(RColorBrewer)


######################################################################
#################### Interchromosomal interactions ###################
######################################################################
CEDIF<-read.csv("CE_HIC.differ",header =F,sep="\t") #see Hi-C_analysis_with_ARIMA.sh
CRDIF<-read.csv("CR_HIC.differ",header =F,sep="\t") #see Hi-C_analysis_with_ARIMA.sh


colnames(CEDIF)<-c("COUNT","CHR1","POS1","CHR2","POS2")
colnames(CRDIF)<-c("COUNT","CHR1","POS1","CHR2","POS2")

CEDIF<-CEDIF[CEDIF$CHR2!="7",]
CRDIF<-CRDIF[CRDIF$CHR2!="7",]

CEDIF$CHRSIZE1<-0
CEDIF$CHRSIZE2<-0
CRDIF$CHRSIZE1<-0
CRDIF$CHRSIZE2<-0

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

CEDIF[CEDIF$CHR1==1,]$CHRSIZE1<-ce1
CEDIF[CEDIF$CHR1==2,]$CHRSIZE1<-ce2
CEDIF[CEDIF$CHR1==3,]$CHRSIZE1<-ce3
CEDIF[CEDIF$CHR1==4,]$CHRSIZE1<-ce4
CEDIF[CEDIF$CHR1==5,]$CHRSIZE1<-ce5
#CEDIF[CEDIF$CHR1==6,]$CHRSIZE1<-ce6
#CEDIF[CEDIF$CHR2==1,]$CHRSIZE2<-ce1
CEDIF[CEDIF$CHR2==2,]$CHRSIZE2<-ce2
CEDIF[CEDIF$CHR2==3,]$CHRSIZE2<-ce3
CEDIF[CEDIF$CHR2==4,]$CHRSIZE2<-ce4
CEDIF[CEDIF$CHR2==5,]$CHRSIZE2<-ce5
CEDIF[CEDIF$CHR2==6,]$CHRSIZE2<-ce6



CRDIF[CRDIF$CHR1==1,]$CHRSIZE1<-cr1
CRDIF[CRDIF$CHR1==2,]$CHRSIZE1<-cr2
CRDIF[CRDIF$CHR1==3,]$CHRSIZE1<-cr3
CRDIF[CRDIF$CHR1==4,]$CHRSIZE1<-cr4
CRDIF[CRDIF$CHR1==5,]$CHRSIZE1<-cr5
#CRDIF[CRDIF$CHR1==6,]$CHRSIZE1<-cr6
#CRDIF[CRDIF$CHR2==1,]$CHRSIZE2<-cr1
CRDIF[CRDIF$CHR2==2,]$CHRSIZE2<-cr2
CRDIF[CRDIF$CHR2==3,]$CHRSIZE2<-cr3
CRDIF[CRDIF$CHR2==4,]$CHRSIZE2<-cr4
CRDIF[CRDIF$CHR2==5,]$CHRSIZE2<-cr5
CRDIF[CRDIF$CHR2==6,]$CHRSIZE2<-cr6




CEDIF$POSNORM1<-abs((CEDIF$POS1*100000-CEDIF$CHRSIZE1/2)/CEDIF$CHRSIZE1)
CEDIF$POSNORM2<-abs((CEDIF$POS2*100000-CEDIF$CHRSIZE2/2)/CEDIF$CHRSIZE2)

CEDIF$POSCLASS1<-"CENTER"
CEDIF[CEDIF$POSNORM1>=0.25,]$POSCLASS1<-"ARM"
CEDIF$POSCLASS2<-"CENTER"
CEDIF[CEDIF$POSNORM2>=0.25,]$POSCLASS2<-"ARM"

CRDIF$POSNORM1<-abs((CRDIF$POS1*100000-CRDIF$CHRSIZE1/2)/CRDIF$CHRSIZE1)
CRDIF$POSNORM2<-abs((CRDIF$POS2*100000-CRDIF$CHRSIZE2/2)/CRDIF$CHRSIZE2)

CRDIF$POSCLASS1<-"CENTER"
CRDIF[CRDIF$POSNORM1>=0.25,]$POSCLASS1<-"ARM"
CRDIF$POSCLASS2<-"CENTER"
CRDIF[CRDIF$POSNORM2>=0.25,]$POSCLASS2<-"ARM"

CEGROUP<-c(sum(CEDIF[(CEDIF$POSCLASS1=="ARM" &CEDIF$POSCLASS2=="ARM"),  ]$COUNT),sum(CEDIF[(CEDIF$POSCLASS1=="ARM" &CEDIF$POSCLASS2=="CENTER")|(CEDIF$POSCLASS1=="CENTER" &CEDIF$POSCLASS2=="ARM"),  ]$COUNT),sum(CEDIF[(CEDIF$POSCLASS1=="CENTER" &CEDIF$POSCLASS2=="CENTER"),]$COUNT))
CRGROUP<-c(sum(CRDIF[(CRDIF$POSCLASS1=="ARM" &CRDIF$POSCLASS2=="ARM"),  ]$COUNT),sum(CRDIF[(CRDIF$POSCLASS1=="ARM" &CRDIF$POSCLASS2=="CENTER")|(CRDIF$POSCLASS1=="CENTER" &CRDIF$POSCLASS2=="ARM"),  ]$COUNT),sum(CRDIF[(CRDIF$POSCLASS1=="CENTER" &CRDIF$POSCLASS2=="CENTER"),]$COUNT))


chisq.test(CEGROUP, p = c(1/4, 1/2, 1/4))
chisq.test(CRGROUP, p = c(1/4, 1/2, 1/4))

#without X chromosome

CEDIFNX<-CEDIF[CEDIF$CHR2!="6",]
CRDIFNX<-CRDIF[CRDIF$CHR2!="6",]

CEGROUPNX<-c(sum(CEDIFNX[(CEDIFNX$POSCLASS1=="ARM" &CEDIFNX$POSCLASS2=="ARM"),  ]$COUNT),sum(CEDIFNX[(CEDIFNX$POSCLASS1=="ARM" &CEDIFNX$POSCLASS2=="CENTER")|(CEDIFNX$POSCLASS1=="CENTER" &CEDIFNX$POSCLASS2=="ARM"),  ]$COUNT),sum(CEDIFNX[(CEDIFNX$POSCLASS1=="CENTER" &CEDIFNX$POSCLASS2=="CENTER"),]$COUNT))
CRGROUPNX<-c(sum(CRDIFNX[(CRDIFNX$POSCLASS1=="ARM" &CRDIFNX$POSCLASS2=="ARM"),  ]$COUNT),sum(CRDIFNX[(CRDIFNX$POSCLASS1=="ARM" &CRDIFNX$POSCLASS2=="CENTER")|(CRDIFNX$POSCLASS1=="CENTER" &CRDIFNX$POSCLASS2=="ARM"),  ]$COUNT),sum(CRDIFNX[(CRDIFNX$POSCLASS1=="CENTER" &CRDIFNX$POSCLASS2=="CENTER"),]$COUNT))
chisq.test(CEGROUPNX, p = c(1/4, 1/2, 1/4))
chisq.test(CRGROUPNX, p = c(1/4, 1/2, 1/4))

######################################################################
############### Figure 3 B, Circos plot, Hi-C ########################
######################################################################

#Generates Circos plots for Hi-C reads mapped to different chromosomes


#library(circlize)
#library(RColorBrewer)



CR_dr<-data.frame(name=c("I","II","III","IV","V","X"), start=rep(1,6), end=c(cr1,cr2,cr3,cr4,cr5,cr6))
CE_dr<-data.frame(name=c("I","II","III","IV","V","X"), start=rep(1,6), end=c(ce1,ce2,ce3,ce4,ce5,ce6))


CRDIF2<-CRDIF
CRDIF2$CHR1<-gsub("1","I",CRDIF2$CHR1)
CRDIF2$CHR1<-gsub("2","II",CRDIF2$CHR1)
CRDIF2$CHR1<-gsub("3","III",CRDIF2$CHR1)
CRDIF2$CHR1<-gsub("4","IV",CRDIF2$CHR1)
CRDIF2$CHR1<-gsub("5","V",CRDIF2$CHR1)
CRDIF2$CHR1<-gsub("6","X",CRDIF2$CHR1)
CRDIF2$CHR2<-gsub("1","I",CRDIF2$CHR2)
CRDIF2$CHR2<-gsub("2","II",CRDIF2$CHR2)
CRDIF2$CHR2<-gsub("3","III",CRDIF2$CHR2)
CRDIF2$CHR2<-gsub("4","IV",CRDIF2$CHR2)
CRDIF2$CHR2<-gsub("5","V",CRDIF2$CHR2)
CRDIF2$CHR2<-gsub("6","X",CRDIF2$CHR2)

CEDIF2<-CEDIF
CEDIF2$CHR1<-gsub("1","I",CEDIF2$CHR1)
CEDIF2$CHR1<-gsub("2","II",CEDIF2$CHR1)
CEDIF2$CHR1<-gsub("3","III",CEDIF2$CHR1)
CEDIF2$CHR1<-gsub("4","IV",CEDIF2$CHR1)
CEDIF2$CHR1<-gsub("5","V",CEDIF2$CHR1)
CEDIF2$CHR1<-gsub("6","X",CEDIF2$CHR1)

CEDIF2$CHR2<-gsub("1","I",CEDIF2$CHR2)
CEDIF2$CHR2<-gsub("2","II",CEDIF2$CHR2)
CEDIF2$CHR2<-gsub("3","III",CEDIF2$CHR2)
CEDIF2$CHR2<-gsub("4","IV",CEDIF2$CHR2)
CEDIF2$CHR2<-gsub("5","V",CEDIF2$CHR2)
CEDIF2$CHR2<-gsub("6","X",CEDIF2$CHR2)

#filter based on the coverage
CRDIF2_FILT<-CRDIF2[CRDIF2$COUNT>100,]
CEDIF2_FILT<-CEDIF2[CEDIF2$COUNT>200,]


#C. remanei
circos.clear()
dev.off()
circos.par(start.degree=85,gap.degree=3)
circos.genomicInitialize(CR_dr,major.by = 6000000,plotType = "labels",labels.cex=1.5)

circos.track(ylim = c(0, 1), bg.col = "#DFA675", bg.border = NA, track.height = 0.067)

for (i in 1:nrow(CRDIF2_FILT)) {
  circos.link(
    as.character(CRDIF2_FILT$CHR1[i]),
    c(CRDIF2_FILT$POS1[i] * 100000, CRDIF2_FILT$POS1[i] * 100000 + 100000),
    as.character(CRDIF2_FILT$CHR2[i]),
    c(CRDIF2_FILT$POS2[i] * 100000, CRDIF2_FILT$POS2[i] * 100000 + 100000),
    col = "#DFA67550", border=FALSE
  )
}

recordPlot()->plot_cr


#C. elegans
circos.clear()
dev.off()
circos.par(start.degree=80,gap.degree=3)
circos.genomicInitialize(CE_dr,major.by = 6000000,plotType = "labels",labels.cex=1.64)

circos.track(ylim = c(0, 1), bg.col = "#DFA675", bg.border = NA, track.height = 0.075)


for (i in 1:nrow(CEDIF2_FILT)) {
  circos.link(
    as.character(CEDIF2_FILT$CHR1[i]),
    c(CEDIF2_FILT$POS1[i] * 100000, CEDIF2_FILT$POS1[i] * 100000 + 100000),
    as.character(CEDIF2_FILT$CHR2[i]),
    c(CEDIF2_FILT$POS2[i] * 100000, CEDIF2_FILT$POS2[i] * 100000 + 100000),
    col = "#DFA67550", border=FALSE
  )
}

recordPlot()->plot_ce




############################################
#### Intrachromosomal interactions #########
############################################


CESAME<-read.csv("CE_HIC_100.bed",header =F,sep="\t")
CRSAME<-read.csv("CR_HIC_100.bed",header =F,sep="\t")

colnames(CESAME)<-c("CHR","START","END","Value")
colnames(CRSAME)<-c("CHR","START","END","Value")

CESAME$Species<-"C.elegans"
CRSAME$Species<-"C.remanei"

CESAME$CHR<-gsub("_pilon","",CESAME$CHR)
CESAME$CHR<-gsub("chr","",CESAME$CHR)
CESAME<-CESAME[CESAME$CHR!="M",]
CRSAME<-CRSAME[CRSAME$CHR!="MtDNA",]

SAME<-rbind(CESAME,CRSAME)
SAME$POS<-(SAME$START+SAME$END)/2
SAME<-SAME[!is.na(SAME$Value),]
SAME<-as.data.frame(SAME)
#EXCLUDE ONE OUTLIER (IT HAS TO FEW SITES + TOO BIX MEDIAN)
SAME<-SAME[SAME$Value!=4257171,]
SAME$Value<-as.numeric(as.character(SAME$Value))


#add vline in the plot
CROMOSOME <-
  data.frame(CHR = rep(c("I","I","II", "II", "III", "III","IV","IV", "V","V","X", "X"), 2),
             Boundary = c(cr1 / 4, cr1 * 3/ 4, cr2 / 4, cr2 * 3 / 4, cr3 / 4, cr3 * 3 / 4, cr4 / 4, cr4 * 3 / 4, cr5 / 4, cr5 * 3 / 4, cr6 / 4, cr6 * 3 / 4, ce1 / 4, ce1 * 3 / 4, ce2 / 4, ce2 * 3/ 4, ce3 / 4, ce3 * 3 / 4, ce4 / 4, ce4 * 3 / 4, ce5 / 4, ce5 * 3 / 4, ce6 / 4, ce6 * 3 / 4), Species=c(rep("C.remanei",12),rep("C.elegans",12)))



######################################################################
######################### Figure 3 A  ################################
######################################################################

ggplot(SAME, aes_string(x = "POS", y = (SAME$Value), color= "Species", fill="Species", ordered = FALSE)) + facet_grid(Species~CHR, space="free",scale="free") +theme_bw() + labs(title ="", x = "Genome Position (Mb)", y = "Median distance between Hi-C reads") + theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0.5)) + theme(legend.position = "none") + scale_colour_manual(values = c(rep(c("#FF7C01","grey40"),3),"#4daf4a","#e31a1c", "#ff7f00","#984ea3","#a65628")) + scale_fill_manual(values = c(rep(c("#FF7C01","grey40"),3),"#377eb8","#4daf4a","#e31a1c", "#ff7f00","#984ea3","#a65628")) + scale_x_continuous(labels = function(x) x/1000000) +theme(strip.background = element_rect(colour="white", fill="white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + geom_point(size=0.75,col="#489DA9",alpha=0.75) + geom_smooth(color="grey30",linetype="solid",se=F,size=0.6) + theme(strip.text= element_text(size = 10)) + geom_vline(data=CROMOSOME,aes(xintercept =Boundary), colour="grey50",size=0.5,linetype="dashed")
#ggsave("CE_CR_HIC_same.pdf",s, width=7, height=3, units="in", scale=1.15)



#add arm-center classes
SAME$CHRSIZE<-0
SAME[SAME$CHR=="I" & SAME$Species=="C.elegans",]$CHRSIZE<-ce1
SAME[SAME$CHR=="II" & SAME$Species=="C.elegans",]$CHRSIZE<-ce2
SAME[SAME$CHR=="III" & SAME$Species=="C.elegans",]$CHRSIZE<-ce3
SAME[SAME$CHR=="IV" & SAME$Species=="C.elegans",]$CHRSIZE<-ce4
SAME[SAME$CHR=="V" & SAME$Species=="C.elegans",]$CHRSIZE<-ce5
SAME[SAME$CHR=="X" & SAME$Species=="C.elegans",]$CHRSIZE<-ce6


SAME[SAME$CHR=="I" & SAME$Species=="C.remanei",]$CHRSIZE<-cr1
SAME[SAME$CHR=="II" & SAME$Species=="C.remanei",]$CHRSIZE<-cr2
SAME[SAME$CHR=="III" & SAME$Species=="C.remanei",]$CHRSIZE<-cr3
SAME[SAME$CHR=="IV" & SAME$Species=="C.remanei",]$CHRSIZE<-cr4
SAME[SAME$CHR=="V" & SAME$Species=="C.remanei",]$CHRSIZE<-cr5
SAME[SAME$CHR=="X" & SAME$Species=="C.remanei",]$CHRSIZE<-cr6

SAME$POSNORM<-abs((SAME$POS-SAME$CHRSIZE/2)/SAME$CHRSIZE)

SAME$POSCLASS<-"CENTER"
SAME[SAME$POSNORM>=0.25,]$POSCLASS<-"ARM"

#######################################################
########From supplementary Table S4 ###################
#######################################################


wilcox.test(
  SAME[SAME$Species == "C.elegans" &
         SAME$POSCLASS == "ARM" , ]$Value,
  SAME[SAME$Species == "C.elegans" &
         SAME$POSCLASS == "CENTER" , ]$Value, conf.int = T

)

wilcox.test(
  SAME[SAME$Species == "C.remanei" &
         SAME$POSCLASS == "ARM" , ]$Value,
  SAME[SAME$Species == "C.remanei" &
         SAME$POSCLASS == "CENTER" , ]$Value, conf.int = T

)


cohensD(Value ~POSCLASS, data=SAME[SAME$Species=="C.elegans",], method="corrected")
cohensD(Value ~POSCLASS, data=SAME[SAME$Species=="C.remanei",], method="corrected")


###########################################################
##### Proportions of different types of interactions ######
###########################################################


CRDIF%>% group_by(CHR1,CHR2,POSCLASS1,POSCLASS2) %>% summarise(Counts = sum(COUNT)) ->CRDIFTABLE
CEDIF%>% group_by(CHR1,CHR2,POSCLASS1,POSCLASS2) %>% summarise(Counts = sum(COUNT)) ->CEDIFTABLE


CRDIFTABLE<-as.data.frame(CRDIFTABLE)
CEDIFTABLE<-as.data.frame(CEDIFTABLE)



CRDIF%>% group_by(CHR1,CHR2) %>% summarise(Counts = sum(COUNT)) ->CRDIFTABLES
CEDIF%>% group_by(CHR1,CHR2) %>% summarise(Counts = sum(COUNT)) ->CEDIFTABLES


CRDIFTABLES<-as.data.frame(CRDIFTABLES)
CEDIFTABLES<-as.data.frame(CEDIFTABLES)

crprop<-c();for (i in 1:6){print(i); a<-sum(CRDIFTABLES[CRDIFTABLES$CHR1==i|CRDIFTABLES$CHR2==i,]$Counts/sum(CRDIFTABLES$Counts)); crprop<-c(crprop,a)}
ceprop<-c();for (i in 1:6){print(i); a<-sum(CEDIFTABLES[CEDIFTABLES$CHR1==i|CEDIFTABLES$CHR2==i,]$Counts/sum(CEDIFTABLES$Counts)); ceprop<-c(ceprop,a)}



CRDIF2_FILT%>% group_by(CHR1,CHR2) %>% summarise(Counts = sum(COUNT)) ->CRDIFTABLEF
CEDIF2_FILT%>% group_by(CHR1,CHR2) %>% summarise(Counts = sum(COUNT)) ->CEDIFTABLEF


CRDIFTABLEF<-as.data.frame(CRDIFTABLEF)
CEDIFTABLEF<-as.data.frame(CEDIFTABLEF)

crpropf<-c();for (i in c("I","II","III","IV","V","X")){print(i); a<-sum(CRDIFTABLEF[CRDIFTABLEF$CHR1==i|CRDIFTABLEF$CHR2==i,]$Counts/sum(CRDIFTABLEF$Counts)); crpropf<-c(crpropf,a)}
cepropf<-c();for (i in c("I","II","III","IV","V","X")){print(i); a<-sum(CEDIFTABLEF[CEDIFTABLEF$CHR1==i|CEDIFTABLEF$CHR2==i,]$Counts/sum(CEDIFTABLEF$Counts)); cepropf<-c(cepropf,a)}
