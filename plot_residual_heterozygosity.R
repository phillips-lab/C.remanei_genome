setwd("/Users/anastasia/Documents/Phillips_lab/drafts/CR_genome/resid_het/")
library(ggplot2)


data<-read.csv("PX506_resid_het.counts.txt",sep="\t",header=F)
#Chromosome sizes, C. remanei
cr1=17247545
cr2=19935723
cr3=17877849
cr4=25790997
cr5=22502457
cr6=21501900

CHROMOSOME <-
  data.frame(CHR = rep(c("I","I","II", "II", "III", "III","IV","IV", "V","V","X", "X"),1),
             Boundary = c(cr1 / 4, cr1 * 3/ 4, cr2 / 4, cr2 * 3 / 4, cr3 / 4, cr3 * 3 / 4, cr4 / 4, cr4 * 3 / 4, cr5 / 4, cr5 * 3 / 4, cr6 / 4, cr6 * 3 / 4))

                          #, ce1 / 4, ce1 * 3 / 4, ce2 / 4, ce2 * 3/ 4, ce3 / 4, ce3 * 3 / 4, ce4 / 4, ce4 * 3 / 4, ce5 / 4, ce5 * 3 / 4, ce6 / 4, ce6 * 3 / 4), Species=c(rep("C.remanei",12),rep("C.elegans",12)))

colnames(data)<-c("CHR","START","END","COUNT")


res<-ggplot(data[data$CHR!="MtDNA" ,], aes_string(x = "START", y = "COUNT", color= "CHR", fill="CHR", ordered = FALSE)) + facet_grid(~CHR,scales  ="free") +theme_bw() + labs(title ="", x = "Genome Position (Mb)", y = "Number of polymorphic SNPs (100 Kb)") + theme(axis.text.x = element_text(angle=0,vjust=0.5, hjust = 0.5), plot.title = element_text(hjust = 0.5)) + theme(legend.position = "none") + scale_colour_manual(values = rep("#FF7C01",6),name="") + scale_fill_manual(values = rep("#FF7C01",6),name="") + scale_x_continuous(labels = function(x) x/1000000) +theme(strip.background = element_rect(colour="white", fill="white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + geom_point(size=0.75,alpha=0.75) + theme(strip.text= element_text(size = 10)) + geom_vline(data=CHROMOSOME,aes(xintercept =Boundary), colour="grey50",size=0.5,linetype="dashed")
ggsave("CR_residual_heterozygosity.pdf",res, width=7, height=3, units="in", scale=1.15)
                                            
