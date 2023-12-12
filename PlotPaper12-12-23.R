#### Plot Papier Model ####

setwd("~/Paper/SupergeneModel/ReviewMolEcol/Code/")
library(ggplot2)
library(cowplot)
library(tidyverse)
library(RColorBrewer)

#####
## Data for h=1, u=1e-6 --> Figure 2-3
data2 <- read.table("Data/SensitivityAnalysis2DataFigure2-3.txt",
                          header=T,sep=";", skipNul = T)

source("inputGenotypes.R") ### Source file containing the Genotype
data2$chromosomeType = "0"
colnames(data2)=c("m", "r", "MutLoad", "HetFit", "BadHomoFit","Dominance","rep","state","deme","chromosome"    
,"frequency","chromosomeType")

#### Change Genotype number to Genotype name, based on the source file
data2$chromosomeType[data2$chromosome %in% chromBAlist2] = "BA"
data2$chromosomeType[data2$chromosome %in% chromABlist2] = "AB"
data2$chromosomeType[data2$chromosome %in% chromablist2] = "ab"
data2$chromosomeType[data2$chromosome %in% chrombalist2] = "ba"
data2$chromosomeType[data2$chromosome %in% chromOtherlist2] = "other"
data2$chromosomeType[data2$chromosome %in% chromInvOthlist2] = "InvOther"

### Define as factor
data2$chromosomeType <- factor(data2$chromosomeType, 
                               levels=c("BA","AB","ab","ba", "other", "InvOther"))

### Merge the frequency of recombinant phenotype
FreqChromosomeType <- data2 %>% group_by(m,r,MutLoad,HetFit,BadHomoFit,Dominance,rep,state,deme,chromosomeType) %>% summarise(frequency = sum(frequency))

FreqChromosomeTypeForm=FreqChromosomeType
FreqChromosomeTypeForm$allele=paste(FreqChromosomeTypeForm$chromosomeType, FreqChromosomeTypeForm$state, FreqChromosomeTypeForm$deme, sep="_") #Allele_state_deme
FreqChromosomeTypeForm$chromosomeType=NULL
FreqChromosomeTypeForm$state=NULL
FreqChromosomeTypeForm$deme=NULL
FreqChromosomeTypeForm_Wide=spread(FreqChromosomeTypeForm, allele, frequency)
FreqChromosomeTypeForm_Wide$RM0_1=FreqChromosomeTypeForm_Wide$ba_0_1+FreqChromosomeTypeForm_Wide$BA_0_1+FreqChromosomeTypeForm_Wide$InvOther_0_1
FreqChromosomeTypeForm_Wide$RM0_2=FreqChromosomeTypeForm_Wide$ba_0_2+FreqChromosomeTypeForm_Wide$BA_0_2+FreqChromosomeTypeForm_Wide$InvOther_0_2
FreqChromosomeTypeForm_Wide$RM1_1=FreqChromosomeTypeForm_Wide$ba_1_1+FreqChromosomeTypeForm_Wide$BA_1_1+FreqChromosomeTypeForm_Wide$InvOther_1_1
FreqChromosomeTypeForm_Wide$RM1_2=FreqChromosomeTypeForm_Wide$ba_1_2+FreqChromosomeTypeForm_Wide$BA_1_2+FreqChromosomeTypeForm_Wide$InvOther_1_2
FreqChromosomeTypeForm_Wide$FitCost=1-FreqChromosomeTypeForm_Wide$MutLoad
my_palette = brewer.pal(n = 9, "Blues")[4:9]

#### Two potential outcome, we keep only the one where the inversion does not evolve in Pop2 ###
FreqChromosomeTypeForm_Wide=subset(FreqChromosomeTypeForm_Wide, FreqChromosomeTypeForm_Wide$ba_0_2<0.1)

## Figure3 ##
base=ggplot(FreqChromosomeTypeForm_Wide[FreqChromosomeTypeForm_Wide$m==0.1,])

PlotRM0xRecomb0=base+ geom_point(aes(x=RM0_1, y=other_0_1, color=as.factor(r)))+
  labs(x="Inversion frequency", y="Frequency of valley phenotype\n (recombinant) in pop 1")+
  ggtitle("Recombinant haplotype")+
  scale_color_manual(name = "Recombination rate", values = my_palette)

PlotRM0xSecondPeak=base+geom_point(aes(x=RM0_1, y=ab_0_1, color=as.factor(r)), alpha=0.5)+
  labs(x="Inversion frequency", y="Frequency of ancestral \n second haplotypes")+
  ggtitle("Immigrant haplotype")+
  theme(axis.title.x = element_text(face="bold"))+
  scale_color_manual(name = "Recombination rate", values = my_palette)

PlotRM0xMainPeak=base+geom_point(aes(x=RM0_1, y=AB_0_1, color=as.factor(r)), alpha=0.7)+
  labs(x="Inversion frequency", y="Frequency of ancestral \n main haplotype")+
  ggtitle("Local haplotype")+
  scale_color_manual(name = "Recombination rate", values = my_palette)

legend=get_legend(PlotRM0xMainPeak+     guides(color = guide_legend(nrow = 1)) +
                    theme(legend.position = "top",
                          legend.direction = "horizontal",
                          legend.justification="center" ,
                          legend.box.just = "bottom"))

my_themeTopPlot=theme(legend.position="none", 
                      axis.title=element_blank(),
                      axis.text.x=element_blank(),
                      plot.title = element_text(color="grey35", face="bold"),
                      plot.margin = unit(c(5, 5, 8, 5), "pt"))

my_themeBotPlot=theme(legend.position="none", 
                      axis.title.y=element_blank(),
                      plot.title = element_text(color="grey35", face="bold"))

PlotRM0xRecomb0_SSleg=PlotRM0xRecomb0 + my_themeTopPlot
PlotRM0xMainPeak_SSleg=PlotRM0xMainPeak + my_themeTopPlot
PlotRM0xSecondPeak_SSleg=PlotRM0xSecondPeak + my_themeBotPlot

gg_all_Inv = plot_grid(PlotRM0xMainPeak_SSleg, PlotRM0xRecomb0_SSleg, 
                       PlotRM0xSecondPeak_SSleg,labels=c('A', 'B', 'C'), ncol=1, rel_heights = c(1,1,1.36))

labelAx=ggdraw() + draw_label("Frequency", angle=90, fontface="bold", size=15, vjust=1)
gg_all=plot_grid(labelAx,gg_all_Inv,ncol=2, rel_widths = c(0.20,2)) 
Figure3=plot_grid( legend, gg_all, ncol=1, rel_heights = c(0.3,4))
Figure3


#### Freqency Scenario. Figure 2A ###
dom=1.0
FreqChromosomeTypeForm_Wide_Sub=subset(FreqChromosomeTypeForm_Wide, (    FreqChromosomeTypeForm_Wide$HetFit==0.2 & #Change the value of the parameter to produce the different figures.
                                                                           FreqChromosomeTypeForm_Wide$BadHomoFit==0.5 &
                                                                           FreqChromosomeTypeForm_Wide$m==0.1 &
                                                                           FreqChromosomeTypeForm_Wide$Dominance==dom))
FreqChromosomeTypeForm_Wide_Sub$scenar=NA
for (i in 1:nrow(FreqChromosomeTypeForm_Wide_Sub))
{
  if (FreqChromosomeTypeForm_Wide_Sub$BA_0_1[i] < 0.05 & FreqChromosomeTypeForm_Wide_Sub$ba_0_2[i] < 0.05 )
  {
    FreqChromosomeTypeForm_Wide_Sub$scenar[i]="No inversion evolved"
  }
  else if (FreqChromosomeTypeForm_Wide_Sub$ba_1_2[i] > 0.95)
  {
    FreqChromosomeTypeForm_Wide_Sub$scenar[i]="Inversion fixed in pop 2"
  }
  else if (FreqChromosomeTypeForm_Wide_Sub$BA_1_1[i]>0.95)
  {
    FreqChromosomeTypeForm_Wide_Sub$scenar[i]="Inversion fixed in pop 1"
  }
  else if (FreqChromosomeTypeForm_Wide_Sub$ba_0_2[i]>0.05 & FreqChromosomeTypeForm_Wide_Sub$BA_0_1[i]>0.05 & FreqChromosomeTypeForm_Wide_Sub$BA_1_1[i]<0.05 & FreqChromosomeTypeForm_Wide_Sub$ba_1_2[i] < 0.05)
  {
    FreqChromosomeTypeForm_Wide_Sub$scenar[i]= "Inversion evolved on both haplotype during phase M \n but disapear during phase NM"
  }
  else if (FreqChromosomeTypeForm_Wide_Sub$ba_0_2[i]>0.05 & FreqChromosomeTypeForm_Wide_Sub$BA_1_1[i]<0.05 & FreqChromosomeTypeForm_Wide_Sub$ba_1_2[i] < 0.05)
  {
    FreqChromosomeTypeForm_Wide_Sub$scenar[i]= "Inversion evolved on recessive haplotype during phase M \n but disapear during phase NM"
  }
  else if (FreqChromosomeTypeForm_Wide_Sub$BA_0_1[i]>0.05 & FreqChromosomeTypeForm_Wide_Sub$BA_1_1[i]<0.05 & FreqChromosomeTypeForm_Wide_Sub$ba_1_2[i] < 0.05)
  {
    FreqChromosomeTypeForm_Wide_Sub$scenar[i]= "Inversion evolved on dominant haplotype during phase M \n but disapear during phase NM"
  }
  else if (FreqChromosomeTypeForm_Wide_Sub$BA_1_1[i]<0.95 & FreqChromosomeTypeForm_Wide_Sub$BA_1_1[i]>0.05 )
  {
    FreqChromosomeTypeForm_Wide_Sub$scenar[i]="Supergene in pop 1"
  }
  else if (FreqChromosomeTypeForm_Wide_Sub$ba_1_2[i]<0.95 & FreqChromosomeTypeForm_Wide_Sub$ba_1_2[i]>0.05 )
  {
    FreqChromosomeTypeForm_Wide_Sub$scenar[i]="Supergene in pop 2"
  }
  else
  {
    FreqChromosomeTypeForm_Wide_Sub$scenar[i]="blabla"
  }
}

FreqChromosomeTypeForm_Wide_Sub$scenar=as.factor(FreqChromosomeTypeForm_Wide_Sub$scenar)
FreqChromosomeTypeForm_Wide_Sub$FitCost=as.numeric(FreqChromosomeTypeForm_Wide_Sub$FitCost)
FreqChromosomeTypeForm_Wide_SubSub=FreqChromosomeTypeForm_Wide_Sub[FreqChromosomeTypeForm_Wide_Sub$MutLoad %in% c(0.0,0.1,0.20,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),]
FreqChromosomeTypeForm_Wide_SubSub=FreqChromosomeTypeForm_Wide_Sub

Figure2=ggplot(FreqChromosomeTypeForm_Wide_SubSub,
                      aes(x=FitCost,
                          fill=scenar))+
  geom_histogram(position= "fill", binwidth = 0.05, colour="black")+  #Change binwidth for different precision.
  scale_y_continuous(labels = scales::percent)+
  labs(y="Scenario Frequency", x='MutationLoad')+
  facet_grid(.~r)+
  
  ggtitle(label = "Recombination Rate")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.title=element_blank())+
  scale_x_continuous(breaks=pretty(FreqChromosomeTypeForm_Wide_SubSub$MutLoad, 5))+
  scale_fill_manual(values=c("#8960b3","#61ad65","#b9495e","#b68b39","#0000ff", "#5555ff", "#9999ff" , "#9ab9ff"))

Figure2 # Figure 2 and S?

#Individual simulation plot. (Figure 2B-G) ##
data2 <- read.table("Data/SingleTimeSeries2SimulorangeDBlue.20.txt",header=T,sep=";")
data2 <- data2[data2$time %% 10==0,]
data2$frequency <- as.numeric(gsub("f", "e", as.character(data2$frequency)))
data2$chromosomeType = "0"
data2$chromosomeType[data2$chromosome %in% chromBAlist2] = "BA"
data2$chromosomeType[data2$chromosome %in% chrombalist2] = "ba"
data2$chromosomeType[data2$chromosome %in% chromInvOthlist2] = "Inverted recombinants"
data2$chromosomeType[data2$chromosome %in% chromABlist2] = "AB"
data2$chromosomeType[data2$chromosome %in% chromablist2] = "ab"
data2$chromosomeType[data2$chromosome %in% chromOtherlist2] = "Recombinants"

data2$chromosomeType <- factor(data2$chromosomeType, 
                               levels=c("BA","ba", "Inverted recombinants", "AB","ab","Recombinants")) ##ModifPJ

FreqChromosomeType <- data2 %>% group_by(time,deme,chromosomeType) %>% summarise(frequency = sum(frequency))


PlotSimulation = ggplot(data=FreqChromosomeType)+
  geom_line(aes(time, frequency, color=chromosomeType,lty=chromosomeType))+
  scale_color_manual(name="Allele", 
                     breaks=c("BA","ba", "Inverted recombinants", "AB","ab","Recombinants"),
                     values=c("red","green","grey50","orange","blue","grey"))+
  scale_linetype_manual(name="Allele", 
                        breaks=c("BA","ba", "Inverted recombinants", "AB","ab","Recombinants"),
                        values=c(1,1,1,1,1,1))+
  ylab("Frequency")+xlab("Generation")+
  ylim(-1e-4,1+1e-4)+
  theme(legend.position="bottom")+
  
  facet_wrap(~deme)
PlotSimulation

