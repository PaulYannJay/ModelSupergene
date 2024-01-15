#### Collection of scripts to produce the main plots of the paper####

setwd("~/Paper/SupergeneModel/ReviewMolEcol/Code/")
library(ggplot2)
library(cowplot)
library(ggpubr)
library(tidyverse)
library(RColorBrewer)
library(viridis)
ThemeSobr=  theme(
  panel.border = element_blank(),  
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  text = element_text(size=12),
  axis.line = element_line(colour = "grey"),
  legend.spacing.y= unit(0, 'cm'),
  axis.title = element_text(face="bold"),
  plot.title=element_text(size=12, face="bold",hjust=0.5, vjust=2),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  axis.text = element_text(size = 10)
)


## Data for h=1, u=1e-6 --> Figure 2-3
data2 <- read.table("Data/SensitivityAnalysis2DataFigure2-3.txt",
                          header=T,sep=";", skipNul = T)

source("inputGenotypes.R") ### Source file containing the Genotype
data2$chromosomeType = "0"
colnames(data2)=c("m", "r", "HomoFit", "HetFit", "BadHomoFit","Dominance","rep","state","deme","chromosome"    
,"frequency","chromosomeType")

#### Change haplotype number to haplotype name, based on the source file
data2$chromosomeType[data2$chromosome %in% chromBAlist2] = "BA"
data2$chromosomeType[data2$chromosome %in% chromABlist2] = "AB"
data2$chromosomeType[data2$chromosome %in% chromablist2] = "ab"
data2$chromosomeType[data2$chromosome %in% chrombalist2] = "ba"
data2$chromosomeType[data2$chromosome %in% chromOtherlist2] = "other"
data2$chromosomeType[data2$chromosome %in% chromInvOthlist2] = "InvOther"

### Define as factor
data2$chromosomeType <- factor(data2$chromosomeType, 
                               levels=c("BA","AB","ab","ba", "other", "InvOther"))

### Merge the frequency of recombinant ("other") haplotypes
FreqChromosomeType <- data2 %>% group_by(m,r,HomoFit,HetFit,BadHomoFit,Dominance,rep,state,deme,chromosomeType) %>% summarise(frequency = sum(frequency))

FreqChromosomeTypeForm=FreqChromosomeType
FreqChromosomeTypeForm$allele=paste(FreqChromosomeTypeForm$chromosomeType, FreqChromosomeTypeForm$state, FreqChromosomeTypeForm$deme, sep="_") #Allele_state_deme
FreqChromosomeTypeForm$chromosomeType=NULL
FreqChromosomeTypeForm$state=NULL
FreqChromosomeTypeForm$deme=NULL
FreqChromosomeTypeForm_Wide=spread(FreqChromosomeTypeForm, allele, frequency)

FreqChromosomeTypeForm_Wide$RM0_1=FreqChromosomeTypeForm_Wide$ba_0_1+FreqChromosomeTypeForm_Wide$BA_0_1+FreqChromosomeTypeForm_Wide$InvOther_0_1 #Frequency of inversion (Recombination modifier, RM) in population 1 in phase 0
FreqChromosomeTypeForm_Wide$RM0_2=FreqChromosomeTypeForm_Wide$ba_0_2+FreqChromosomeTypeForm_Wide$BA_0_2+FreqChromosomeTypeForm_Wide$InvOther_0_2
FreqChromosomeTypeForm_Wide$RM1_1=FreqChromosomeTypeForm_Wide$ba_1_1+FreqChromosomeTypeForm_Wide$BA_1_1+FreqChromosomeTypeForm_Wide$InvOther_1_1
FreqChromosomeTypeForm_Wide$RM1_2=FreqChromosomeTypeForm_Wide$ba_1_2+FreqChromosomeTypeForm_Wide$BA_1_2+FreqChromosomeTypeForm_Wide$InvOther_1_2

FreqChromosomeTypeForm_Wide$FitCost=1-FreqChromosomeTypeForm_Wide$HomoFit #The fitess cost of the inversion ("mutation load") is the inverse of the fitness of the inversion homozygotes.


#### Two potential outcomes, we keep only the one where the inversion does not evolve in Pop2 (so when the inversion evolves in pop 1) ###
FreqChromosomeTypeForm_Wide=subset(FreqChromosomeTypeForm_Wide, FreqChromosomeTypeForm_Wide$ba_0_2<0.1)

################################### 
## Figure3 ##
################################### 

FreqChromosomeTypeForm_WideSum=FreqChromosomeTypeForm_Wide %>% group_by(FitCost, r, m) %>% summarise_all(mean)
base=ggplot(FreqChromosomeTypeForm_WideSum[FreqChromosomeTypeForm_WideSum$m==0.1,])
PlotRM0xRecomb0=base+ geom_point(aes(x=RM0_1, y=other_0_1, color=FitCost, shape=as.factor(r)), size=3, alpha=0.8)+
  labs(x="Inversion frequency", y="Frequency of Ab & aB")+
  ggtitle("Recombinant haplotypes (aB-Ab)")+
  scale_shape("Recombination rate", solid = T)+
  scale_color_viridis(name = "Mutation load", option = "A", begin = 0.0, end = 0.8, direction = -1)+ThemeSobr

PlotRM0xSecondPeak=base+geom_point(aes(x=RM0_1, y=ab_0_1, color=FitCost, shape=as.factor(r)), alpha=0.8, size=3)+
  labs(x="Inversion frequency (BA)", y="Frequency of ab")+
  ggtitle("Immigrant haplotype (ab)")+
  theme(axis.title.x = element_text(face="bold"))+
  scale_shape("Recombination rate", solid = T)+
  scale_color_viridis(name = "Mutation load", option = "A", begin = 0.0, end = 0.8, direction = -1)+ThemeSobr

PlotRM0xMainPeak=base+geom_point(aes(x=RM0_1, y=AB_0_1, color=FitCost, shape=as.factor(r)), alpha=0.8, size=3)+
  labs(x="Inversion frequency", y="Frequency of AB")+
  ggtitle("Local haplotype (AB)")+
  scale_shape("Recombination rate", solid = T)+
  scale_color_viridis(name = "Mutation load", option = "A", begin = 0.0, end = 0.8, direction = -1)+ThemeSobr

PlotRM1xRecomb1=base+ geom_point(aes(x=RM1_1, y=other_1_1, color=FitCost, shape=as.factor(r)), size=3, alpha=0.8)+
  labs(x="Inversion frequency", y="Frequency of aB-Ab")+
  ggtitle("Recombinant haplotypes (aB-Ab)")+
  scale_shape("Recombination rate", solid = T)+
  scale_color_viridis(name = "Mutation load", option = "A", begin = 0.0, end = 0.8, direction = -1)+ThemeSobr

PlotRM1xSecondPeak1=base+geom_point(aes(x=RM1_1, y=ab_1_1, color=FitCost, shape=as.factor(r)), alpha=0.8, size=3)+
  labs(x="Inversion frequency (BA)", y="Frequency of ab")+
  ggtitle("Immigrant haplotype (ab)")+
  theme(axis.title.x = element_text(face="bold"))+
  scale_shape("Recombination rate", solid = T)+
  scale_color_viridis(name = "Mutation load", option = "A", begin = 0.0, end = 0.8, direction = -1)+ThemeSobr

PlotRM1xMainPeak1=base+geom_point(aes(x=RM1_1, y=AB_1_1, color=FitCost, shape=as.factor(r)), alpha=0.8, size=3)+
  labs(x="Inversion frequency", y="Frequency of AB")+
  ggtitle("Local haplotype (AB)")+
  scale_shape("Recombination rate", solid = T)+
  scale_color_viridis(name = "Mutation load", option = "A", begin = 0.0, end = 0.8, direction = -1)+ThemeSobr

legend=cowplot::get_legend(PlotRM0xMainPeak +
                             theme(legend.position = "top",
                                   legend.direction = "horizontal",
                                   legend.justification="center" ,
                                   legend.box.just = "bottom",
                                   legend.text = element_text(size=8)))

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
PlotRM1xRecomb1_SSleg=PlotRM1xRecomb1 + my_themeTopPlot
PlotRM1xMainPeak1_SSleg=PlotRM1xMainPeak1 + my_themeTopPlot
PlotRM1xSecondPeak1_SSleg=PlotRM1xSecondPeak1 + my_themeBotPlot
Title1=ggdraw() + draw_label("Migratory phase", angle=0, fontface="bold", size=15, vjust=1)
Title2=ggdraw() + draw_label("Non-migratory phase", angle=0, fontface="bold", size=15, vjust=1)
plotsABC <- align_plots(PlotRM0xMainPeak_SSleg, PlotRM0xRecomb0_SSleg, 
                       PlotRM0xSecondPeak_SSleg,align = 'hv', axis = 'rltp')
plotsDEF <- align_plots(PlotRM1xMainPeak1_SSleg, PlotRM1xRecomb1_SSleg, 
                        PlotRM1xSecondPeak1_SSleg,align = 'hv', axis = 'rltp')

gg_all_Inv_Mig = plot_grid(Title1,plotsABC[[1]],plotsABC[[2]],plotsABC[[3]],labels=c('', 'A', 'B', 'C'), ncol=1, rel_heights = c(0.1,1,1,1.36))
gg_all_Inv_NoMig = plot_grid(Title2, plotsDEF[[1]], plotsDEF[[2]], plotsDEF[[3]],labels=c('','D', 'E', 'F'), ncol=1, rel_heights = c(0.1,1,1,1.36))

labelAx=ggdraw() + draw_label("Frequency", angle=90, fontface="bold", size=15, vjust=1)
gg_all=plot_grid(labelAx,gg_all_Inv_Mig,gg_all_Inv_NoMig,ncol=3, rel_widths = c(0.20,2,2)) 
Figure3=plot_grid( legend, gg_all, ncol=1, rel_heights = c(0.3,4))
Figure3
save_plot("~/Paper/SupergeneModel/ReviewMolEcol/Code/Figure3_V3.png", Figure3, nrow = 3, base_aspect_ratio = 2.5)
save_plot("~/Paper/SupergeneModel/ReviewMolEcol/Code/Figure3_V3.svg", Figure3, nrow = 3, ncol=2)

################################### 
#### Scenario frequency. Figure 2A ###
################################### 

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
FreqChromosomeTypeForm_Wide_SubSub=FreqChromosomeTypeForm_Wide_Sub[FreqChromosomeTypeForm_Wide_Sub$HomoFit %in% c(0.0,0.1,0.20,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),]
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
  scale_x_continuous(breaks=pretty(FreqChromosomeTypeForm_Wide_SubSub$HomoFit, 5))+
  scale_fill_manual(values=c("#8960b3","#61ad65","#b9495e","#b68b39","#0000ff", "#5555ff", "#9999ff" , "#9ab9ff"))

Figure2 

################################### 
##### Simulation plots. (Figure 2B-G) ##
################################### 
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

