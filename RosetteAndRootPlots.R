#Generating plots for rosette area and root length assays

library(rstatix)
library(ggpubr)
library(dplyr)

#set the colour palette for plots  
mypalette = c("#0073C2FF", "#00A36C", "#868686FF","#FFD700","#00FFFF","#009999")


#Rosette Data Plotting

#read in data 

leafAreaTable<-read.csv("~/PlantPhenotyping/RosetteArea/areaOfLeaves.csv")

ggboxplot(leafAreaTable, x = "label", y = "Area", color = "label", palette = mypalette,  outlier.shape = NA) +  
  geom_jitter(width = 0.2,aes(color = label)) +  ylab("Rosette Area (cm^2)") + scale_y_continuous(limits = c(0,NA)) +
  stat_compare_means()   +   # Add global p-value. Function does this automatically. label.y specification is positioning relative to axes
  theme(text = element_text(size = 10))


#Root length Data Plotting

#read in data 

rootLengthTable<-read.csv("~/PlantPhenotyping/RootLength2.csv")


kt<- kruskal.test(length..cm. ~ Genotype, data = rootLengthTable)


if(kt$p.value < 0.05){
  #if kruskal wallis test is significant to 0.05, perform and plot the pairwise post-hoc comparison 
  

  rootLengthTable %>% dunn_test(length..cm. ~ Genotype,p.adjust.method="fdr") 
  #pwc = pairwise comparison
  pwc<- rootLengthTable %>% dunn_test(length..cm. ~ Genotype,p.adjust.method="fdr") 
  #position over appropriate banners
  pwc <- pwc %>% add_xy_position(x = "Genotype")
  #round it
  pwc$p.adj<-round(pwc$p.adj, digits = 3)

ggboxplot(rootLengthTable, x = "Genotype", y = "length..cm.", color = "Genotype", palette = mypalette,  outlier.shape = NA) +  
  geom_jitter(width = 0.2,aes(color = Genotype)) +  ylab("Root Length (cm)") + scale_y_continuous(limits = c(0,NA)) +
  stat_compare_means(label.y= max(rootLengthTable$length..cm.)+(max(rootLengthTable$length..cm.)))   +   # Add global p-value. Function does this automatically. label.y specification is positioning relative to axes
  stat_pvalue_manual(pwc, label= "p = {p.adj}",step.increase = 0.09) + # Add pairwise comparisons p-value change height of label by vjust = -0.2
  theme(text = element_text(size = 10))
}

rootLengthTable %>% count(Genotype)
#Just AOX plates
justAOX <- rootLengthTable %>% filter(Type == 1)
#Just IVD plates
justIVD <- rootLengthTable %>% filter(Type == 2)
#Just mCherry plates
justmCherry <- rootLengthTable %>% filter(Type == 3)
#Just Col-0 background plates
justCol0bck1 <- rootLengthTable %>% filter(Type == 4)
justCol0bck2 <- rootLengthTable %>% filter(Type == 5)
justCol0bck<-rbind(justCol0bck1,justCol0bck2)

ggboxplot(justAOX, x = "Genotype", y = "length..cm.", color = "Genotype", palette = mypalette,  outlier.shape = NA) +  
  geom_jitter(width = 0.2,aes(color = Genotype)) +  ylab("Root Length (cm)") + scale_y_continuous(limits = c(0,NA)) +
  stat_compare_means()   +   # Add global p-value. Function does this automatically. label.y specification is positioning relative to axes
  theme(text = element_text(size = 10))

ggboxplot(justIVD, x = "Genotype", y = "length..cm.", color = "Genotype", palette = mypalette,  outlier.shape = NA) +  
  geom_jitter(width = 0.2,aes(color = Genotype)) +  ylab("Root Length (cm)") + scale_y_continuous(limits = c(0,NA)) +
  stat_compare_means()   +   # Add global p-value. Function does this automatically. label.y specification is positioning relative to axes
  theme(text = element_text(size = 10))

ggboxplot(justmCherry, x = "Genotype", y = "length..cm.", color = "Genotype", palette = mypalette,  outlier.shape = NA) +  
  geom_jitter(width = 0.2,aes(color = Genotype)) +  ylab("Root Length (cm)") + scale_y_continuous(limits = c(0,NA)) +
  stat_compare_means()   +   # Add global p-value. Function does this automatically. label.y specification is positioning relative to axes
  theme(text = element_text(size = 10))

ggboxplot(justCol0bck, x = "Genotype", y = "length..cm.", color = "Genotype", palette = mypalette,  outlier.shape = NA) +  
  geom_jitter(width = 0.2,aes(color = Genotype)) +  ylab("Root Length (cm)") + scale_y_continuous(limits = c(0,NA)) +
  stat_compare_means()   +   # Add global p-value. Function does this automatically. label.y specification is positioning relative to axes
  theme(text = element_text(size = 10))
