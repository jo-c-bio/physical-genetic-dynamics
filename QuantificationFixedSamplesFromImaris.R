#mtDNA quantification from fixed cells
#This data comes exported from Imaris- I chose the strictest and easiest definition of 'mtDNA presence'
#Using the pre-set filter of stats-> detailed-> 'nucleus number of vesicles VesicleType=vesicles mtDNA', for each sample

#Data Input

########################################################################
#New Samples - just epidermal so far, 11 x cells. 
#The original file names had stupid labels to prevent Imaris reading them as stacked time series files!
########################################################################

library(dplyr)
setwd("/Users/joannachustecki/Documents/PostDoc23-Data/nucleoidQuantification/round2")

#epidermal cell vols
for(i in c(1:6,8:12)){
  #skip=3 because headers are found on row 4 of each data sheet
  x <- read.csv(paste("e",i,"_Statistics/e",i,"_Cell_Volume.csv",sep=""),skip=3)
  assign(paste("epi",i,"v",sep=""), x )
}

#epidermal cell number of mitos with mtDNA, read in volume data alongside.
for(i in c(1:6,8:12)){
  x <- read.csv(paste("e",i,"_Statistics/e",i,"_Nucleus_Number_of_Vesicles_VesicleType=Vesicles_mtDNA.csv",sep=""),skip=3)
  y <- read.csv(paste("e",i,"_Statistics/e",i,"_Nucleus_Volume.csv",sep=""),skip=3)
  colnames(y) <- paste(colnames(y), ".V", sep="")
  z<-cbind(x,y)
  assign(paste("epi",i,sep=""), z )
}

#epidermal cell volume of mitos
for(i in c(1:6,8:12)){
  x <- read.csv(paste("e",i,"_Statistics/e",i,"_Nucleus_Volume.csv",sep=""),skip=3)
  assign(paste("epi",i,"mv",sep=""), x )
}

#Mean Intensity nucleus (mitos) in channel 2 (mito)
for(i in c(1:6,8:12)){
  x <- read.csv(paste("e",i,"_Statistics/e",i,"_Nucleus_Intensity_Mean_Ch=2_Img=1.csv",sep=""),skip=3)
  assign(paste("epi",i,"mI",sep=""), x )
}

#Mean Intensity vesicles (SYBR) channel 1
for(i in c(1:6,8:12)){
  x <- read.csv(paste("e",i,"_Statistics/e",i,"_Vesicles_mtDNA_Intensity_Mean_Ch=1_Img=1.csv",sep=""),skip=3)
  assign(paste("epi",i,"vI",sep=""), x )
}



##############################################
#Then we can collect the hypocotyl data
##############################################

setwd("/Users/joannachustecki/Documents/PostDoc23-Data/nucleoidQuantification/round2/")

#might need to remove 5 and 7- they have a much higher resolution, 
#10.7501 pixels per micron compared to 4.8272 for the rest. 
#hypocotyl cell vols
for(i in c(1,2,3,4,5,7,8,11)){
  x <- read.csv(paste("h",i,"_Statistics/h",i,"_Cell_Volume.csv",sep=""),skip=3)
  assign(paste("hyp",i,"v",sep=""), x )
}
hypCol0v<- read.csv(paste("controlCol_Statistics/controlCol_Cell_Volume.csv",sep=""),skip=3)
hypmCherryv<- read.csv(paste("controlmCherry_Statistics/controlmCherry_Cell_Volume.csv",sep=""),skip=3)

#hypocotyl cell number of mitos with mtDNA, read in volume data alongside. 
for(i in c(1,2,3,4,5,7,8,11)){
  x <- read.csv(paste("h",i,"_Statistics/h",i,"_Nucleus_Number_of_Vesicles_VesicleType=Vesicles_mtDNA.csv",sep=""),skip=3)
  y <- read.csv(paste("h",i,"_Statistics/h",i,"_Nucleus_Volume.csv",sep=""),skip=3)
  colnames(y) <- paste(colnames(y), ".V", sep="")
  z<-cbind(x,y)
  assign(paste("hyp",i,sep=""), z )
}

#hypocotyl cell  volume of mitos
for(i in c(1,2,3,4,5,7,8,11)){
  x <- read.csv(paste("h",i,"_Statistics/h",i,"_Nucleus_Volume.csv",sep=""),skip=3)
  assign(paste("hyp",i,"mv",sep=""), x )
}
hypmCherrymv<- read.csv(paste("controlmCherry_Statistics/controlmCherry_Nucleus_Volume.csv",sep=""),skip=3)

#Mean Intensity nucleaus (mitos) in channel 2 (mito)
for(i in c(1,2,3,4,5,7,8,11)){
  x <- read.csv(paste("h",i,"_Statistics/h",i,"_Nucleus_Intensity_Mean_Ch=2_Img=1.csv",sep=""),skip=3)
  assign(paste("hyp",i,"mI",sep=""), x )
}
hypmCherrymI<- read.csv(paste("controlmCherry_Statistics/controlmCherry_Nucleus_Intensity_Mean_Ch=2_Img=1.csv",sep=""),skip=3)


#Mean Intensity vesicles (SYBR) channel 1
for(i in c(1,2,3,4,5,7,8,11)){
  x <- read.csv(paste("h",i,"_Statistics/h",i,"_Vesicles_mtDNA_Intensity_Mean_Ch=1_Img=1.csv",sep=""),skip=3)
  assign(paste("hyp",i,"vI",sep=""), x )
}
hypCol0vI<- read.csv(paste("controlCol_Statistics/controlCol_Vesicles_mtDNA_Intensity_Mean_Ch=1_Img=1.csv",sep=""),skip=3)


realSampleNames <- c("epi1","epi2","epi3","epi4","epi5","epi6","epi8","epi9","epi10","epi11","epi12",
                 "hyp1","hyp2","hyp3","hyp4","hyp5","hyp7","hyp8","hyp11")
graphSampleNames <- c("epi1","epi2","epi3","epi4","epi5","epi6","epi7","epi8","epi9","epi10","epi11",
                     "hyp1","hyp2","hyp3","hyp4","hyp5","hyp6","hyp7","hyp8")
numbermitoslist<- list(epi1,epi2,epi3,epi4,epi5,epi6,epi8,epi9,epi10,epi11,epi12,
                       hyp1,hyp2,hyp3,hyp4,hyp5,hyp7,hyp8,hyp11)
volcellslist<- list(epi1v,epi2v,epi3v,epi4v,epi5v,epi6v,epi8v,epi9v,epi10v,epi11v,epi12v,
                    hyp1v,hyp2v,hyp3v,hyp4v,hyp5v,hyp7v,hyp8v,hyp11v)
mitovolslist<- list(epi1mv,epi2mv,epi3mv,epi4mv,epi5mv,epi6mv,epi8mv,epi9mv,epi10mv,epi11mv,epi12mv,
                    hyp1mv,hyp2mv,hyp3mv,hyp4mv,hyp5mv,hyp7mv,hyp8mv,hyp11mv)
mitoIntensitylist<- list(epi1mI,epi2mI,epi3mI,epi4mI,epi5mI,epi6mI,epi8mI,epi9mI,epi10mI,epi11mI,epi12mI,
                    hyp1mI,hyp2mI,hyp3mI,hyp4mI,hyp5mI,hyp7mI,hyp8mI,hyp11mI)
vesicleIntensitylist<- list(epi1vI,epi2vI,epi3vI,epi4vI,epi5vI,epi6vI,epi8vI,epi9vI,epi10vI,epi11vI,epi12vI,
                   hyp1vI,hyp2vI,hyp3vI,hyp4vI,hyp5vI,hyp7vI,hyp8vI,hyp11vI)

#lets get number of mitochondria with mtDNA present, absent and in total, as well as cell volume

#########.   For gathering the cell volumes associated with each sample ################
#for those with two variables, we want the biggest ones. The second cell in the surface will be an offshoot.

volChecker <- function(vol.df) {
  if (nrow(vol.df) >= 2) {
    return(max(vol.df$Cell.Volume))
  } else{
    return(as.numeric(vol.df$Cell.Volume))
  }
}

### mitochondria volume filter #######
#Imaris gave out thresholded mitochondria with tiny volumes, that do not accurately represent whole organelles. 
#We will here, filter out any mito with volume < 0.1 µm^3, which we wouldn't expect or an organelle with an average width of 1µm

for(i in 1:length(realSampleNames)){
  numbermitoslist[[i]] <- numbermitoslist[[i]] %>% filter(Nucleus.Volume.V > 0.1)
}

mtDNAPresence<-c()
mtDNAAbsence<-c()
mitoTotal<-c()
mitoProportion <- c()
cellVols <- c()

for(i in 1:length(realSampleNames)){
  mitosNwDNApresence <- sum(numbermitoslist[[i]]$Nucleus.Number.of.Vesicles > 0)
  mtDNAPresence <- c(mtDNAPresence,mitosNwDNApresence)
  
  mitosNwDNAabsent <- sum(numbermitoslist[[i]]$Nucleus.Number.of.Vesicles == 0)
  mtDNAAbsence <- c(mtDNAAbsence,mitosNwDNAabsent)
  
  allMitos <- length(numbermitoslist[[i]]$Nucleus.Number.of.Vesicles)
  mitoTotal <- c(mitoTotal,allMitos)
  
  presenceOverTotal <- round(sum(numbermitoslist[[i]]$Nucleus.Number.of.Vesicles > 0)/length(numbermitoslist[[i]]$Nucleus.Number.of.Vesicles),digits=3)
  mitoProportion<- c(mitoProportion,presenceOverTotal)
  
  CellVolume <- volChecker(volcellslist[[i]])
  cellVols <- c(cellVols,CellVolume)
}

df2<-data.frame(c(rep("mtDNAPresence",length(realSampleNames)),rep("mtDNAAbsence",length(realSampleNames))),
               c(mtDNAPresence, mtDNAAbsence),
               c(realSampleNames,realSampleNames),
               c(graphSampleNames,graphSampleNames),
               c(mitoTotal,mitoTotal),
               c(mitoProportion,mitoProportion),
               c(cellVols,cellVols))

colnames(df2)<-c("Type","NumberMitos","realSampleName","graphSampleNames","mitoTotal","mitoProportion","cellVolumeµm3")

#Export dataframe
write.csv(df2,"ProportionalValuesformtDNApresenceAbsenceFromImaris.csv")
#create a stacked barplot with the number of mitos and the number that have mtDNA
library(ggplot2)

df2$graphSampleNames <- as.character(df2$graphSampleNames )
#Then turn it back into a factor with the levels in the correct order
df2$graphSampleNames  <- factor(df2$graphSampleNames, levels=unique(df2$graphSampleNames ) )

mainFigure<-
ggplot(data=df2, aes(x=graphSampleNames, y=NumberMitos, fill=Type)) +
  geom_bar(stat="identity")+
  #add in a value as a proportion here
  geom_text(aes(y=mitoTotal, label=mitoProportion), vjust=-0.35)+
  
  #add in cell volume here, rotated by 90
  geom_text(aes(y=0, label=round(cellVolumeµm3)),hjust=-0.1,size = 10/.pt,angle = 90)+
  
  #add a title
  ggtitle(paste("Upper value = Proportion of mitochondria with mtDNA present",
                "\n",
                "Lower value = Cell volume in µm^3",sep="")) +
  #axis labels
  xlab("Sample name") +
  
  #axis text size increase
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  labs(y= "Number of mitochondria") 

ggsave(
  paste("presenceabsenceTotals.pdf", sep = ""),
  mainFigure,
  width = 25,
  height = 12,
  units = "cm"
)

#Let's do the same plot but re-ordered for the cell volume high-> low

ggplot(data=df2, aes(x=reorder(graphSampleNames, -cellVolumeµm3), y=NumberMitos, fill=Type)) +
  geom_bar(stat="identity")+
  #add in a value as a proportion here
  geom_text(aes(y=mitoTotal, label=mitoProportion), vjust=-0.35)+
  
  #add in cell volume here, rotated by 90
  geom_text(aes(y=0, label=round(cellVolumeµm3)),hjust=-0.1,size = 10/.pt,angle = 90)+
  
  #add a title
  ggtitle(paste("Upper value = Proportion of mitochondria with mtDNA present",
                "\n",
                "Lower value = Cell volume in µm^3",sep="")) +
  #axis labels
   xlab("Sample name") +
 
  #axis text size increase
   theme(axis.text=element_text(size=12),
      axis.title=element_text(size=14))+
  labs(y= "Number of mitochondria")
#pdf export 10x4.5
#let's try some different stats
#mitochondrial volume versus number of mitos and cell volume

mitosVolume<-c()
mitoTotal<-c()
mitoProportion <- c()
cellVols <- c()

for(i in 1:length(realSampleNames)){
  mitosVol <- sum(mitovolslist[[i]]$Nucleus.Volume)
  mitosVolume <- c(mitosVolume,mitosVol)
  
  allMitos <- length(numbermitoslist[[i]]$Nucleus.Number.of.Vesicles)
  mitoTotal <- c(mitoTotal,allMitos)
  
  CellVolume <- volChecker(volcellslist[[i]])
  cellVols <- c(cellVols,CellVolume)
  
  mitoVolOverCellVol <- round(mitosVol/CellVolume,digits=3)
  mitoProportion<- c(mitoProportion,mitoVolOverCellVol)
}

df3<-data.frame(mitosVolume,
                realSampleNames,
                graphSampleNames,
                mitoTotal,
                mitoProportion,
                cellVols)

colnames(df3)<-c("MitoVolSum","realSampleNames","graphSampleNames","mitoTotal","mitoProportion","cellVolumeµm3")


ggplot(data=df3, aes(x=graphSampleNames, y=mitosVolume, fill="lightblue")) +
  geom_bar(stat="identity")+
  #add in a value as a proportion here
  geom_text(aes(y=MitoVolSum, label=mitoProportion), vjust=-0.35)+
  
  #add in cell volume here
  geom_text(aes(y=0, label=round(cellVolumeµm3)),vjust=-0.5,size = 7/.pt)+
  
  #add a title
  ggtitle(paste("upper value = Proportion of mitochondria as a volume to whole cell",
                "\n",
                "lower value = Cell volume in µm^3",sep=""))




#let's just do a simple correlation plot between area and mito no. 

substr(df2$graphSampleNames, 4,5)
substr(df2$graphSampleNames, 1,3)

library("ggpubr")
#using the Spearman correlation as can't guarantee these values are normally distributed

# Shapiro-Wilk normality test 
shapiro.test(df3$mitoTotal)
shapiro.test(df3$cellVolumeµm3)
#QQPlot for normality
ggqqplot(df3$mitoTotal)
ggqqplot(df3$cellVolumeµm3)

df3$g_label<-substr(df3$graphSampleNames, 1,3)
df3$g_number<-substr(df3$graphSampleNames,4,5)

volNumberCorrelation<-
        ggscatter(df3, x = "mitoTotal", y = "cellVolumeµm3", 
          add = "reg.line", 
          fill="g_label",
          palette = c("blue", "orange"),
          label = "g_number",
          shape = 21,
          cor.coef = TRUE, cor.method = "spearman",
          main="Spearman Correlation",
          xlab="Total number of mitochondria") +
          labs(y=expression(paste("Cell volume ", µm^3))) 

ggsave(
  paste("volNumberCorrelation.pdf", sep = ""),
  volNumberCorrelation,
  width = 18,
  height = 12,
  units = "cm"
)

#### Supplementary figures


#mitochondrial volumes, with mCherry control

mitovolslistMV<- append(mitovolslist, list(hypmCherrymv))
graphsampleNamesMV<- c(graphSampleNames,"hypmCherry")


mitosVols<-c()
typeVollist<-c()


for(i in 1:length(graphsampleNamesMV)){
  mitosVol<- mitovolslistMV[[i]]$Nucleus.Volume
  mitosVols<-c(mitosVols,mitosVol)
  
  typeVol<- rep(graphsampleNamesMV[i], length(mitosVol))
  typeVollist<-c(typeVollist,typeVol)
  
}

df4<-data.frame(mitosVols,
                typeVollist,
                log(mitosVols)
)

ggboxplot(df4,
          x = "typeVollist",
          y = "mitosVols",
          color = "typeVollist", outlier.shape = NA) +
  #include outlier shape= NA to remove double-plotted points
  geom_jitter(width = 0.2, aes(color = typeVollist), alpha = 0.3 ) +
  theme(text = element_text(size = 10), legend.position="none")  + labs(y= "Mitochondrial Volume (µm^3)", x=NA)

mitoVolGraph<-
  ggboxplot(df4,
          x = "typeVollist",
          y = "log.mitosVols.",
          color = "typeVollist", outlier.shape = NA) +
  #include outlier shape= NA to remove double-plotted points
  geom_jitter(width = 0.2, aes(color = typeVollist), alpha = 0.3 ) +
  theme(text = element_text(size = 10), legend.position="none",axis.title.x=element_blank())  + labs(y= "Log Mitochondrial Volume (µm^3)")
#weird- hyp 5 and 6 have a higher (double the) resolution than all the other files.
#hyp1 has a sort function on top of its nuclei for "Nucleus Number of Voxels" above 2.05 ie 0.0386238 µm^3 or -3.253887 when logged. 

ggsave(
  paste("mitochondrialVolumes.pdf", sep = ""),
  mitoVolGraph,
  width = 25,
  height = 12,
  units = "cm"
)



### let's do another graph- this time SYBR intensity
vesicleIntensitylistCH1I<- append(vesicleIntensitylist, list(hypCol0vI))
sampleNamesCH1I<- c(graphSampleNames,"hypCol0")


CH1Vesicleintensities<-c()
typeIntlist<-c()

for(i in 1:length(sampleNamesCH1I)){
  vesicleIntensities<- vesicleIntensitylistCH1I[[i]]$Vesicles.mtDNA.Intensity.Mean
  CH1Vesicleintensities<-c(CH1Vesicleintensities,vesicleIntensities)
  
  typeInt<- rep(sampleNamesCH1I[i], length(vesicleIntensities))
  typeIntlist<-c(typeIntlist,typeInt)
  
}


df5<-data.frame(CH1Vesicleintensities,
                typeIntlist
)

intensitiesGraph<-
          ggboxplot(df5,
          x = "typeIntlist",
          y = "CH1Vesicleintensities",
          color = "typeIntlist", outlier.shape = NA) +
  #include outlier shape= NA to remove double-plotted points
  geom_jitter(width = 0.2, aes(color = typeIntlist), alpha = 0.3 ) +
  theme(text = element_text(size = 10), legend.position="none", axis.title.x=element_blank())  + labs(y= "SYBR sphere/vesicle Intensity (A.U)")

#There are some videos (h3,4) that have a sheet of SYBR over the top of the cell, so pick up lots of spots there. 

#not going to quantify number of vesicles as it does NOT reflect mtDNA molecules only. 

ggsave(
  paste("SYBRVesicleIntensities.pdf", sep = ""),
  intensitiesGraph,
  width = 25,
  height = 12,
  units = "cm"
)

