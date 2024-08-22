#Proportion plot comparing Fixed Imaris output with Live Trackmate output
#Compared presence/absence mtDNA between the two data sets

library(dplyr)
library(ggplot2)

Imaris<- read.csv("/Users/joannachustecki/Documents/PostDoc23-Data/nucleoidQuantification/round2/ProportionalValuesformtDNApresenceAbsenceFromImaris.csv")

vidList<- c("3in1000SYBR-MS-45min002-croppedd.xml",
            "3in1000SYBR-MS-45min003-cropped.xml"	,
            "sample145mnSYBR005-cropped.xml"	,
            "sample345mnSYBR-cropped.xml"	,
            "sample445mnSYBR003-cropped.xml"	,
            "seedling1-45minsybr001-cropped-cell1.xml",
            "seedling1-45minsybr001-cropped-cell2.xml"	,
            "seedling3-45minsybr003-3DdriftCorrect-cell1.xml",
            "seedling3-45minsybr003-3DdriftCorrect-cell2.xml",
            "seedling4-45minsybr002-cropped.xml"	
)

overallVids2a.m<-c()
overallVids2a.m.p<-c()
overallVids2a.m.a<-c()


for(v in 1:length(vidList)){

Strategy2a.m<- read.csv( paste( "/Users/joannachustecki/Documents/PostDoc23-Data/nucleoidQuantification/currentDecentTimelapseSYBR/cropped/retracked28-2-24/trajsatLeastOneAdjacentmtDNA_Maxima_",vidList[v],".csv",sep=""))

with<-count(Strategy2a.m %>% filter(count2a.m....1 == TRUE))
without<-count(Strategy2a.m %>% filter(count2a.m....1 == FALSE))



total<-length(Strategy2a.m$count2a.m....1)

#proportion here is always mtDNA presence/ mtDNA total
mitoProportion <- round(with/total, digits = 3)


row.dataP<- c("mtDNAPresence", paste("vidList",v,sep=""), with, total, mitoProportion)
row.dataA<- c("mtDNAAbsence", paste("vidList",v,sep=""), without, total, mitoProportion)
overallVids2a.m.p<- rbind(overallVids2a.m.p, row.dataP)
overallVids2a.m.a<- rbind(overallVids2a.m.a, row.dataA)
}

overallVids2a.m<- rbind(overallVids2a.m.p,overallVids2a.m.a)

colnames(overallVids2a.m)<-c("Type","realSampleName", "NumberMitos", "mitoTotal","mitoProportion")


meldedTable <- rbind(Imaris[,c("Type","realSampleName","NumberMitos", "mitoTotal", "mitoProportion")],overallVids2a.m )


meldedTable$mitoProportion<-as.numeric(meldedTable$mitoProportion)
meldedTable$NumberMitos<-as.numeric(meldedTable$NumberMitos)
meldedTable$Type<-as.character(meldedTable$Type)
meldedTable$realSampleName<-as.character(meldedTable$realSampleName)
meldedTable$mitoTotal<-as.numeric(meldedTable$mitoTotal)

ggplot(data=meldedTable, aes(x=realSampleName, y=NumberMitos, fill=Type)) +
  geom_bar(stat="identity")+
  #add in a value as a proportion here
  geom_text(aes(y=mitoTotal, label=mitoProportion), vjust=-0.35)+
  
  #add in cell volume here
  #geom_text(aes(y=0, label=round(cellVolumeµm3)),vjust=-0.5,size = 7/.pt)+
  
  #rotate x axis 
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5, hjust=1))+
  #add a title
  ggtitle(paste("upper value = Proportion of mitochondria with mtDNA present",sep=""))     

