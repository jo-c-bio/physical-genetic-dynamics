#Oxygen electrode data

mypalette = c("#0073C2FF", "#00A36C", "#868686FF","#FFD700")


#read all into a list
#These are direct exports from the hansatech software, as .csv format. 
setwd("~/OxygenElectrodeProper26-7-24/")
temp = list.files(pattern="\\.csv$")
myfiles=  lapply(temp, read.csv, skip = 27, header = F, row.names=NULL)

#rename columns with this
x<-myfiles[[1]][1,]

#apply to all in list
newlist <- lapply(myfiles, "colnames<-", x)

#remove leftover row of names
myfiles2 <- lapply(newlist, tail, -1)
#remove last two rows 
myfiles3 <- lapply(myfiles2, head, -2)

library(plyr)
library(stringr)
min(unlist(lapply(myfiles3, nrow)))
myfiles4<-do.call(rbind.fill, myfiles3)

#fill genoptype variable in, have to fill blanks with NAs before using this
myfiles4 <- myfiles4 %>%
  mutate_at("Label", str_replace, "Stop", "")
myfiles4<-apply(myfiles4, 2, function(x) gsub("^$|^ $", NA, x))
myfiles4<-data.frame(myfiles4)
myfiles5<-myfiles4 %>%
         tidyr::fill(Label)
myfiles6<-myfiles5[,1:4]

myfiles6$Oxygen.1<-as.numeric(myfiles6$Oxygen.1)
myfiles6$Time<-as.numeric(myfiles6$Time)



#plot over time
ggplot(data=myfiles6, aes(x=Time, y=Oxygen.1)) +
  geom_line(aes(colour = Label))

#These values will vary depending on your experiment!
#need to enter fresh weights
#0.06g is 60mg
#We need a value for each dataframe in nmol per min per gram of FW. 
OCR<- function(y,FW){
  nmol = max(myfiles6[myfiles6$"Label"==y,"Oxygen.1"]) - min(myfiles6[myfiles6$"Label"==y,"Oxygen.1"])
  #print(nmol)
  mins = max(myfiles6[myfiles6$"Label"==y,"Time"])/60
  return(nmol/mins/FW)
}

#need to update these with the accurate FW measurements
OCRTable<-cbind(c(
            OCR("AOX1",0.06),
            OCR("AOX2",0.062),
            OCR("AOX3",0.057),
            #OCR("AOX4",0.06),
            OCR("Col1",0.059),
            OCR("Col2",0.059),
            OCR("Col3",0.06),
           # OCR("Col4",0.062),
            OCR("IVD1",0.059),
            OCR("IVD2",0.06),
            OCR("IVD3",0.06)
           # OCR("IVD4",0.059)
            ),
            #remove the 1-4 from sample names here
            #gsub('[0-9]+', '', unique(myfiles6$Label)))
            c(rep("AOX",3),rep("col",3),rep("IVD",3)))
#make sure here that the unique string call matches with the order you've
#put the function calls in. 
colnames(OCRTable)<- c("OCR", "Genotype")
OCRTable<-data.frame(OCRTable)
OCRTable$OCR<-as.numeric(OCRTable$OCR)

ggboxplot(OCRTable, x = "Genotype", y = "OCR", color = "Genotype", palette = mypalette,  outlier.shape = NA) +  
  geom_jitter(width = 0.2,aes(color = Genotype)) +  ylab("OCR (nmol/min/FW)") + scale_y_continuous(limits = c(0,NA)) +
  stat_compare_means()   +   # Add global p-value. Function does this automatically. label.y specification is positioning relative to axes
  theme(text = element_text(size = 10))

ggplot(OCRTable, aes(Genotype, OCR, fill = Genotype)) +
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")+ stat_compare_means() 



