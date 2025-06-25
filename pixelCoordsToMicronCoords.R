#getMaxima data table processing. 

#slightly awkward output for the coordinates of the maxima points out of imagej macro.
#Awkward because they are in pixel coordinates instead of micron coordinates like Trackmate uses.
#the resolution is the same across all ten videos, so we will just divide every number by that and export again. 

setwd("~/cropped/retracked28-2-24/getMaxima")
#read all the .csv files in the getMaxima directory
filelist <- list.files(pattern = ".csv", recursive = TRUE, full.names = TRUE)

#read all of these in recursively
datalist <- lapply(filelist, read.csv)

#set the resolution, need to check all vids are the same before this. Can use getMetadata.ijm to check
resolution<-as.numeric(4.8272)


#Add new columns that are the coordinates in microns, not pixels
datalist.map <-  
  Map(function(x) {
  mutate(x,
         new_col_one = X / resolution,
         new_col_two = Y / resolution)
}, datalist)


#set names to the datalist, these are just their filenames
 names(datalist.map)<- filelist
 
#ImageJ didn't export frame name alongside coordinates, so we'll add them in here, 
#frame is stored within the filename itself, so we'll use that. 
 for(i in 1:length(datalist.map)){
   p<-strsplit(filelist[i],"\\/")[[1]][3]
   q<-strsplit(p,"-")[[1]][1]
   
   col_frame<-rep(q,nrow(datalist.map[[i]]))
   datalist.map[[i]]$new_col <- col_frame
   
 }

#Export these using their file names. will be .csv.txt output. 
 sapply(names(datalist.map), 
        function (x) write.table(datalist.map[[x]], file=paste(x, "txt", sep=".") )   )
 