##install.packages("plotrix")
##install.packages("sp")
library(sp)
library(plotrix)
library(TrackMateR)
library(dplyr)
library(igraph)
library(brainGraph)
library(ggpubr)
library(ggplot2)


#Functions

readDataframe <- function(filename) {
  tmObj <-
    readTrackMateXML(XMLpath = paste(directory, filename, sep = ""))
  x <- tmObj[1]
  dfBig <- data.frame(x)
  
  #debugging incase Detected multiple spots per frame for one or more tracks.
  #This can happen during manual track editing.
  a <- cbind(dfBig$trace, dfBig$frame)
  which(duplicated(a) == TRUE)
  if (length(which(duplicated(a) == TRUE)) >= 1) {
    print("check this against track viewer")
  }
  
  #Let's make a new dataframe with the relevant data as the existing one is huge
  #Track position, spot intensity, time, contrast, we'll keep name because that will help us trace back to individual mitos
  dfSmall <-
    data.frame(
      dfBig$name,
      dfBig$trace,
      dfBig$t,
      dfBig$frame,
      dfBig$x,
      dfBig$y,
      dfBig$mean_intensity,
      dfBig$contrast ,
      dfBig$radius
    )
  #This is contrast for channel 1, which is the SYBR channel
  #To get it into the format we used for our code
  colnames(dfSmall) <-
    c("ID",
      "traj",
      "secs",
      "t",
      "x",
      "y",
      "mean_intensity",
      "contrast",
      "radius")
  dfSmall$traj <- as.numeric(dfSmall$traj)
  dfSmall$traj <- dfSmall$traj + 1
  dfSmall$t <- dfSmall$t + 1
  df <- dfSmall
  
  return(df)
}


getCorrespondingMitos <- function(currentDFj, subsetTMj) {

  subsetTMj$maximamtDNA <- rep(NA, nrow(subsetTMj))
  subsetTMj$numberPassing <- rep(NA, nrow(subsetTMj))
  #we're drawing the circle encompassing the mito as it's used in Trackmate.
  #x,y&radius come from Trackmate output.
  for (m in 1:nrow(subsetTMj)) {
    currentMito <-
      draw.circle(subsetTMj$x[m], subsetTMj$y[m], subsetTMj$radius[m])
    #Then, go through all the maximal points we found for that frame, test; are any
    #within the mito we just isolated from the trackmate data?
    
    testOutput <- c()
    
    for (maxPoint in 1:nrow(currentDFj)) {
      pointTestOut <-
        point.in.polygon(
          currentDFj$new_col_one[maxPoint],
          #X IN MICRONS
          currentDFj$new_col_two[maxPoint],
          #Y IN MICRONS
          currentMito$x,
          currentMito$y,
          mode.checked = FALSE
        )
      
 #     print(pointTestOut)
      
      testOutput <- c(testOutput, pointTestOut)
      
    }
    if (sum(testOutput) >= 1) {
      subsetTMj$maximamtDNA[m] <- 1
      subsetTMj$numberPassing[m] <- sum(testOutput >= 1)
    } else {
      subsetTMj$maximamtDNA[m] <- 0
    }
  }
  return(subsetTMj)
}

generateTableWithMaximaData <- function(datalist,df) {
  for (j in 1:length(datalist)) {
    print(j)
    #go through the datalist, aka list of findMaxima outputs, by frame.
    #convoluted as the frames didn;t sit in 1-10 order as they were characters in the filenames
    currentFramej <- datalist[[j]][6][[1]][1]
    #grab the whole dataframe for the maxima points at the current frame
    currentDFj <- datalist[[j]]
    #which trajectories correspond to this frame in the Trackmate output.
    subsetTMj = df[df$t == currentFramej, ]
    
    
    subsetTMj.m <- getCorrespondingMitos(currentDFj, subsetTMj)
    
    #datasetwithMaximas <- rbind(datasetwithMaximas, subsetTMj.m)
    datalist.withMaximas <- append(list(subsetTMj.m), datalist.withMaximas)
  }
  return(datalist.withMaximas)
}

#if numberPassing > 1 go back and look at the data, how many maxima does that mito have in it?
#it seems that most of the time, 2 maxima in a mito is fine, there is another mito nearby and it has also been detected. 








#ARGS

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

for(i in 1:length(vidList)) {
  #Drafting the use of getMaxima coordinates as the threshold to mtDNA
  #this grabs the filename, without .xml at the end
  file <- strsplit(vidList[i], "\\.")[[1]][1]
  filename <- vidList[i]
  directory <-
    "/Users/joannachustecki/Documents/PostDoc23-Data/nucleoidQuantification/currentDecentTimelapseSYBR/cropped/retracked28-2-24/"

  df <- readDataframe(filename)
  
  #read all the .csv files in the getMaxima directory
  filelistMaxima <-
    list.files(
      path = paste(directory, "/getMaxima/", file, ".nd2/", sep = ""),
      # $ at the end of pattern means end of string
      pattern = ".csv.txt$",
      recursive = TRUE,
      full.names = TRUE
    )
  
  #read all of these in recursively
  datalist <- lapply(filelistMaxima, read.csv, sep = " ")
  
  #re-set params
  datalist.withMaximas <- list()
  datasetwithMaximas <- c()
  
  #must have a plot for the draw.circles to work on, even though we're just saving it as a variable
  plot(100, 100)
  
  #run fucntions
  datalist.withMaximas <- generateTableWithMaximaData(datalist, df)
  datasetwithMaximas <-
    do.call(rbind.data.frame, datalist.withMaximas)
  
  #let's get these files exported for completion.
  #Export these using their file names. will be .csv.txt output.
  names(datalist.withMaximas) <- filelistMaxima
  sapply(names(datalist.withMaximas),
         function (x)
           write.csv(datalist.withMaximas[[x]], file = paste(x, "-withMaxima.csv", sep =  "")))
  
  write.csv(datasetwithMaximas,paste(directory, "/getMaxima/", file, ".nd2/",file, "-datasetwithMaximas.txt", sep = ""))
  #end loop
  
}



