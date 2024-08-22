#Contrast threshold passing script, but edited to use maxima as the evidence for mtDNA presence


#need to change the files in the directories to have this name
#vidList<-rep(paste("vid",1:10,sep=""))
vidList <- c(
  "3in1000SYBR-MS-45min002-croppedd.xml",
  "3in1000SYBR-MS-45min003-cropped.xml"	,
  "sample145mnSYBR005-cropped.xml",
  "sample345mnSYBR-cropped.xml"	,
  "sample445mnSYBR003-cropped.xml"	,
  "seedling1-45minsybr001-cropped-cell1.xml",
  "seedling1-45minsybr001-cropped-cell2.xml"	,
  "seedling3-45minsybr003-3DdriftCorrect-cell1.xml",
  "seedling3-45minsybr003-3DdriftCorrect-cell2.xml",
  "seedling4-45minsybr002-cropped.xml"
) 



library(TrackMateR)
library(dplyr)
library(igraph)
library(brainGraph)


directory <-
  "/Users/joannachustecki/Documents/PostDoc23-Data/nucleoidQuantification/currentDecentTimelapseSYBR/cropped/retracked28-2-24/"


matrix.totals <- c()
for (v in 1:length(vidList)) {
  print(v)
  filename <- vidList[v]
  file <- strsplit(as.character(vidList[v]), "\\.")[[1]][1]
  directory <-
    "/Users/joannachustecki/Documents/PostDoc23-Data/nucleoidQuantification/currentDecentTimelapseSYBR/cropped/retracked28-2-24/"
  df <- read.csv(paste(directory,"/getMaxima/",file, ".nd2/",file,"-datasetwithMaximas.txt",sep = ""))
  
  
  #what are the proportions of trajs with 1, 2 or more frame times that cross the contrast threshold?
  #by definition a trajectory is at least 2 frames, so the data is dependant of how many pass that number once you get over 2 frames.
  thosewithnomtdna <- c()
  withatleast1mtDNA <- c()
  withatleast2mtDNA <- c()
  withmorethan2mtDNA <- c()
  withAllPassing <- c()
  
  thosewithnomtdna2 <- c()
  withatleast1mtDNA2 <-c()
  withatleast2mtDNA2 <- c()
  withmorethan2mtDNA2 <- c()
  withAllPassing2 <- c()
  
  d <- 0
  trajsInOrder.1 <- c()
  for (i in 1:max(df$traj)) {
    #Subset by trajectory
    subset <- df[df$traj == i, ]
    if (nrow(subset) > 0) {
      d <- d + 1
      #These go through the subset of a single trajectory ove all frames it is visible in, 
      #and if these conditions are met, say, it has two instances where contrast or maxima passes, it prints true for that traj.
      #They produce TRUE or FALSE for each trajectory, and add it to a growing list.
      
      thosewithnomtdna <-
        c(thosewithnomtdna,
          sum(subset$maximamtDNA == 0) == nrow(subset))

      withatleast1mtDNA <-
        c(withatleast1mtDNA,
          any(subset$maximamtDNA == 1))
      
      withatleast2mtDNA <-
        c(withatleast2mtDNA,
          sum(subset$maximamtDNA == 1) >= 2)
      
      withmorethan2mtDNA <-
        c(withmorethan2mtDNA,
          sum(subset$maximamtDNA == 1) > 2)
      
      withAllPassing <-
        c(withAllPassing,
          sum(subset$maximamtDNA == 1) == nrow(subset))
      
      #Here we collect the trajectory names
      trajsInOrder.1 <- c(trajsInOrder.1, i)
      
    }
  }
  
  max(df$traj)
  #no of trajs that have more than 0 frames
  d
  sum(thosewithnomtdna)
  sum(withatleast1mtDNA)
  sum(withatleast2mtDNA)
  sum(withmorethan2mtDNA)
  sum(withAllPassing)
  
  
  #Let's align these with the trajectories so we can see which 'have' or 'have not' mtDNA by Strategy 1.
  thosewithnomtdna.df <- data.frame(trajsInOrder.1, thosewithnomtdna)
  withatleast1mtDNA.df <-
    data.frame(trajsInOrder.1, withatleast1mtDNA)
  withatleast2mtDNA.df <-
    data.frame(trajsInOrder.1, withatleast2mtDNA)
  withmorethan2mtDNA.df <-
    data.frame(trajsInOrder.1, withmorethan2mtDNA)
  withAllPassing.df <-
    data.frame(trajsInOrder.1, withAllPassing)
  
  
  #export
  write.csv(
    withatleast1mtDNA.df,
    paste(
      directory,
      "trajsWithAtLeastOnemtDNA_Maxima_",
      filename,
      ".csv",
      sep = ""
    )
  )
  
  
  #export
  write.csv(
    withAllPassing.df,
    paste(
      directory,
      "trajsWithAllPointsPassing_Maxima_",
      filename,
      ".csv",
      sep = ""
    )
  )
  
  #export
  write.csv(
    thosewithnomtdna.df,
    paste(
      directory,
      "trajsWithNOPointsPassing_Maxima_",
      filename,
      ".csv",
      sep = ""
    )
  )
  
  #Strategy 2a: How many are there that have at least one ADJACENT contrast threshold passing frames?
  

  #let's count them too
  count2a.m <- c()
  #let's store all the info as a dataframe as it may come in useful
  alladjacentpassingdataforTrajs2a <- c()
  #lets also build a distribution of trajectory length so we can see what we're working with.
  distTrajLength <- c()
  
  trajsInOrder.2a <- c()
  
  for (i in 1:max(df$traj)) {
    #subset data by traj
    subset <- df[df$traj == i, ]
    currentlist <- c()
    #let's store all the info as a dataframe as it may come in useful
    nof2adjacentpassing <- c()
    p <- c()
    #is the trajectory actually populated?
    if (nrow(subset) > 0) {
      distTrajLength <- c(distTrajLength, nrow(subset))
      #got through the subset, are any frames that are adjacent both crossing the contrast threshold?
      for (j in 1:nrow(subset)) {
        #this sum equation puts NA at the end because the last value in the dataframe doesn't have a value to go to next
        #list of 1 or 0 depending if the adjacent frames passed
        currentlist <- c(currentlist,
                         (
                           sum(
                             subset$maximamtDNA[j] == 1 &
                               subset$maximamtDNA[j + 1] == 1,
                             na.rm = TRUE
                           )
                         ))
      }
      #look at the list, count the number of instances where frames passed the requirement to have an adjacent frame pass also
      count2a.m <- c(count2a.m , (sum(currentlist)))
      
      #build a list of trajs we're working on
      trajsInOrder.2a <- c(trajsInOrder.2a, i)
      
      # if you need a more detailed look into the data frame
      nof2adjacentpassing <- c(nof2adjacentpassing, currentlist)
      p <- (data.frame(i, nof2adjacentpassing))
      alladjacentpassingdataforTrajs2a <-
        rbind(alladjacentpassingdataforTrajs2a, p)
    }
  }
  
  hist(distTrajLength, breaks = seq(0, max(df$t), 1))
  
  
  sum(count2a.m  > 1)
  
  #Let's align these with the trajectories so we can see which 'have' or 'have not' mtDNA by Strategy 2a.
  #we have here the trajectories, and a true or false if they do have >1
  atLeastOneAdjacent <- data.frame(trajsInOrder.2a, count2a.m  >= 1)
  #note; count2a on its own doesn't make much sense, it would be how many instances of two adjacent frames you have per trajectory, when we're really interested
  #in whether these trajectories have >= 1 instances of two frames adjacently having passing frames
  
  #export
  write.csv(
    atLeastOneAdjacent,
    paste(
      directory,
      "trajsatLeastOneAdjacentmtDNA_Maxima_",
      filename,
      ".csv",
      sep = ""
    )
  )
  
  #Strategy 2b: How many are there that have three ADJACENT contrast threshold passing frames?
  

  #let's count them too
  count2b.m <- c()
  #let's store all the info as a dataframe as it may come in useful
  alladjacentpassingdataforTrajs2b <- c()
  
  trajsInOrder.2b <- c()
  
  for (i in 1:max(df$traj)) {
    #subset data by traj
    subset <- df[df$traj == i, ]
    currentlist <- c()
    #let's store all the info as a dataframe as it may come in useful
    nof3adjacentpassing <- c()

    #is the trajectory actually populated?
    if (nrow(subset) > 0) {
      #go through the subset, are any frames that are adjacent both crossing the contrast threshold?
      for (j in 1:nrow(subset)) {
        #this sum equation puts NA at the end because the last value in the dataframe doesn't have a value to go to next
        #list of 1 or 0 depending if the adjacent frames passed
        currentlist <- c(currentlist,
                         (
                           sum(
                             subset$maximamtDNA[j] == 1 &
                               subset$maximamtDNA[j + 1] == 1 &
                               subset$maximamtDNA[j + 2] == 1,
                             na.rm = TRUE
                           )
                         ))
      }
      #look at the list, count any that had => 1 occurences of adjacent passing
      count2b.m <- c(count2b.m, (sum(currentlist)))
      
      #build a list of trajs we're working on
      trajsInOrder.2b <- c(trajsInOrder.2b, i)
      
      # if you need a more detailed look into the data frame
      nof3adjacentpassing <- c(nof3adjacentpassing, currentlist)
      p <- (data.frame(i, nof3adjacentpassing))
      alladjacentpassingdataforTrajs2b <-
        rbind(alladjacentpassingdataforTrajs2b, p)
    }
  }
  
  
  sum(count2b.m > 1)
  
  #Let's align these with the trajectories so we can see which 'have' or 'have not' mtDNA by Strategy 2b.
  #we have here the trajectories, then how many adjacent frames the traj has that pass contrast, and a true or false if they do have >1
  atLeastThreeAdjacent <- data.frame(trajsInOrder.2b, count2b.m >= 1)
  
  #export
  write.csv(
    atLeastThreeAdjacent,
    paste(
      directory,
      "trajsatLeastThreeAdjacentmtDNA_Maxima_",
      filename,
      ".csv",
      sep = ""
    )
  )

  
  #Okay, we have three strategies. Now, how do we build networks from this?
  
  
  
  #let's make a data table that accumulates
  
  currentlist <- matrix(c(d,
                          sum(withatleast1mtDNA) ,
                          sum(count2a.m > 1) ,
                          sum(count2b.m > 1),
                          sum(withAllPassing)))

  
  matrix.totals <- cbind(matrix.totals, currentlist)
  
  
}



colnames(matrix.totals) <- c(vidList)
rownames(matrix.totals) <-
  c("TotalValidTrajs",
    "withatleast1mtDNA",
    "withatLeastOneAdjacent",
    "withatLeastThreeAdjacent",
    "withAllMustBePassing"
  )

#Make a barplot

Type<-rep(rownames(matrix.totals),10)

library("tidyr")
library("ggplot2")
gth<- data.frame(matrix.totals)

my_data2<-gather(gth)
my_data2<-data.frame(my_data2,Type)
my_data2

#lets include a list of pseudonyms for the videos to save labelling space

namesF <- function() {
  names<-c()
  for (i in 1:10) {
    names <- c(names, rep(paste("cell", i, sep = ""), 4))
  }
  return(names)
}
names<-namesF()

my_data2<-cbind(my_data2,names)

#order as 1 - 10
my_data2$names = factor(my_data2$names, 
                        levels=c("cell1","cell2","cell3","cell4","cell5","cell6","cell7","cell8","cell9","cell10"))

TotalsGraph<-ggplot(my_data2, aes(fill=Type, y=value, x=names)) + 
              geom_bar(position="dodge", stat="identity") + 
              theme(axis.text.x = element_text(angle = 20, vjust = 0.5, size= 15) , 
                    axis.title.x=element_blank(),
                    axis.text.y = element_text(size = 15),
                    axis.title.y=element_text(size = 15),
                    legend.text=element_text(size=13)) + 
              ylab("Trajectory Number")
TotalsGraph
ggsave(
  paste("presenceAbsenceDefinitionsTotalsMaximaMethodv2.pdf", sep = ""),
  TotalsGraph,
  width = 20,
  height = 7,
  units = "cm"
)

