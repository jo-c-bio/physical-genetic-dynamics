#This script reads in all data as trackmate files, and builds per trajectory physical statistics, such as speed, intermitochondrial distance etc.
#It also reads corresponding contrast or maxima definition files, and compares the physical statistics of those defined as with or without
#mtDNA (over six definitons).
#Plots are made at the end of the script.



#These can all be bought in as args, like allfiles in a directory, but all need corresponding contrast values.
#if bringing in filename as entire files directory, will need to alter that below
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

contrastList<-
  c(as.numeric(0.05),
    as.numeric(0.18),
    as.numeric(0.16),
    as.numeric(0.21),
    as.numeric(0.22),
    as.numeric(0.17),
    as.numeric(0.16),
    as.numeric(0.07),
    as.numeric(0.04),
    as.numeric(0.11))




# this script is organised where each physical statistic is generate in its own function,
#and then physicalStatisticsGeneration() is run, where the definitions are implemented.
#This make take a while for large datasets
#Graph generation at end of script. 

library(TrackMateR)
library(dplyr)
library(graphics)
library(sp)

directory<-"~/cropped/retracked28-2-24/"

#function for reading and formatting trackmate data frames
readDataframe <- function(xmlFilename) {
  tmObj <-readTrackMateXML(XMLpath = paste(directory, xmlFilename, sep = ""))
  x <- tmObj[1]
  dfBig <- data.frame(x)
  
  #Let's make a new dataframe with the relevant data as the existing one is huge
  #Track position, spot intensity, time, contrast, we'll keep name because that will help us trace back to individual mitos
  dfSmall <- data.frame(dfBig$name,dfBig$trace,dfBig$t, dfBig$frame,dfBig$x, dfBig$y,dfBig$mean_intensity, dfBig$contrast ) 
  #This is contrast for channel 1, which is the SYBR channel
  #To get it into the format we used for our code
  colnames(dfSmall) <- c("ID","traj","secs", "t","x","y","mean_intensity", "contrast")
  dfSmall$traj <- as.numeric(dfSmall$traj)
  dfSmall$traj <- dfSmall$traj + 1
  dfSmall$t <- dfSmall$t + 1
  df <- dfSmall
  print(paste("trajectorymaxis",max(dfSmall$traj) ))
  return(df)
}

# loop through trajectories
speedsFunction <- function(df, ntraj, frameTime) {
  message("Computing speeds...")
  speeds = NULL
  meanSpeeds <- c()
  trajs<-c()
  sumTrajDistance <- c()
  trajDuration <- c()
  for (t in ntraj) {
    # grab trajectories from this frame
    #So here, we'll need to divide the distance by the frame time.
    subset = df[df$traj == t, ]
    
    #had to put this catch in as I have some videos where they skip over trajectory numbers, eg video 1 goes 28,29,31...
    #cannot fix it at source in the trackmate file.
    if (nrow(subset) > 0) {
      #in order to collect mean speed per trajectory
      speedsperTraj <- c()
      trajDistances <- c()
      
      for (i in 2:nrow(subset)) {
        if (subset$t[i] == subset$t[i - 1] + 1) {
          pair.dist.2 = (subset$x[i] - subset$x[i - 1]) ** 2 + (subset$y[i] - subset$y[i - 1]) ** 2
          this.speed = sqrt(pair.dist.2) / frameTime
          speeds = c(speeds, this.speed)
          speedsperTraj <- c(speedsperTraj, this.speed)
          
          #we can also collect trajectory distance travelled.
          trajDistances<- c(trajDistances,pair.dist.2)
        }
      }
      meanSpeeds <- c(meanSpeeds, mean(speedsperTraj))
      trajs<-c(trajs,t)
      sumTrajDistance<- c(sumTrajDistance, sum(trajDistances))
      #let's also collect frame duration of that trajectory
      trajDuration <- c(trajDuration, nrow(subset))
    }
  }
  sp<-data.frame(trajs,meanSpeeds,sumTrajDistance,trajDuration)
  #return(speeds)
  return(sp)
  
}


minDistancesfunction <- function(maxt,ntraj,threshold,df) {
  amlist = data.frame(firstframe = NULL, t1 = NULL, t2 = NULL)
  am = matrix(nrow = ntraj, ncol = ntraj)
  colocal.time = matrix(0, nrow = ntraj, ncol = ntraj)
  dist.frame = data.frame(frame = NULL, mean.min.dist = NULL)
  vlist = NULL
  coloc.traj.means <-c()
  
  df.min.dists <- data.frame()
  
  # loop through frames
  message("Computing distances...")
  for (t in 1:maxt) {
    #Jc adding: empty all dists list
    allDists <- c()
    # store total number of trajectories recorded up to this frame (to account for singletons)
    vlist = c(vlist, length(unique(df$traj[df$t <= t])))
    # grab trajectories from this frame
    subset = df[df$t == t, ]
    min.dists = NULL
    # loop through pairs of trajectories
    for (i in 1:nrow(subset)) {
      min.dist.2 = -1
      for (j in 1:nrow(subset)) {
        if (i != j) {
          # compute distance for this pair -- recording minimum dists
          pair.dist.2 = (subset$x[i] - subset$x[j]) ** 2 + (subset$y[i] -
                                                              subset$y[j]) ** 2
          #JC: records minimum distance between the current trajectory and all other trajectories
          if (min.dist.2 == -1 |
              pair.dist.2 < min.dist.2) {
            min.dist.2 = pair.dist.2
          }
          #so here, we can flag for that single trajectory, and it's minimum distance to all other trajectories at that time point.
          subset[i, "min.dist.everymito"] <- sqrt(min.dist.2)
          #JC added: if we want hists of all dists, let's take all recorded distances between individuals, sqrt here because it's key for Euc. distance but Iian doe it at the end of this code instead.
          #allDists<-c(allDists,sqrt(pair.dist.2))
          # if this distance is below our threshold
          if (pair.dist.2 < threshold ** 2) {
            # if we haven't recorded this pair yet, do so
            if (is.na(am[subset$traj[i], subset$traj[j]])) {
              amlist = rbind(amlist,
                             data.frame(
                               firstframe = t,
                               t1 = subset$traj[i],
                               t2 = subset$traj[j]
                             ))
              am[subset$traj[i], subset$traj[j]] = am[subset$traj[j], subset$traj[i]] = t
            }
            # add colocalisation time to matrix
            if (i < j) {
              colocal.time[subset$traj[i], subset$traj[j]] = colocal.time[subset$traj[j], subset$traj[i]] = colocal.time[subset$traj[j], subset$traj[i]] +  1
            }
          }
        }
      }
      # summarise distance statistics
      min.dists = c(min.dists, sqrt(min.dist.2))
    }
    dist.frame = rbind(dist.frame, data.frame(frame = t, mean.min.dist = mean(min.dists)))
    
    #JC again, adding:
    #distsplot<-hist(allDists)
    #distsplot
    #we want to bind the subsets to each other so we can get a full reading of trajectory, by time with the minimum distances.
    #table with min dists
    df.min.dists <- rbind(df.min.dists, subset)
    
  }
  
  df.min.dists.mean <- trajAndDistMean.maker(df.min.dists)
  df.min.dists.mean. <- arrange(df.min.dists.mean,traj)
  
  
  #colocalisation time is currently in a traj-traj matrix,so we'll simplify it and
  #bind it to the min dists df to be able to use the same return function.
  #we remove the 0 as mean(0,1,2,3) is 0. We removed when averaging in the original code too. 
  colocal.time[colocal.time == 0] <- NA
  coloc.traj.means <-rowMeans(colocal.time,na.rm=TRUE)
  coloc.traj.means <- data.frame(c(1:nrow(colocal.time)),coloc.traj.means )
  colnames(coloc.traj.means) <- c("traj","coloc.traj.mean")
  
  #need to keep only those trajs that appear in the dists lists. (trajs non continuous, this is a weird catch i can't fix at source)
  dd<-coloc.traj.means$traj[(coloc.traj.means$traj %in% df.min.dists.mean.$traj)]
  coloc.traj.means. <- coloc.traj.means[dd,]
  
  print("passed so far")
  export<- cbind(df.min.dists.mean.,coloc.traj.means.)
  return(list(export, df.min.dists))
}

#function to get mean per traj of intermitochondrial distances
trajAndDistMean.maker<- function(df){
  trajAndDistMean.df<-data.frame()
  for(nt in unique(df$traj)) {
    trajAndDistMean <- c ( nt, 
                           mean(df[df$traj==nt, "min.dist.everymito"]))
    trajAndDistMean.df <- rbind(trajAndDistMean.df,trajAndDistMean)
  }
  colnames(trajAndDistMean.df)<-c("traj","minimumDistMeanPerTraj")
  return(trajAndDistMean.df)
}



createCHMapPlot<-function(df, MaxTraj = ntraj, directory, xmlFilename){
  message("Computing areas...")
  arealist<-c()
  logarealist<-c()
  trajlist<-c()
  wholecell<-data.frame(df$x,df$y)
  pdf(paste(directory,"/",xmlFilename,"ConvexHullMap.pdf",sep=""))
  plot(wholecell,xlab = "x", ylab="y",main = xmlFilename,cex=0.5)
  
  #not sure yet If i need a catch like this
  filenumbersthatexist<-c()
  for (i in 1:MaxTraj){
    if(i %in% df$traj){
      filenumbersthatexist<-c(filenumbersthatexist,i)
    }
  }
  #Need this catch as sometimes the simulation just doesnt produce files for certain numbers
  for(j in filenumbersthatexist){
    trajfile<- df %>% filter(traj == j)
    #had to put this if catch in to exclude any trajectories that were only one coordinate long, as they would generate an area of 0
    
    if(nrow(trajfile)>=3){
      trajfile<-trajfile[,c("x","y")]
      hpts <- chull(trajfile)
      hpts <- c(hpts, hpts[1])
      #Need to include this catch as some simulation tejectories only have two furthest points- meaning it makes a line, not a polygon. Polygon needs at least 4 coordinates
      #Or you get the error Polygon(trajfile[hpts, ], hole = F) : less than 4 coordinates in polygon
      if(length(hpts) > 3){
        
        rn1<- sample(0:1,1)
        rn2<- sample(0:1,1)
        rn3<- sample(0:1,1)
        
        #plot the convex hull as a polygon shape (this depends on the "sp" package being installed)
        chull.poly <- Polygon(trajfile[hpts, ], hole=F)
        #find its area
        chull.area <- chull.poly@area
        #Fill the middle of the polygon transparently
        polygon(trajfile[hpts,],col=rgb(rn1,rn2,rn3,1/4))
        #keep all the areas for each trajectory in a list
        arealist<-c(arealist,chull.area/5)
        logarealist<-c(logarealist,log(chull.area/5))
        #keep all the trajectory names in a list
        trajlist<-c(trajlist,j)
        #add a line to the plot for every new trajectory you want to add
        lines(trajfile[hpts, ], col=rgb(rn1,rn2,rn3,3/4))
      } 
    } else {
      
      arealist<-c(arealist,NA)
      logarealist<-c(logarealist,NA)
      trajlist<-c(trajlist,j)
    }
  }
  dev.off()
  
  convexHulldf<- data.frame(trajlist,arealist,logarealist)
  return(convexHulldf)
  #write.csv(data.frame(trajlist,arealist,logarealist), paste(directory,"/areaandlogarealist",xmlFilename,".csv",sep=""))
}

#this is a function to plot convex hull maps for those with and then without mtDNA seperately to visualise the relationship. 
convexHullPlotWithWithout<-function(df,dd,ntraj,directory,xmlFilename){
  #we need the wholedata frame, not just summary stats so we use df not tempdf
  df$label <-"withOutmtDNA"
  df$label[which(df$traj %in% dd)]<-"withmtDNA"
  
  #seperate them into two datafiles, with and without
  WOxmlFilename<-paste(xmlFilename,"withOutmtDNA", sep="")
  WxmlFilename<-paste(xmlFilename,"withmtDNA", sep="")
  
  dfWO = df %>% filter(label == "withOutmtDNA")
  dfW = df %>% filter(label == "withmtDNA")
  #sent them to be plotted by the same function as we use for all the trajectories. 
  convexHulldf<-createCHMapPlot(dfWO, ntraj, directory, WOxmlFilename)
  convexHulldf<-createCHMapPlot(dfW, ntraj, directory, WxmlFilename)
}


#this is the main function using all others, to generate one data table that gives statistics per trajectory.

physicalStatisticsGeneration<-function(StrategyLongName,StrategyShortName){
  overallDF<-c()
  overallOki<-c()
  justoki<-c()
  
  strategy<-StrategyShortName
  
  for(i in 1:length(vidList)){
    print(i)
    xmlFilename<- vidList[i]
    directory<-"~/cropped/retracked28-2-24/"
    strategyFile<- read.csv( paste(directory, StrategyLongName,xmlFilename,".csv",sep="")) 
    contrastThreshold <- contrastList[i]
    df <- readDataframe(xmlFilename)
    # overall stats for this dataset
    maxt = max(df$t)
    ntraj = max(df$traj)
    #we'll take the difference in secs between the first two value as the frametime, it should be consistent over the video. Also round.
    frameTime= df$secs[2]-df$secs[1]
    #Here we impose a colocalisation threshold for what makes an interaction, in microns
    threshold = 1.6
    
    #mean speed per trajectory
    meanSpeedsperTraj <- speedsFunction(df,c(1:ntraj),frameTime)
    
    #Previously with intermitochondrial distance, we took it per frame, which makes sense to do.
    #here, we are doing the average over each trajectory, which makes less sense, as the mitos are moving through the rest of the population
    #however, the presence/absence is done by traj, so we need a per traj value. 
    meanMinDistsAndColocperTraj.toSplit <- minDistancesfunction(maxt,ntraj,threshold,df)
    meanMinDistsAndColocperTraj<-meanMinDistsAndColocperTraj.toSplit[[1]]
    oki<-meanMinDistsAndColocperTraj.toSplit[[2]]
    meanMinDistsAndColocperTraj$logminimumDistMeanPerTraj <- log(meanMinDistsAndColocperTraj$minimumDistMeanPerTraj)
    meanMinDistsAndColocperTraj$logcoloc.traj.mean <- log(meanMinDistsAndColocperTraj$coloc.traj.mean)
    
    #here we are generating a list per cell of convex hull areas, in µm^2, for each traj
    convexHulldf<-createCHMapPlot(df, ntraj, directory, xmlFilename)
    
    #Then we can combine the stats to one data frame, as should all be relating to per traj stats, keep traj names for now to double check. 
    #stats combine
    statsCombine <- cbind(meanSpeedsperTraj,log(meanSpeedsperTraj$meanSpeeds),
                          log(meanSpeedsperTraj$sumTrajDistance),
                          log(meanSpeedsperTraj$trajDuration),
                          meanMinDistsAndColocperTraj,
                          convexHulldf)
    
    #cell lettered
    cell<-print(rep(letters[i],nrow(statsCombine)))
    
    tempdf<-data.frame(cell,statsCombine)
    print(paste("tempDF Is ",nrow(tempdf)))
    
    #Split the rows of the data frame by if contrast is over the threshold
    #passes if at least one of the frames passes the threshold.
    subsetY = df[df$contrast >= contrastThreshold,]
    subsetN = df[df$contrast <= contrastThreshold,]
    #Find out which unique trajectories pass the threshold
    yesTrajs<-unique(subsetY$traj)
  
    
    #Find out which unique trajectories pass the threshold, according to the strategy we are using. 
    trueSubset =  strategyFile[strategyFile[,3]==TRUE,2]
    yesTrajs = trueSubset[!is.na(trueSubset)]
    #dont need but might
    falseSubset =  strategyFile[strategyFile[,3]==FALSE,2]
    subsetN.a = falseSubset[!is.na(falseSubset)]
    
    print(paste("lengths are ",length(subsetY) + length(subsetN)))
    
    
    #create a new column with the mtDNA presence/absence status for each trajectory
    tempdf$label <-"withOutmtDNA"
    dd<-meanSpeedsperTraj$trajs[(meanSpeedsperTraj$trajs %in% yesTrajs)]
    tempdf$label[ which(meanSpeedsperTraj$trajs %in% dd)]<-"withmtDNA"
    
    #this little section is for calculating minimum mitochondrial distance between each mito in each traj at one frame. 
    #first need to label them by positive or negative
    #  
    oki$label <-"withOutmtDNA"
    ddd<-oki$traj[(oki$traj %in% yesTrajs)]
    oki$label[ which(oki$traj %in% ddd)]<-"withmtDNA"
    #
    okiYes<-
      oki %>%
      filter(label == "withmtDNA") %>% #filter by if +ve
      group_by(t) %>% #group by frames
      summarise(min.dist.everymito_perframe = mean(min.dist.everymito)) #mean the min dists
    okiYes$label <-"withmtDNA"
    #
    okiNo<-
      oki %>%
      filter(label == "withOutmtDNA") %>% #filter by if -ve
      group_by(t) %>% #group by frames
      summarise(min.dist.everymito_perframe = mean(min.dist.everymito)) #mean the min dists
    okiNo$label <-"withOutmtDNA"
    # For a data table of every intermitochondrial distance averaged PER FRAME
    okiAll<-rbind(okiYes, okiNo)
    cellOki<-print(rep(letters[i],nrow(okiAll)))
    okiAll<-cbind(okiAll, cellOki)
    #
    overallOki<- rbind(overallOki,okiAll)
    
    #  For a data table of every single, non-averaged intermitochondrial distance averaged
    justokicell<-print(rep(letters[i],nrow(oki)))
    oki.c<-cbind(oki,justokicell)
    justoki<-rbind(justoki,oki.c)
    #
    
    overallDF<-rbind(overallDF,tempdf)
    
    #sub-function just to plot the convex hulls of with/withoutmtDNA 
    convexHullPlotWithWithout(df,dd,ntraj,directory,xmlFilename)
    
  }
  
  return(list(overallDF,overallOki,justoki,strategy))
  #return(list(overallDF,strategy))
}

plotFunction<-function(dataframe, specificdata,specificdata_char, model_char,yaxis,pvaltable){
  plott <- ggplot(data = dataframe ,aes(x=label, y= specificdata, colour = cell)) +
    geom_jitter(alpha=0.35)+
    stat_summary(
      geom = "point",
      fun = "mean",
      col = "white",
      size = 3,
      shape = 23,
      fill = "red",
      alpha=0.8
    )+ 
    ggtitle(paste(
      "Bonf. adjusted Anova output= Pr(>Chisq) = ",signif(pvaltable[model_char,"Bonferroni"],3), "\n",
      "withOutmtDNA n = ",length(which(!is.na(dataframe[dataframe$label == "withOutmtDNA",specificdata_char]))), ", ",
      "withmtDNA n  = ",length(which(!is.na(dataframe[dataframe$label == "withmtDNA",specificdata_char]))),sep="")) + 
    ylab(yaxis)+
    theme(
      legend.text = element_text(size = 12), 
      legend.title = element_text(size = 14),
      axis.text=element_text(size=12),
      axis.title.x=element_blank(),
      axis.title=element_text(size=14)
    )
  
  return(plott)
}


####### Here we are running these defined functions. ########


#ok.n<-physicalStatisticsGeneration(strategyLongName,StrategyShortName)
  StrategyShortNames<-c("Strategy1",
                       "Strategy2a",
                       "Strategy2b",
                       "Strategy1.m",
                       "Strategy2a.m",
                       "Strategy2b.m")
  strategyLongNames<-c("trajsWithAtLeastOnemtDNA_Contrast_",
                       "trajsatLeastOneAdjacentmtDNA_Contrast_",
                       "trajsatLeastThreeAdjacentmtDNA_Contrast_",
                       "trajsWithAtLeastOnemtDNA_Maxima_",
                       "trajsatLeastOneAdjacentmtDNA_Maxima_",
                       "trajsatLeastThreeAdjacentmtDNA_Maxima_")
ok.all<-list()
  for(s in 1:length(StrategyShortNames)){
    ok.n<-physicalStatisticsGeneration(strategyLongNames[s],StrategyShortNames[s])
    ok.all<-append(ok.all,ok.n)
  }


for(i in seq(from=1, to=24, by=4)){
  ok.n<- ok.all[i:(i+3)]
  print(c(i,i+1,i+2,i+3))



#log new table
ok.n[[2]]$log.min.dist.everymito_perframe<-log(ok.n[[2]]$min.dist.everymito_perframe)
ok.n[[2]]<-ok.n[[2]] %>%  filter(log.min.dist.everymito_perframe!="-Inf")
#change cell column name
names(ok.n[[2]])[names(ok.n[[2]]) == 'cellOki'] <- 'cell'
#define strategy 
strategy <- ok.n[[4]]
strategyName<-  ok.n[[4]]


library(ggplot2)
library(ggpubr)
library(car)
library(lme4)


### revamped stats- using blocking factor anova, using cell as the blocks, as it describes the known variance. 

ok<-ok.n[[1]]

#two ways to do it , 
#but I think first, using aov() won't work correctly as  you need to wrap them in Error(), 
#and your data has to be balanced, i.e. have equal replication across all treatment.
#so, using the lmer() function to specify the correct block as random 


summary(aov(meanSpeeds~as.factor(label)+cell, data=ok))
Anova(lmer(meanSpeeds ~ label + (1 | cell), ok, REML=FALSE))

summary(aov(log.meanSpeedsperTraj.meanSpeeds.~as.factor(label)+cell, data=ok))
Anova(lmer(log.meanSpeedsperTraj.meanSpeeds. ~ label + (1 | cell), ok, REML=FALSE))

summary(aov(minimumDistMeanPerTraj~as.factor(label)+cell, data=ok))
Anova(lmer(minimumDistMeanPerTraj ~ label + (1 | cell), ok, REML=FALSE))

summary(aov(logminimumDistMeanPerTraj~as.factor(label)+cell, data=ok))
Anova(lmer(logminimumDistMeanPerTraj ~ label + (1 | cell), ok, REML=FALSE))

summary(aov(min.dist.everymito_perframe~as.factor(label)+cell, data=ok.n[[2]]))
Anova(lmer(min.dist.everymito_perframe ~ label + (1 | cell), ok.n[[2]], REML=FALSE))

summary(aov(log.min.dist.everymito_perframe~as.factor(label)+cell, data=ok.n[[2]]))
Anova(lmer(log.min.dist.everymito_perframe ~ label + (1 | cell),  ok.n[[2]], REML=FALSE))

summary(aov(coloc.traj.mean~as.factor(label)+cell, data=ok))
Anova(lmer(coloc.traj.mean ~ label + (1 | cell), ok, REML=FALSE))

summary(aov(logcoloc.traj.mean~as.factor(label)+cell, data=ok))
Anova(lmer(logcoloc.traj.mean ~ label + (1 | cell), ok, REML=FALSE))

summary(aov(arealist~as.factor(label)+cell, data=ok))
Anova(lmer(arealist ~ label + (1 | cell), ok, REML=FALSE))

summary(aov(logarealist~as.factor(label)+cell, data=ok))
Anova(lmer(logarealist ~ label + (1 | cell), ok, REML=FALSE))

summary(aov(sumTrajDistance~as.factor(label)+cell, data=ok))
Anova(lmer(sumTrajDistance ~ label + (1 | cell), ok, REML=FALSE))

summary(aov(log.meanSpeedsperTraj.sumTrajDistance.~as.factor(label)+cell, data=ok))
Anova(lmer(log.meanSpeedsperTraj.sumTrajDistance. ~ label + (1 | cell), ok, REML=FALSE))

summary(aov(trajDuration~as.factor(label)+cell, data=ok))
Anova(lmer(trajDuration ~ label + (1 | cell), ok, REML=FALSE))

summary(aov(log.meanSpeedsperTraj.trajDuration.~as.factor(label)+cell, data=ok))
Anova(lmer(log.meanSpeedsperTraj.trajDuration. ~ label + (1 | cell), ok, REML=FALSE))

#Bonferroni correction 

# Compute adjusted p-values

p.unadj.rv<-c(  Anova(lmer(meanSpeeds ~ label + (1 | cell), ok, REML=FALSE))[[3]],
                
                Anova(lmer(log.meanSpeedsperTraj.meanSpeeds. ~ label + (1 | cell), ok, REML=FALSE))[[3]],
                
                Anova(lmer(minimumDistMeanPerTraj ~ label + (1 | cell), ok, REML=FALSE))[[3]],
                
                Anova(lmer(logminimumDistMeanPerTraj ~ label + (1 | cell), ok, REML=FALSE))[[3]],
                
                Anova(lmer(min.dist.everymito_perframe ~ label + (1 | cell), ok.n[[2]], REML=FALSE))[[3]],
                 
                Anova(lmer(log.min.dist.everymito_perframe ~ label + (1 | cell),  ok.n[[2]], REML=FALSE))[[3]],
                
                Anova(lmer(coloc.traj.mean ~ label + (1 | cell), ok, REML=FALSE))[[3]],
                
                Anova(lmer(logcoloc.traj.mean ~ label + (1 | cell), ok, REML=FALSE))[[3]],
                
                Anova(lmer(arealist ~ label + (1 | cell), ok, REML=FALSE))[[3]],
                
                Anova(lmer(logarealist ~ label + (1 | cell), ok, REML=FALSE))[[3]],
                
                Anova(lmer(sumTrajDistance ~ label + (1 | cell), ok, REML=FALSE))[[3]],
                
                Anova(lmer(log.meanSpeedsperTraj.sumTrajDistance. ~ label + (1 | cell), ok, REML=FALSE))[[3]],
                
                Anova(lmer(trajDuration ~ label + (1 | cell), ok, REML=FALSE))[[3]],
                
                Anova(lmer(log.meanSpeedsperTraj.trajDuration. ~ label + (1 | cell), ok, REML=FALSE))[[3]])


p.bonf.rv   <- p.adjust(p.unadj.rv, method = "bonferroni")
p.hommel.rv <- p.adjust(p.unadj.rv, method = "hommel")

# Compare the results
DFadjustedP.rv   <- data.frame(Unadjusted = p.unadj.rv,
                               Bonferroni = p.bonf.rv,
                               Hommel     = p.hommel.rv)
rownames(DFadjustedP.rv) <- c("outspeed",
                              "outlogspeed",
                              "outdistsperTraj",
                              "outlogdistsperTraj",
                              "outdistsperFrame",
                              "outlogdistsperFrame",
                              "outcoloc",
                              "outlogcoloc",
                              "outarea",
                              "outlogarea",
                              "outTrajLength",
                              "outlogTrajLength",
                              "outTrajDuration",
                              "outlogTrajDuration")
DFadjustedP.rv

#here we plot the graphs for the revamped (rv) stats

ms.plot.rv<-plotFunction( ok, ok$meanSpeeds,"meanSpeeds", "outspeed","Mean speed per trajectory (µm/s)",DFadjustedP.rv)
ggplot(ok, aes(x=meanSpeeds, fill= label)) + geom_density(alpha = 0.2)

#logspeed
lms.plot.rv<-plotFunction( ok, ok$log.meanSpeedsperTraj.meanSpeeds.,"log.meanSpeedsperTraj.meanSpeeds.", "outlogspeed","log mean speed per trajectory (µm/s)",DFadjustedP.rv)
lms.hist.rv<-ggplot(ok, aes(x=log.meanSpeedsperTraj.meanSpeeds., fill= label)) + geom_density(alpha = 0.2)


#mindists PER TRAJ
mdpt.plot.rv<-plotFunction( ok, ok$minimumDistMeanPerTraj,"minimumDistMeanPerTraj", "outdistsperTraj","mean intermito min. distances per traj (µm)",DFadjustedP.rv)
ggplot(ok, aes(x=minimumDistMeanPerTraj, fill= label)) + geom_density(alpha = 0.2)


#logmindists PER TRAJ
lmdpt.plot.rv<-plotFunction( ok, ok$logminimumDistMeanPerTraj,"logminimumDistMeanPerTraj", "outlogdistsperTraj","logmean intermito min. distances per traj (µm)",DFadjustedP.rv)
ggplot(ok, aes(x=logminimumDistMeanPerTraj, fill= label)) + geom_density(alpha = 0.2)


#mindists PER FRAME
mdpf.plot.rv<-plotFunction( ok.n[[2]], ok.n[[2]]$min.dist.everymito_perframe,"min.dist.everymito_perframe", "outdistsperFrame","mean intermito min. distances per frame (µm)",DFadjustedP.rv)
ggplot( ok.n[[2]], aes(x=min.dist.everymito_perframe, fill= label)) + geom_density(alpha = 0.2)


#logmindists PER FRAME
lmdpf.plot.rv<-plotFunction( ok.n[[2]], ok.n[[2]]$log.min.dist.everymito_perframe,"log.min.dist.everymito_perframe", "outlogdistsperFrame","logmean intermito min. distances per frame (µm)",DFadjustedP.rv)
ggplot( ok.n[[2]], aes(x=log.min.dist.everymito_perframe, fill= label)) + geom_density(alpha = 0.2)


#colocTime
ct.plot.rv<-plotFunction( ok, ok$coloc.traj.mean,"coloc.traj.mean", "outcoloc","mean colocalisation times per trajectory (frames)",DFadjustedP.rv)
ggplot(ok, aes(x=coloc.traj.mean, fill= label)) + geom_density(alpha = 0.2)


#logcolocTime
lct.plot.rv<-plotFunction( ok, ok$logcoloc.traj.mean,"coloc.traj.mean", "outlogcoloc","log mean colocalisation times per trajectory (frames)",DFadjustedP.rv)
ggplot(ok, aes(x=logcoloc.traj.mean, fill= label)) + geom_density(alpha = 0.2)


#area
a.plot.rv<-plotFunction( ok, ok$arealist,"arealist", "outarea","area covered per trajectory (µm^2)",DFadjustedP.rv)
ggplot(ok, aes(x=arealist, fill= label)) + geom_density(alpha = 0.2)


#logarea
la.plot.rv<-plotFunction( ok, ok$logarealist,"logarealist", "outlogarea","log area covered per trajectory (µm^2)",DFadjustedP.rv)
ggplot(ok, aes(x=logarealist, fill= label)) + geom_density(alpha = 0.2)


#trajectory distance travelled µm^2
tdi.plot.rv<-plotFunction( ok, ok$sumTrajDistance,"sumTrajDistance", "outTrajLength","2D distance trajectory (µm)",DFadjustedP.rv)
ggplot(ok, aes(x=sumTrajDistance, fill= label)) + geom_density(alpha = 0.2)


#trajectory distance travelled µm^2
ltdi.plot.rv<-plotFunction( ok, ok$log.meanSpeedsperTraj.sumTrajDistance.,"log.meanSpeedsperTraj.sumTrajDistance.", "outlogTrajLength","log 2D distance trajectory (µm)",DFadjustedP.rv)
ggplot(ok, aes(x=log.meanSpeedsperTraj.sumTrajDistance., fill= label)) + geom_density(alpha = 0.2)

#trajectory duration in frames
tdu.plot.rv<-plotFunction( ok, ok$trajDuration,"trajDuration", "outTrajDuration","trajectory duration in frames",DFadjustedP.rv)
ggplot(ok, aes(x=trajDuration, fill= label)) + geom_density(alpha = 0.2)

#trajectory duration in frames
ltdu.plot.rv<-plotFunction( ok, ok$log.meanSpeedsperTraj.trajDuration.,"log.meanSpeedsperTraj.trajDuration.", "outlogTrajDuration","log trajectory duration in frames",DFadjustedP.rv)
ggplot(ok, aes(x=log.meanSpeedsperTraj.trajDuration., fill= label)) + geom_density(alpha = 0.2)

histFunction<-function(dataframe, specificdata, xaxis){
  ggplot(dataframe, aes(x=specificdata, fill= label)) + geom_density(alpha = 0.2) + xlab(xaxis)+
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=15))
}

lms.hist.rv<-histFunction(ok, ok$log.meanSpeedsperTraj.meanSpeeds.,"log mean speed per trajectory (µm/s)")
lct.hist.rv<-histFunction(ok, ok$logcoloc.traj.mean,"log mean colocalisation times per trajectory (frames)")
lmd.hist.rv<-histFunction(ok, ok$logarealist,"log area covered per trajectory (µm^2)")
lmdpert.hist.rv<-histFunction(ok, ok$logminimumDistMeanPerTraj,"logmean intermito min. distances per traj (µm)")
lmdperf.hist.rv<-histFunction( ok.n[[2]], ok.n[[2]]$log.min.dist.everymito_perframe,"logmean intermito min. distances per frame (µm)")
la.hist.rv<-histFunction(ok, ok$logarealist,"log area per trajectory (µm^2)")
ltdi.hist.rv<-histFunction(ok, ok$log.meanSpeedsperTraj.sumTrajDistance.,"log 2D distance per trajectory (µm)")
ltdu.hist.rv<-histFunction(ok, ok$log.meanSpeedsperTraj.trajDuration.,"log trajectory duration (frames)")


ggsave(paste(directory,"plot-logmeanspeed-",strategyName,"-rv",".pdf", sep=""),
       lms.plot.rv, units = "cm", width=18,height=12)
ggsave(paste(directory,"hist-logmeanspeed-",strategyName,"-rv",".pdf", sep=""),
       lms.hist.rv, units = "cm", width=18,height=12)

ggsave(paste(directory,"plot-logarea-",strategyName,"-rv",".pdf", sep=""),
       la.plot.rv, units = "cm", width=18,height=12)
ggsave(paste(directory,"hist-logarea-",strategyName,"-rv",".pdf", sep=""),
       la.hist.rv, units = "cm", width=18,height=12)

ggsave(paste(directory,"plot-logintermitoperframe-",strategyName,"-rv",".pdf", sep=""),
       mdpf.plot.rv, units = "cm", width=18,height=12)
ggsave(paste(directory,"hist-logintermitoperframe-",strategyName,"-rv",".pdf", sep=""),
       lmdperf.hist.rv, units = "cm", width=18,height=12)

ggsave(paste(directory,"plot-logcoloc-",strategyName,"-rv",".pdf", sep=""),
       lct.plot.rv, units = "cm", width=18,height=12)
ggsave(paste(directory,"hist-logcoloc-",strategyName,"-rv",".pdf", sep=""),
       lct.hist.rv, units = "cm", width=18,height=12)

ggsave(paste(directory,"plot-logtrajduration-",strategyName,"-rv",".pdf", sep=""),
       ltdu.plot.rv, units = "cm", width=18,height=12)
ggsave(paste(directory,"hist-logtrajduration-",strategyName,"-rv",".pdf", sep=""),
       ltdu.hist.rv, units = "cm", width=18,height=12)

ggsave(paste(directory,"plot-logtrajdistance-",strategyName,"-rv",".pdf", sep=""),
       ltdi.plot.rv, units = "cm", width=18,height=12)
ggsave(paste(directory,"hist-logtrajdistance-",strategyName,"-rv",".pdf", sep=""),
       ltdi.hist.rv, units = "cm", width=18,height=12)

#all the other histograms, taken from code above
ggplot(ok, aes(x=meanSpeeds, fill= label)) + geom_density(alpha = 0.2)
ggplot(ok, aes(x=minimumDistMeanPerTraj, fill= label)) + geom_density(alpha = 0.2) 
ggplot(ok, aes(x=coloc.traj.mean, fill= label)) + geom_density(alpha = 0.2)
ggplot(ok, aes(x=arealist, fill= label)) + geom_density(alpha = 0.2)
ggplot(ok, aes(x=sumTrajDistance, fill= label)) + geom_density(alpha = 0.2)
ggplot(ok, aes(x=log.meanSpeedsperTraj.sumTrajDistance., fill= label)) + geom_density(alpha = 0.2)
ggplot(ok, aes(x=trajDuration, fill= label)) + geom_density(alpha = 0.2)
ggplot(ok, aes(x=log.meanSpeedsperTraj.trajDuration., fill= label)) + geom_density(alpha = 0.2)



}
