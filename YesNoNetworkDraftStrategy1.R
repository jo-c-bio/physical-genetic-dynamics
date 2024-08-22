
#need to change the files in the directories to have this name
#vidList<-rep(paste("vid",1:10,sep=""))
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


#These can all be bought in as args, like allfiles in a directory, but all need corresponding contrast values.
#if bringing in filename as entire files directory, will need to alter that below

library(TrackMateR)
library(dplyr)
library(igraph)
library(brainGraph)
library(ggplot2)
library(ggpubr)
#install.packages("lme4")
#install.packages("Matrix")
library(lme4)
library(Matrix)
library(car)
library(RColorBrewer)


directory<-"/Users/joannachustecki/Documents/PostDoc23-Data/nucleoidQuantification/currentDecentTimelapseSYBR/cropped/retracked28-2-24/"

#empty dataframes to collect stats for over all videos
overallDegreeYesNoframe<-c()
overallLCCYesNoFrame<-c()
degree1sCountOverFrames<-c()
degree0sCountOverFrames<-c()
degree3sCountOverFrames<-c()
totalTypeNodesOverFrames<-c()

readDataframe <- function(filename) {
  tmObj <-readTrackMateXML(XMLpath = paste(directory, filename, sep = ""))
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
  return(df)
}


#socialNetworkGeneration<-function(){
 for(v in 1:length(vidList)){
#for(v in 10){
  print(v)
  filename<- vidList[v]
  contrastThreshold <- contrastList[v]
  directory<-"/Users/joannachustecki/Documents/PostDoc23-Data/nucleoidQuantification/currentDecentTimelapseSYBR/cropped/retracked28-2-24/"
  df <- readDataframe(filename)
  #read in the mtDNA presence/absence categories
  
  #Strategy1<- read.csv( paste(directory, "trajsWithAtLeastOnemtDNA_Contrast_",filename,".csv",sep="")) 
  #Strategy2a<- read.csv( paste(directory, "trajsatLeastOneAdjacentmtDNA_Contrast_",filename,".csv",sep=""))
  #Strategy2b<- read.csv( paste(directory, "trajsatLeastThreeAdjacentmtDNA_Contrast_",filename,".csv",sep=""))
  #Strategy1.m<- read.csv( paste(directory, "trajsWithAtLeastOnemtDNA_Maxima_",filename,".csv",sep="")) 
  Strategy2a.m<- read.csv( paste(directory, "trajsatLeastOneAdjacentmtDNA_Maxima_",filename,".csv",sep=""))
  #Strategy2b.m<- read.csv( paste(directory, "trajsatLeastThreeAdjacentmtDNA_Maxima_",filename,".csv",sep=""))
  #StrategyList<-list(Strategy1,Strategy2a,Strategy2b,Strategy1.m,Strategy2a.m,Strategy2b.m)
  #StrategyListNames<-c("Strategy1","Strategy2a","Strategy2b","Strategy1.m","Strategy2a.m","Strategy2b.m")
  StrategyList<-list(Strategy2a.m)
  StrategyListNames<-c("Strategy2a.m")
  # overall stats for this dataset
  maxt = max(df$t)
  ntraj = max(df$traj)
  #we'll take the difference in secs between the first two values as the frametime, it should be consistent over the video. Also round.
  frameTime= df$secs[2]-df$secs[1]
  #Here we impose a colocalisation threshold for what makes an interaction, in microns
  threshold = 1.6
  
  #which frames do we want figures output for?
  frames.to.plot = NULL
  #These commands are for when we automate the process.
  # if(length(args) > 2) {
  #   for(i in 3:length(args)) {
  #    frames.to.plot = c(frames.to.plot, as.numeric(args[i]))
  #  }
  # }
  frames.to.plot = c(1, 10,20, 50 ,80, 109)
  

  # initialise data structures. adj matrix and list, matrix for colocalisation times, data frame to store mean min distances, list to store number of vertices at each frame
  amlist = data.frame(firstframe = NULL, t1 = NULL, t2 = NULL)
  am = matrix(nrow = ntraj, ncol = ntraj)
  colocal.time = matrix(0, nrow = ntraj, ncol = ntraj)
  dist.frame = data.frame(frame = NULL, mean.min.dist = NULL)
  vlist = NULL
  mylist = list()
  
  # loop through frames
  message("Computing distances...")
  for(t in 1:maxt) {
    #Jc adding: empty all dists list
    allDists<-c()
    # store total number of trajectories recorded up to this frame (to account for singletons)
    vlist = c(vlist, length(unique(df$traj[df$t <= t])))
    mylist <- append(mylist, list(unique(df$traj[df$t <= t])))
    
    # grab trajectories from this frame
    subset = df[df$t == t,]
    min.dists = NULL
    # loop through pairs of trajectories
    for(i in 1:nrow(subset)) {
      min.dist.2 = -1
      for(j in 1:nrow(subset)) {
        if(i != j) {
          # compute distance for this pair -- recording minimum dists
          pair.dist.2 = (subset$x[i]-subset$x[j])**2 + (subset$y[i]-subset$y[j])**2
          if(min.dist.2 == -1 | pair.dist.2 < min.dist.2) { min.dist.2 = pair.dist.2 }
          # if this distance is below our threshold
          if(pair.dist.2 < threshold**2) {
            # if we haven't recorded this pair yet, do so
            if(is.na(am[subset$traj[i],subset$traj[j]])) {
              
              amlist = rbind(amlist, data.frame(firstframe = t, t1=subset$traj[i], t2=subset$traj[j]))
              am[subset$traj[i],subset$traj[j]] = am[subset$traj[j],subset$traj[i]] = t
              
        
            }
            # add colocalisation time to matrix
            if(i < j) {
              colocal.time[subset$traj[i],subset$traj[j]] = colocal.time[subset$traj[j],subset$traj[i]] = colocal.time[subset$traj[j],subset$traj[i]]+1
            }
          }
        }
      }
      # summarise distance statistics
      min.dists = c(min.dists, sqrt(min.dist.2))
    }
    dist.frame = rbind(dist.frame, data.frame(frame = t, mean.min.dist = mean(min.dists)))
    
  }
  
  # output the adjacency list. a fake self-edge on the highest-labelled node is appended as the final entry to characterise the number of nodes
  output.filename = paste(directory, filename, "-amlist.csv", sep="")
  amlist.out = rbind(amlist, data.frame(firstframe=0, t1=length(unique(df$traj)), t2=length(unique(df$traj))))
  write.csv(amlist.out, output.filename, row.names=FALSE)
  
  
  message("Outputting trajectory statistics...")
  traj.frame = data.frame(traj.time = NULL, traj.degree = NULL, traj.rate = NULL)
  # loop through individual trajectories, reporting length of trajectory and degree (number of encounters)
  for(traj in 1:ntraj) {
    subset = df[df$traj == traj,]
    traj.frame = rbind(traj.frame, data.frame(traj.time = nrow(subset), traj.degree = length(which(am[traj,] != 0)), traj.rate = length(which(am[traj,] != 0))/nrow(subset)))
  }
  
  # initialise data frame to report statistics
  stats.frame = data.frame(frame = NULL, mean.min.dist = NULL, mean.degree = NULL, sd.degree = NULL, max.degree = NULL, diameter = NULL, betweenness = NULL, efficiency = NULL, num.edges = NULL, num.vertices = NULL, cc.num = NULL, cc.mean.size = NULL, cc.max.size = NULL, singletons = NULL)
  
  # loop through frames
  message("Calculating network statistics...")
  
  #The whole same thing but for .mtDNA graph
  for(frame.number in 1:maxt) {
    # subset out the interactions recorded up to this frame
    amlist.frame = amlist[amlist$firstframe <= frame.number,]
    # construct graph structure from edge list. coerce to character labels first, otherwise igraph fills in any missing numeric labels -- which we don't want, as these may correspond to trajectories that have yet to appear
    edgelist.frame = as.matrix(data.frame(t1=as.character(amlist.frame$t1), t2=as.character(amlist.frame$t2)))
    g = graph_from_edgelist(edgelist.frame, directed=F)
    # figure out if we need to add any trajectories that exist, but haven't interacted, so don't appear in the edge list 
    
    to.add = vlist[frame.number]-length(V(g))
    #g = add.vertices(g, to.add)
    # ^ we won;t just add n number of random vertices, but labelled ones (code below)
    
    #instead of just adding 'NA' to the graph for singletons, let's use their real traj names
    #What are the unique trajectories of the graph as it stands?
    p<- c(edgelist.frame[,1],edgelist.frame[,2])
    a<-as.numeric(unique(p))
    #vlist aka What are the total number of trajectories recorded up to this frame (to account for singletons)?
    b<-as.numeric(mylist[[frame.number]] )
    b[!(b %in% a)]
    if(length(b[!(b %in% a)]) !=  to.add){print("error in singletons nomination")}
    
    #add labelled singletons
    g <- g + vertex(b[!(b %in% a)])
    
    # add statistics to data frame
    g.frame = data.frame(frame = frame.number, mean.min.dist = dist.frame$mean.min.dist[frame.number], mean.degree = mean(degree(g)), sd.degree = sd(degree(g)), max.degree = max(degree(g)), diameter = diameter(g), betweenness = mean(betweenness(g)), efficiency = efficiency(g, type="global"), num.edges = length(E(g)), num.vertices = length(V(g)), cc.num = components(g)$no, cc.mean.size = mean(components(g)$csize), cc.max.size = max(components(g)$csize), singletons = sum(degree(g)==0))
    stats.frame = rbind(stats.frame, g.frame)
      
      #let's colour the graph by mtDNA presence/absence. This does not alter stats table, 
      #that has already been written
      for(p in 1:length(StrategyList)){
        
        s<-StrategyList[[p]]
        strategyListName<- StrategyListNames[p]
        colnames(s)<- c("X","trajsInOrder", "count")
        #This match function goes through the strategy dataframe, find the nodes names the graph at that frametime uses,
        # and prints out the corresponding trajectory data from strategy
        #we need to sort out the singletons- they need presence/absence too. 
        
        
        #this is for main network
        o<-data.frame(s[match(V(g)$name , s$trajsInOrder),], V(g)$name)
        trueSubset =  o[o[,3]==TRUE,2]
        trueSubset = trueSubset[!is.na(trueSubset)]
        falseSubset =  o[o[,3]==FALSE,2]
        falseSubset = falseSubset[!is.na(falseSubset)]
        #print(falseSubset)
        
        g <- set.vertex.attribute(
          graph = g,
          name = "yesno",
          #must use matching here, or the listed 1,2,3...n nodes will be weirdly labelled by out of numeric order traj names.
          index = match(falseSubset,V(g)$name),
          value = 0
        )
        g <- set.vertex.attribute(
          graph = g,
          name = "yesno",
          index =  match(trueSubset,V(g)$name),
          value = 1
        )
        
        
        # if we want to output the structure at this frame, do so
        if(frame.number %in% frames.to.plot | frame.number == maxt) {
        pdf(paste(directory,filename,"-plot-", frame.number, "_colourByPresence_",strategyListName,".pdf", sep=""))
        plot(g, vertex.color = V(g)$yesno,vertex.size=3, vertex.label=NA) #vertex.label=V(g)$name
        dev.off()
        
        #fiddling about with painting by degree
        gx<-g
        V(gx)$degree <- degree(gx)
        uni_all <- seq( min(V(gx)$degree), max(V(gx)$degree))
        colors <- data.frame( color = hcl.colors(length(uni_all), rev = T),
                              levels = uni_all)
        # Use match to get the index of right color, matching on levels
        V(gx)$color <- colors$color[match(V(gx)$degree, colors$levels)]
        # use degree as labels
        V(gx)$xoro <- ifelse(V(g)$yesno == "0", ".", "x")
        V(gx)$label.cex = 0.8
        pdf(paste(directory,filename,"-plot-", frame.number, "_colourByDegree_",strategyListName,".pdf", sep=""))
        plot(gx, vertex.color = V(gx)$color ,vertex.size=4, vertex.label= V(gx)$xoro) #vertex.label=V(g)$name
        dev.off()
        }
        #now let's get a degree comparison frame for each plot
        degreeYesNoframe<-data.frame(rep(frame.number,length(degree(g))),
                                     degree(g), 
                                     V(g)$yesno,
                                     rep(letters[v],length(degree(g)))) #add cell lettering to each also, so we can keep track of which cell thse values refer to.
        
        #^degree and attribute are printed out in the same order. but if not sure, can check labelled network image to validate
        #plot(g, vertex.color = V(g)$yesno,vertex.size=3, vertex.label.cex = .7)
        
        #Then let's get a local clustering coefficient
        LCCYesNoFrame<- data.frame(rep(frame.number,length(transitivity(g, type = "local"))),
                                   transitivity(g, type = "local"), 
                                   V(g)$yesno,
                                   log(transitivity(g, type = "local")),
                                   rep(letters[v],length(transitivity(g, type = "local")))) #add cell lettering to each also, so we can keep track of which cell thse values refer to.
        
        
        #bind these stats table in order to made a table to run lmer across frames of interest, over all videos. 
        overallDegreeYesNoframe <- rbind(overallDegreeYesNoframe,degreeYesNoframe)
        overallLCCYesNoFrame<- rbind(overallLCCYesNoFrame,LCCYesNoFrame)
        
        
        yMax <- 35
        #Be careful: magic number here in y axis limits^
        if(max(degreeYesNoframe$degree.g.) > yMax){print("alter y max")}
  
        degcompplot<- ggboxplot(degreeYesNoframe, x = "V.g..yesno", y = "degree.g.", color = "V.g..yesno") +  
                      geom_jitter(width = 0.2,aes(color = V.g..yesno)) +  stat_compare_means() +
                      theme(text = element_text(size = 10)) + ylim(-1, yMax)
        LCCcompplot<- ggboxplot(LCCYesNoFrame, x = "V.g..yesno", y = "transitivity.g..type....local..", color = "V.g..yesno") +  
                      geom_jitter(aes(color = V.g..yesno)) +stat_compare_means(vjust=-0.6,hjust=1)
        # if we want to output the plot at this frame, do so
      if(frame.number %in% frames.to.plot | frame.number == maxt) {   
        ggsave(paste(directory,filename,"-plot-", frame.number, "_Deg_",strategyListName,".pdf", sep=""),
               degcompplot, units = "cm")
        ggsave(paste(directory,filename,"-plot-", frame.number, "_LCC_",strategyListName,".pdf", sep=""),
               LCCcompplot, units = "cm")
      }
        
        
        #let's test the number of nodes that have degree 1 between +ve and -ve nodes for mtDNA. 
        #This data frame builds over frame time. 
        degree1s <- degreeYesNoframe %>% filter(degree.g. == 1)
        degree0s <- degreeYesNoframe %>% filter(degree.g. == 0)
        degree3s <- degreeYesNoframe %>% filter(degree.g. == 3)
        totalTypeNodes <-  degreeYesNoframe %>% count(V.g..yesno)
        totalTypeNodes <-cbind(totalTypeNodes,rep(unique(degree1s$rep.frame.number..length.degree.g...),2))
        totalTypeNodes <- cbind(totalTypeNodes,rep(letters[v],2))
        #degree1sCountFrame<-cbind(degree1sCount,rep(unique(degree1s$rep.frame.number..length.degree.g...),2))
        
        totalTypeNodesOverFrames<- rbind(totalTypeNodesOverFrames, totalTypeNodes)
        degree1sCountOverFrames <- rbind(degree1sCountOverFrames, degree1s)
        degree0sCountOverFrames <- rbind(degree0sCountOverFrames, degree0s)
        degree3sCountOverFrames <- rbind(degree3sCountOverFrames, degree3s)
      
        
      }
    
  }
  
  # output stats table
  output.filename = paste(directory,filename, "-stats.csv", sep="")
  write.csv(stats.frame, output.filename, row.names=FALSE)
  
}
#}
#!!!

colnames(degree1sCountOverFrames) <-c("frame","degree","yesorno","cell")
colnames(degree0sCountOverFrames) <-c("frame","degree","yesorno","cell")
colnames(degree3sCountOverFrames) <-c("frame","degree","yesorno","cell")
colnames(totalTypeNodesOverFrames) <- c("yesorno","nodeN","frame","cell")

#let's plot out number of nodes with degree 1 and then number of nodes with degree 1  
#proportional to the number of nodes in that frame with 0 or 1 status.


specificDegreesFunction<- function(degSpecificDF, degreeNo,MN){
  allCellsdf3Proportional<-c()
  #for each cell
  for(c in 1:length(vidList)){
    celln <- letters[c]
    file <- vidList[c]
    #get degree and total node type in that cell
    df <-  degSpecificDF %>% filter(cell==celln)
    df2 <-  totalTypeNodesOverFrames %>% filter(cell==celln)
    degreeXplot<-ggplot(df, aes(x = as.factor(frame), fill = as.factor(yesorno))) +
            geom_bar(position = "dodge")+labs(y = paste("Nodes of degree",degreeNo,sep=""), x = "Frame", title = paste(file))
    
    ggsave(paste(directory,file,"-plot-Degreeof",degreeNo,"_",strategyListName,".pdf", sep=""),
           degreeXplot, units = "cm")
    
    #let's make these values proportional to the number of nodes in that frame with 0 or 1 mtdna status. 
    degreeXsCount <- df %>%  group_by(frame) %>% count(degree,yesorno)
    nodesCount <- df2 %>%  group_by(frame)
    
    #this is a specific bug- if frame= 1 sometimes the network has not build up far enough
    #to include even one node that is mtDNA -ve. So we remove the first two frames.
    print(nrow(degreeXsCount))
    print(nrow(nodesCount))
   #if you don;t have double frame entries in both data sets, remove and that are only singles!
    if( nrow(degreeXsCount)-length(unique(degreeXsCount$frame))  != length(unique(degreeXsCount$frame))){
      toremove<- length(unique(degreeXsCount$frame))  - (nrow(degreeXsCount)-length(unique(degreeXsCount$frame))) + MN
      toremove<- c(1:toremove)
      print("need")
      x<-degreeXsCount$frame
      y<-nodesCount$frame
      
      SameElements <- function(a, b) return(identical(sort(a), sort(b)))
      SameElements(x,y)
      
  
      nodesCount <- nodesCount[ ! nodesCount$frame %in% toremove, ]
      degreeXsCount <- degreeXsCount[ ! degreeXsCount$frame %in% toremove, ]
  
      print(nrow(degreeXsCount))
      print(nrow(nodesCount))
      print("okay")
    }
    df3<-cbind(degreeXsCount,nodesCount)
    proportions<- data.frame(df3$cell,df3$frame...1,df3$yesorno...3,df3$n/df3$nodeN)
    colnames(proportions)<-c("cell","frame","yesorno","proportion")
    
    allCellsdf3Proportional <-rbind(allCellsdf3Proportional,proportions)
    
    degreeXProportionplot<-ggplot(proportions, aes(x = as.factor(frame), y= proportion, fill = as.factor(yesorno))) +
      geom_col(position = "dodge")+labs(y = paste("Nodes of degree",degreeNo), x = "Frame", title = paste(file))+
      labs(y = "Proportion", x = "Frame",
           title = paste(filename,"-", s,"\n","Nodes of degree ",degreeNo,", proportional to number of nodes per type",sep=""))
  
      ggsave(paste(directory,file,"-plot-Degreeof",degreeNo,"Proportional_",strategyListName,".pdf", sep=""),
             degreeXProportionplot, units = "cm")
    
  }
  return(allCellsdf3Proportional)
}

#including MN, to remove n frames that don't generate because of a lack of nodes at that degree value. 
degree1sProportionalCountOverFrames <- specificDegreesFunction(degree1sCountOverFrames , "1", 1)
degree0sProportionalCountOverFrames <- specificDegreesFunction(degree0sCountOverFrames , "0", 5)
degree3sProportionalCountOverFrames <- specificDegreesFunction(degree3sCountOverFrames , "3", 19)
#^must add big magic no.

#how to plot these nicely?
degree0sProportionalCountOverFrames$cellframe <-  paste(degree0sProportionalCountOverFrames$cell,
                                                        degree0sProportionalCountOverFrames$frame,
                                                        sep="")
slimmedDown0s<-degree0sProportionalCountOverFrames[ degree0sProportionalCountOverFrames$frame %in% c(20,50,80,109), ]

#use the forcats package to force sensible x axis labels
slimmedDown0plot<- 
  ggplot(slimmedDown0s, aes(x = fct_inorder(as.factor(cellframe)), y= proportion, fill = as.factor(yesorno))) +
  geom_col(position = "dodge") + labs(y = "Proportion", x = "Frame",
                                      title = paste("Nodes of degree 0, proportional to number of nodes per type",sep=""))
slimmedDown0plot
ggsave(paste(directory,"plot-Degreeof0Proportional_",strategyListName,".pdf", sep=""),
       slimmedDown0plot, width = 30,
       height = 10, units = "cm")


#same but for degree 1s.
degree1sProportionalCountOverFrames$cellframe <-  paste(degree1sProportionalCountOverFrames$cell,
                                                        degree1sProportionalCountOverFrames$frame,
                                                        sep="")
slimmedDown1s<-degree1sProportionalCountOverFrames[ degree1sProportionalCountOverFrames$frame %in% c(20,50,80,109), ]

#use the forcats package to force sensible x axis labels
slimmedDown1plot<-
  ggplot(slimmedDown1s, aes(x = fct_inorder(as.factor(cellframe)), y= proportion, fill = as.factor(yesorno))) +
  geom_col(position = "dodge") + labs(y = "Proportion", x = "Frame",
                                      title = paste("Nodes of degree 1, proportional to number of nodes per type",sep=""))
slimmedDown1plot
ggsave(paste(directory,"plot-Degreeof",1,"Proportional_",strategyListName,".pdf", sep=""),
       slimmedDown1plot,  width = 30,
       height = 10,units = "cm")


#same but for degree 3s.
degree3sProportionalCountOverFrames$cellframe <-  paste(degree3sProportionalCountOverFrames$cell,
                                                        degree3sProportionalCountOverFrames$frame,
                                                        sep="")
slimmedDown3s<-degree3sProportionalCountOverFrames[ degree3sProportionalCountOverFrames$frame %in% c(20,50,80,109), ]

#use the forcats package to force sensible x axis labels
slimmedDown3plot<-
ggplot(slimmedDown3s, aes(x = fct_inorder(as.factor(cellframe)), y= proportion, fill = as.factor(yesorno))) +
  geom_col(position = "dodge") + labs(y = "Proportion", x = "Frame",
                                      title = paste("Nodes of degree 3, proportional to number of nodes per type",sep=""))

slimmedDown3plot
ggsave(paste(directory,"plot-Degreeof",3,"Proportional_",strategyListName,".pdf", sep=""),
       slimmedDown3plot,  width = 30,
       height = 10, units = "cm")


#build, for each video, a table of degree and LCC values. 
#However, this is not yet running over all strategy types. 
#BUT we don;t want it to be - we need to pick one and move on. 
#chosen MAxima use- strategy 2a- "at least one adjacent"

#we want to do a linear mixed effect model across all the cells.
#we will start with degree and LCC, for networks built until the final frame of the video. 

#first, rename columns for easier filtering
colnames(overallDegreeYesNoframe)<- c("frame","stat","mtDNA","cell")
colnames(overallLCCYesNoFrame)<- c("frame","LCC","mtDNA","stat","cell")

#lets' plot these all data with averages, 
#comparisonPlot fucntion with EXPORT an image too.
comparisonPlot<- function(inputData, anovaInput.pvalue, plotType,statLabel){
  pVal<- anovaInput.pvalue
  
  cplot<-ggplot(data = inputData ,aes(x=mtDNA, y= stat, colour = cell)) +
    geom_jitter(alpha=0.35)+
    stat_summary(
      geom = "point",
      fun = "mean",
      col = "white",
      size = 3,
      shape = 23,
      fill = "red",
      alpha=0.6
    ) + 
    ggtitle(paste("Anova output= Pr(>ChiSq) =",pVal, "\n",
                  "n0 = ",count(inputData[inputData$mtDNA == 0,]),",",
                  "n1 = ",count(inputData[inputData$mtDNA == 1,]),sep=" "))+ 
    ylab(statLabel)+
    theme(
      legend.text = element_text(size = 12), 
      legend.title = element_text(size = 14),
      axis.text=element_text(size=12),
      axis.title=element_text(size=14)
    )
  
  cplot
  
   ggsave(paste(directory,"AllVideos-plot-", plotType, "_Deg_",statLabel,".pdf", sep=""),
        cplot, units = "cm", width= 20, height=14
      )
  
  
}

#for chosen frame (Filtered)

install.packages("lme4")
library(lme4)
library(ggplot2)
library(car)
library(dplyr)

frameFiltered <- 109
comparisonPlot(inputData = overallDegreeYesNoframe %>% filter(frame == frameFiltered), 
               #old anova method
                anovaInput = Anova( lmer(stat ~ mtDNA + (1 | cell), 
                                       overallDegreeYesNoframe %>% filter(frame == frameFiltered), 
                                        REML=FALSE) )[[3]], 
        
               plotType = paste("FilteredFrame",frameFiltered,sep=""),
               statLabel= "Node Degree")      

#for all frames (Overall)- gives very small P-value. 
comparisonPlot(inputData = overallDegreeYesNoframe, 
                  #old anova method
                   anovaInput = Anova( lmer(stat ~ mtDNA + (1 | cell), 
                                           overallDegreeYesNoframe , 
                                            REML=FALSE) )[[3]], 
          
               plotType = "AllFrames",
               statLabel= "Node Degree")

fit<- lmer(stat ~ mtDNA + (1 | cell), 
                   overallDegreeYesNoframe, 
                   REML=FALSE) 
hist(residuals(fit))
qqnorm(residuals(fit))
qqline(residuals(fit))

#these give you the slope and intercept of every individual
plot(coef(outdegreeF)$cell)
plot(coef(outdegreeO)$cell)

frames.to.plot = c(1, 10, 20, 50 ,80, 109)
framesName<-paste(1, 10, 20, 50 ,80, 109)
comparisonPlot(
  inputData = overallDegreeYesNoframe %>% filter(frame %in%  frames.to.plot),
  #old anova method
   anovaInput = signif(Anova(
     lmer(
       stat ~ mtDNA + (1 | cell),
       overallDegreeYesNoframe %>% filter(frame == frameFiltered),
       REML = FALSE
     )
   )[[3]], 4),

  
  plotType = paste("FilteredFrame", framesName, sep = ""),
  statLabel = "Node Degree"
)   


#Now for Local clustering coefficients
#regular
colnames(overallLCCYesNoFrame)<-c("frame", "stat"   ,"mtDNA" ,"logstat" , "cell" )
frameFiltered <- 50
comparisonPlot(inputData = overallLCCYesNoFrame %>% filter(frame == frameFiltered), 
               #old anova method
                anovaInput = Anova( lmer(stat ~ mtDNA + (1 | cell), 
                                         overallDegreeYesNoframe %>% filter(frame == frameFiltered), 
                                         REML=FALSE) )[[3]], 
             
               plotType = paste("FilteredFrame",frameFiltered,sep=""),
               statLabel= "LCC")      

frames.to.plot = c(1, 10, 20, 50 ,80, 109)
framesName<-paste(1, 10, 20, 50 ,80, 109)
comparisonPlot(inputData = overallLCCYesNoFrame  %>% filter(frame %in%  frames.to.plot), 
               #old anova method
                anovaInput = Anova( lmer(stat ~ mtDNA + (1 | cell), 
                                         overallDegreeYesNoframe %>% filter(frame == frameFiltered), 
                                         REML=FALSE) )[[3]], 
             
               plotType = paste("FilteredFrame",framesName,sep=""),
               statLabel= "LCC")      
ggplot(overallLCCYesNoFrame, aes(x=stat, fill= as.factor(mtDNA))) + geom_density(alpha = 0.2)
ggplot(overallLCCYesNoFrame, aes(x=logstat, fill= as.factor(mtDNA))) + geom_density(alpha = 0.2)

#Any difference in number of singletons?
#a singleton has degree=0


singletonsDF<- overallDegreeYesNoframe %>% filter(stat == 0) %>% filter(frame %in%  c(100:109))

ggplot(singletonsDF, aes(x = as.factor(frame), fill = as.factor(mtDNA))) +
  geom_bar(position = "dodge")+labs(y = "Nodes of degree 0", x = "Frame",
                                    title = paste("-"))






