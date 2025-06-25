# Physical genetic dynamics of mitochondria in plant cells

This repo is for the analysis of confocal time lapse of any organelle, here mitochondria, alongside other functional biological data, here mtDNA presence/absence. 

Required R packages\
sp\
plotrix\
TrackMateR\
dplyr\
igraph\
brainGraph\
ggpubr\
ggplot2\
lme4\
Matrix\
car\
RColorBrewer\
tidyr\
ggplot2\


Following imaging in one plane and tracking of mitochondria using Trackmate, trackmate files in .xml format will be generated. 
The scripts use the TrackmateR importer (https://quantixed.github.io/TrackMateR/)

The first set of scripts define mtDNA puncta based on Trackmate output contrast:
```sh
contrastThresholdPassing.R
```

or Maxima via choosing a maxima prominence in imageJ, unit switching, and correlating maxima points to trackmate output:\
Find maxima
```sh
getMaxima.ijm
```
Processing
```sh
pixelCoordsToMicronCoords.R
usingMaxima.R
```
Defining trajectories via Maxima
```sh
maximaThresholdPassing.R
```

Once lists of trajectories and their mtDNA presence/absence has been defined, we take these forward to analyse.
These two scripts are the bulk of the analysis.

```sh
PhysicalStatisticsBymtDNADefinition.R
NetworkAnalysisBymtDNADefinition.R
```

Plots relating to quantifications of fixed cells and mitochondrial populations with or without mtDNA can be found here:
This uses output tables from non-automated segmentation of confocal z-stack images. 
```sh
QuantificationFixedSamplesFromImaris.R
```

Supplementary plots for rosette area and root lengths can be found here:
```sh
RosetteAndRootPlots.R 
```
And for oxygen consumption measurements, here:
```sh
OxygenPlotting.R
```
