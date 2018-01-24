# Stability Properties Computation for Flowcytometric Data V1.0
# by Florian Centler and Zishu Liu, UFZ Leipzig, Germany
#
# This script accompanies the paper on "Ecological Stability Properties of
# Microbial Communities Assessed by Flow Cytometry" by Liu et al., 2018
# http://msphere.asm.org/content/3/1/e00564-17
# Please cite this work when using this script.
#
# USER INPUT STARTS HERE ======================================================
#
# 1. Specify the directory and filename of your experimental data.
# 
# Data is expected in the format:
# Time Gate1 Gate2 ...
# 0 0.1 0.5 ...
# 1.5 0.2 0.4 ...
# ...
# 
# Columns can be separated by Tab or Space. Gate abundances are relative, i.e.
# summing over all gates at a time point should result in value 1
#
# Output of this script will be stored in "dir/filename.pdf", e.g.,
# "dir/relativeGateAbundanceOverTime.txt.pdf"

filename <- "experiment1/relativeGateAbundanceOverTime.txt"

# 2. Select time points referring to ...
#
#   the start of the experiment (experimentStart),
#   the end of the experiment (experimentEnd),
#   maximal deviation (d_max) using Euclidean distance (overrideMaxEuclidean),
#   maximal deviation (d_max) using Canberra distance (overrideMaxCanberra), and
#   the disturbance event (tref).
#
# A value of -1 will let the script select these time points automatically,
# assuming that the file contains a single disturbance experiment (i.e. start and
# end of the experiment refer to the first and last entry, respectively).
#
# Otherwise, specify a time value which is actually present in the data, i.e.
# a value which appears in the first column of the input file.

experimentStart <- -1
experimentEnd <- -1
overrideMaxEuclidean <- -1
overrideMaxCanberra <- -1
tref <- -1

# To manually select gates to be plotted instead of automatically selecting the
# two most dominant gates at the start and end of the experiment, uncomment the
# following line(s) (remove the '#' character at the beginning of the line) and
# provide gate names as they appear in the column headings in the input file.
#
#plotGatesStart <- c("G1", "G2")
#plotGatesEnd <- c("G3", "G4")
#
# USER INPUT ENDS HERE ========================================================


library("plotrix")
library("vegan")

titleSize <- 0.95

getMembers <- function(mergeTable, idx) {
  if (idx < 0) {
    return(-idx)
  } else {
   
    aList <- getMembers(mergeTable, mergeTable[idx, 1])
    bList <- getMembers(mergeTable, mergeTable[idx, 2])
   
    return(c(aList, bList))
  }
}

getReferencePhaseEnd <- function(data, begin, end) {
  myData <- data[data[1] >= begin & data[1] <= end,]
  rownames(myData) <- myData[[1]]
  myData <- myData[-1]
  myD <- dist(myData, method="canberra")
  hc <- hclust(myD)
  
  plot(hc, cex=0.6, main="Reference state detection", xlab=NA, sub=NA, font.main=1, cex.main=titleSize)
  
  startTime <- min(as.numeric(hc$labels))
  
  index <- match(startTime, hc$labels)
  
  nMerger <- nrow(hc$merge)
  startIdx <- 0
  
  for (i in 1:nMerger) {
    if (hc$merge[i,1] == -index) {
      startIdx <- hc$merge[i,2]
      break
    } else {
      if (hc$merge[i,2] == -index) {
        startIdx <- hc$merge[i,1]
        break
      }
    }
  }
  
  if (startIdx == 0) {
    print("Error!")
  }  
  
  myList <- getMembers(hc$merge, startIdx)
  
  value <- max(as.numeric(hc$labels[myList]))
  
  return(value)
}

data <- read.table(filename, header=TRUE, stringsAsFactors = FALSE)

if (!exists("data[1,1]")) {
  stop(paste("File '", filename, "' could not be read. Aborting.", sep=""))
}

pdf(file=paste(filename, "pdf", sep="."))
par(mfrow=c(3,3))

# compute defaults if user didn't specify values 

if (experimentStart == -1) {
  experimentStart <- data[1,1]
}
if (experimentEnd == -1) {
  experimentEnd <- data[nrow(data),1]
}
if (tref == -1) {
  tref <- getReferencePhaseEnd(data, experimentStart, experimentEnd) # time point of the disturbance, estimated from change in community structure
}
# Characterizing the reference space

numberOfGates <- ncol(data) - 1
data$referencePhase <- FALSE
data[data[1] >= experimentStart & data[1] <= tref, "referencePhase"] <- TRUE

referenceState <- colMeans(data[data$referencePhase==TRUE, 2:(numberOfGates + 1)])
gateSDs <- apply(data[data$referencePhase==TRUE, 2:(numberOfGates + 1)], 2, sd)

# (Radii will be computed later)

# get two most abundant gates
if (!exists("plotGatesStart")) {
  topGates <- names(sort(referenceState, decreasing=TRUE))[1:2]
} else {
  topGates <- plotGatesStart
}
 
if (!exists("plotGatesEnd")) {
  finalTopGates <- names(sort(data[data[1] == experimentEnd, 2:(numberOfGates + 1)], decreasing=TRUE))[1:2]
} else {
  finalTopGates <- plotGatesEnd
}

# Function for plotting community evolution with respect to two gates
plot2Dprojection <- function(gates, title, xlabel, ylabel, ref_radius) {
  # plot states
  plot(data[data[1]>=experimentStart & data[1]<=experimentEnd, gates], col=c("black", "grey")[as.factor(data[data[1]>=experimentStart & data[1]<=experimentEnd,"referencePhase"])],font.main=1, main =title, xlab=xlabel,ylab=ylabel, cex.main=titleSize)
  # add lines to be able to follow dynamic evolution of the state
  lines(data[data[1]>tref & data[1]<=experimentEnd, gates])
  #ZS: give the direction of evolution from reference state to disturbed state
  
  # mark tref
  #points(data[data[1]==tref, gates], pch=16, col="lightgreen")
  # mark final state
  points(data[data[1]==experimentEnd, gates], pch=19, cex=1.3)
  # mark the reference state
  points(referenceState[gates[1]], y=referenceState[gates[2]],  col="red", pch=21, bg="white",cex=1.3)
  # mark states of maximal deviation from reference state
  points(maxDevStateEuclidean[gates[1]], y=maxDevStateEuclidean[gates[2]],  col = "deepskyblue3", pch = 24, cex=1.2,bg="white")
  points(maxDevStateCanberra[gates[1]], y=maxDevStateCanberra[gates[2]],  col = "brown3", pch = 24, cex=1.2,bg="white")
  arrows(referenceState[gates[1]],referenceState[gates[2]],data[,gates[1]][trefId+1],data[,gates[2]][trefId+1],length = 0.08, angle = 8,lwd=1)
  # mark reference space defining
  draw.circle(referenceState[gates[1]], y=referenceState[gates[2]], radius=ref_radius, border="grey", lty=2)
}

# Function for computing resilience
computeRL <- function(dx0xt, dx0x1) {
  return(2.0*dx0x1/(dx0x1+dx0xt)-1.0)
}

# restricting data set to one experiment
data <- data[data[1]>=experimentStart & data[1]<=experimentEnd,]

# compute deviation from reference state
data$euclidean <- apply(data[2:(numberOfGates+1)], 1, function(x) dist(rbind(referenceState, x), method = "euclidean"))
data$euclidean <- data$euclidean / sqrt(2)
data$canberra <- apply(data[2:(numberOfGates+1)], 1, function(x) dist(rbind(referenceState, x), method = "canberra"))
data$canberra <- data$canberra / numberOfGates
# compute radii for reference space
r_Euclidean <- max(data$euclidean[data$referencePhase==TRUE])
r_Canberra <- max(data$canberra[data$referencePhase==TRUE]) 

if (overrideMaxEuclidean == -1) {
  maxEuclideanId <- order(data$euclidean, decreasing=TRUE)[1]
} else {
  maxEuclideanId <- match(overrideMaxEuclidean, data[[1]])
}
if (overrideMaxCanberra == -1) {
  maxCanberraId <- order(data$canberra, decreasing=TRUE)[1]
} else {
  maxCanberraId <- match(overrideMaxCanberra, data[[1]])
}

maxDevStateEuclidean <- data[maxEuclideanId, 2:(numberOfGates+1)]
maxDevStateCanberra <- data[maxCanberraId, 2:(numberOfGates+1)]

maxDevStateEuclideanTime <- data[maxEuclideanId, 1]
maxDevStateCanberraTime <- data[maxCanberraId, 1]

# compute online version of resilience
data$maxCanberraOnline[1] <- data$canberra[1]
for (i in 2:nrow(data)) {
  if (data$canberra[i]< data$maxCanberraOnline[i-1]) {
    data$maxCanberraOnline[i] <- data$maxCanberraOnline[i-1]
  }else {
    data$maxCanberraOnline[i] <- data$canberra[i]
  }
}

# compute Resilience RL
data$RLeuclidean <- apply(data, 1, function(x) computeRL(x["euclidean"], data[maxEuclideanId, "euclidean"]))
data$RLcanberra <- apply(data, 1, function(x) computeRL(x["canberra"], data[maxCanberraId, "canberra"]))
data$RLcanberraOnline <- apply(data, 1, function(x) computeRL(x["canberra"], x["maxCanberraOnline"]))

# do NMDS plot

trefId <- match(tref, data[[1]])

mds.out <- metaMDS(rbind(data[,2:(numberOfGates+1)], referenceState), distance="bray", autotransform=FALSE, zerodist="add")
plot(mds.out, type="n", main=expression('Community evolve apart s'[ref]), font.main=1, cex.main = titleSize)
lines(mds.out$points[c((trefId+1):(length(mds.out$points[,1])-1)),], col = "black")
# marking points in reference space
if (trefId > 1) {
  points(mds.out$points[c(1:trefId),], col = "grey", pch = 21, bg="white")
} else {
  points(mds.out$points[c(trefId, trefId),], col = "grey", pch = 21, bg="white")
}
# mark reference state
points(mds.out$points[c(length(mds.out$points[,1]), length(mds.out$points[,1])),], col ="red", pch = 21, cex = 1.3, bg = "white")
# mark evolution
points(mds.out$points[c((trefId+1):(length(mds.out$points[,1])-1)),], bg = "white", pch = 21)
# mark direction
arrows(mds.out$points[length(mds.out$points[,1]),1],mds.out$points[length(mds.out$points[,1]),2],mds.out$points[(trefId+1),1],mds.out$points[(trefId+1),2], length = 0.08, angle = 8, lwd=1)
# mark most distant states
points(mds.out$points[c(maxEuclideanId, maxEuclideanId),], col = "deepskyblue3", pch = 24, cex = 1.2, bg = "white")
points(mds.out$points[c(maxCanberraId, maxCanberraId),], col = "brown3", pch = 24, cex = 1.2, bg = "white")
# mark final state
points(mds.out$points[c(length(mds.out$points[,1])-1, length(mds.out$points[,1])-1),], col = "black", pch = 16, cex=1.3)

# Plot evolution for dominant gates

labelStart <- "Relative abundance in"

plot2Dprojection(topGates, expression('Evolution in dominant gates of s'[ref]), paste(labelStart, topGates[1]), paste(labelStart, topGates[2]), r_Euclidean * sqrt(2))
plot2Dprojection(finalTopGates, expression('Evolution in dominant gates of s'[end]), paste(labelStart, finalTopGates[1]), paste(labelStart, finalTopGates[2]), r_Euclidean * sqrt(2))

# computing stability properties
RSeuclidean <- 1.0 - data[maxEuclideanId, "euclidean"]
RScanberra <- 1.0 - data[maxCanberraId, "canberra"]
DisSpeedeuclidean <- data[maxEuclideanId, "euclidean"] / (data[maxEuclideanId, 1] - tref)
DisSpeedcanberra <- data[maxCanberraId, "canberra"] / (data[maxCanberraId, 1] - tref)
ElasticityEuclidean <- (data[maxEuclideanId, "euclidean"] - data[nrow(data), "euclidean"]) / (data[nrow(data), 1] - data[maxEuclideanId, 1])
ElasticityCanberra <- (data[maxCanberraId, "canberra"] - data[nrow(data), "canberra"]) / (data[nrow(data), 1] - data[maxCanberraId, 1])

# plots over time 
plotData <- data[data[1] >= tref,]
plotData[1,2:(numberOfGates+1)] <- referenceState
plotData$euclidean[1] <- 0
plotData$canberra[1] <- 0
plotData$RLcanberraOnline[1] <- 0

# Distance over time
matplot(plotData[[1]], cbind(plotData$euclidean, plotData$canberra), ylim=c(0,1), main = "Resistance and displacement speed\n computing", ylab=expression('Deviation (d) from s'[ref]), xlab="Time", lty=c(1,1), pch=c(21,24), col=c("deepskyblue3","brown3"), type="l", font.main=1, cex.main=titleSize)
abline(h=c(r_Euclidean,r_Canberra),col=c("deepskyblue3","brown3"),lty=c(2,2))

points(plotData$Time, y=plotData$euclidean,pch=16,col="deepskyblue3")
points(data[maxEuclideanId, 1], y=data[maxEuclideanId, "euclidean"], col="deepskyblue3",cex=1.2,pch=24,bg="white")
points(experimentEnd, y=plotData[nrow(plotData), "euclidean"], col="black",cex=1.3,pch=21,bg="black")

points(plotData$Time, y=plotData$canberra,pch=16,col="brown3")
points(data[maxCanberraId, 1], y=data[maxCanberraId, "canberra"], col="brown3",cex=1.2,pch=24,bg="white")
points(experimentEnd, y=plotData[nrow(plotData), "canberra"], col="black",cex=1.3,pch=21,bg="black")

points(plotData$Time[1],plotData$euclidean[1], col="red", pch=21,bg="white",cex=1.3)
points(plotData$Time[1],plotData$canberra[1], col="red", pch=21,bg="white",cex=1.3)
legend("topright", c(paste("Euclidean, RS=", toString(round(RSeuclidean, digits=4)), ", DS=", toString(round(DisSpeedeuclidean, digits=4)), sep=""), paste("Canberra, RS=", toString(round(RScanberra, digits=4)), ", DS=", toString(round(DisSpeedcanberra, digits=4)), sep="")), pch=c(19,19), col=c("deepskyblue3", "brown3"), cex=0.7)

# Resilience over time
plot(data[maxEuclideanId:nrow(data), 1], y=data$RLeuclidean[maxEuclideanId:nrow(data)], ylim=c(0,1), font.main=1,main = "Resilience and elasticity computing\n(based on Euclidean distance)", ylab="Resilience (RL)", xlab="Time", xlim=c(plotData$Time[1], experimentEnd), type="l", col="deepskyblue3", cex.main=titleSize)
abline(h=0,col="black", lty=2)
points(data[maxEuclideanId:nrow(data),c(colnames(data)[1], "RLeuclidean")],pch=19,col="deepskyblue3")
points(data[maxEuclideanId, 1], y=data[maxEuclideanId, "RLeuclidean"], col="deepskyblue3",cex=1.2,pch=24,bg="white")
points(experimentEnd, y=plotData[nrow(plotData), "RLeuclidean"], col="black",cex=1.3,pch=21,bg="black")
legend("topright",paste("RL=", toString(round(plotData$RLeuclidean[nrow(plotData)], digits=4)),", E=", toString(round(ElasticityEuclidean, digits=4))), cex=0.7)

plot(data[maxCanberraId:nrow(data), 1], y=data$RLcanberra[maxCanberraId:nrow(data)], col="brown3", ylim=c(0,1), font.main=1, main = "Resilience and elasticity computing\n(based on Canberra distance)", ylab="Resilience (RL)", xlab="Time", xlim=c(plotData$Time[1], experimentEnd), type="l", cex.main=titleSize)
abline(h=0,col="black", lty=2)
points(data[maxCanberraId:nrow(data),c(colnames(data)[1], "RLcanberra")], pch=19,col="brown3")
points(data[maxCanberraId, 1], y=data[maxCanberraId, "RLcanberra"], col="brown3",cex=1.2,pch=24,bg="white")
points(experimentEnd, y=plotData[nrow(plotData), "RLcanberra"], col="black",cex=1.3,pch=21,bg="black")
legend("topright", paste("RL=", toString(round(plotData$RLcanberra[nrow(plotData)], digits=4)),", E=", toString(round(ElasticityCanberra, digits=4))), cex=0.7)

# Online resilience
plot(plotData[[1]], y=plotData$RLcanberraOnline, ylim=c(0,1), main = "Online resilience computing\n(based on Canberra distance)", ylab="Resilience", xlab="Time", col="black", pch=21, bg="white", font.main=1, type="l", cex.main = titleSize)
abline(h=0,col="black", lty=2)
points(plotData[[1]], y=plotData$RLcanberraOnline, pch=21,bg="white", col="black")
points(plotData$Time[1],plotData$RLcanberraOnline[1], col="red", pch=21,bg="white",cex=1.3)

if (exists("plotGatesStart")) {
  rm("plotGatesStart")
}
if (exists("plotGatesEnd")) {
  rm("plotGatesEnd")
}

dev.off()

