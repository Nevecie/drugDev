# 
# Script to plot 2D diagram of pre-calculated compound similarities with a reference ligand
# against pre-calculated distances from the reference ligand in its binding pocket
# with a subset of pre-selected compounds displayed in a different style.
# 
# Similarity between two compounds is defined as the Tanimoto coefficient of their 
# molecular fingerprints. Marker color reflects ligand efficiency or predicted binding energy.
#
# Datafile format of a set of N compounds C1, ..., CN:
#
# <Name of C1> <Ligand efficiency of C1> <Predicted binding energy of C1> <X1> <Y1>
# ...
# <Name of CN> <Ligand efficiency of CN> <Predicted binding energy of CN> <XN> <YN> 
#
#
# Authors: Steffen Moeller (steffen_moeller@gmx.de), Natalia Nikitina (nikitina@krc.karelia.ru)
#

library("plotfunctions")

args              <- commandArgs(TRUE)
datafileLocations <- "sample_locations.dat" #args[1]
datafileSelected  <- "sample_selected.dat" #args[2] 
datasetname       <- "Sample Dataset" #args[3] 
diagramfname      <- "sample_locations.png" #args[4]
byEnergy          <- FALSE  # Color plot markers by predicted binding energy 

data.locations <- read.table(datafileLocations, row.names=1)
colnames(data.locations) <- c("Efficiency", "Energy", "X", "Y")
efficiency     <- data.locations[,1]
energy         <- data.locations[,2]
pointCoords     <- matrix(ncol=2, nrow=nrow(data.locations))
pointCoords[,1] <- data.locations[,3]
pointCoords[,2] <- data.locations[,4]
rownames(pointCoords) <- rownames(data.locations)

# Assign point colors depending on the values
colorLow <- '#006400'; colorMedium <- '#ADFF2F'; colorHigh <-'#FFA500'; # Your palette
sortedValues <- if(byEnergy) sort(energy) else sort(efficiency)
ns <- length(sortedValues); qm <- median(sortedValues); nm <- max(which(qm >= sortedValues))[[1]];
q1 <- quantile(sortedValues,1/3); q2 <- quantile(sortedValues,2/3); n1 <- max(which(q1 >= sortedValues))[[1]]; n2 <- max(which(q2 >= sortedValues))[[1]]
rbPal  <- colorRampPalette(c(colorLow, colorMedium, colorHigh)); Color  <- rbPal(100)[as.numeric(cut(sortedValues, breaks = 100))]; 
LColor <- c(Color[1], Color[nm], Color[ns])
LData  <- c(sprintf("%.1f ... %.1f kcal/mol", sortedValues[1], sortedValues[n1]),
            sprintf("%.1f ... %.1f kcal/mol", sortedValues[n1]+0.1, sortedValues[n2]),
            sprintf("%.1f ... %.1f kcal/mol", sortedValues[n2]+0.1, sortedValues[ns]))

listSelectedCompounds <- scan(datafileSelected, what="", sep="\n")
NSelected <- length(listSelectedCompounds)
selectedCompounds <- matrix(ncol=5, nrow=NSelected)
for (i in 1:NSelected) {
  compoundName <- listSelectedCompounds[i]
  selectedCompounds[i,1] <- compoundName
  selectedCompounds[i,2] <- data.locations[compoundName,"Energy"]
  selectedCompounds[i,3] <- data.locations[compoundName,"Efficiency"]
  selectedCompounds[i,4] <- pointCoords[grep(listSelectedCompounds[i], rownames(pointCoords)),1][[1]]
  selectedCompounds[i,5] <- pointCoords[grep(listSelectedCompounds[i], rownames(pointCoords)),2][[1]]}
selectedCompounds <- data.frame(selectedCompounds)
selectedCompounds <- selectedCompounds[order(as.numeric(as.character(selectedCompounds$X2)), as.numeric(as.character(selectedCompounds$X3))),]

SColor <- c(Color[as.numeric(rownames(selectedCompounds[1,]))])
SLabels <- vector(mode = "list", length = NSelected)
LSelected <- NULL
for(row in 1:NSelected) {SColor[row]=Color[as.numeric(rownames(selectedCompounds[row,]))]; 
if(row<10) {SLabels[row]=paste(' (',row,')',sep='')} # intToUtf8(96+row) for letters
else {SLabels[row]=paste('(',row,')',sep='')}
LSelected <- c(LSelected,sprintf("%s %s",SLabels[row],selectedCompounds$X1[row]))}

# Plot the data of ligand similarities
graphics.off()
par("mar")
par(mar=c(1,1,1,1))
png(diagramfname, width = 8, height = 8, units = "cm", res=2400, pointsize = 6)
xLimits=range(pointCoords[,1])
yLimits=c(0,1) #range(pointCoords[,2])
plot(pointCoords, xlab="", ylab="", xlim=xLimits, ylim=yLimits, yaxs="i", col=Color, pch=20, cex=1) 
if(byEnergy) {legendTitle="Predicted binding energy"} else {legendTitle<-"Predicted ligand efficiency"}

# Plot the data
points(x=as.numeric(as.character(selectedCompounds$X4)), y=as.numeric(as.character(selectedCompounds$X5)), pch='*', cex=2.5, col=SColor)
text(x=as.numeric(as.character(selectedCompounds$X4)), y=as.numeric(as.character(selectedCompounds$X5)), 
     labels=SLabels, cex=1, pos=sample(c(1:4), replace=T, size=nrow(pointCoords)), col=SColor) 

gradientLegendTitle <- if(byEnergy) "Predicted binding energy (kcal/mol)" else "Ligand efficiency (kcal/mol)  "
text(xLimits[1]+(xLimits[2]-xLimits[1])/5, yLimits[2]-(yLimits[2]-yLimits[1])/20, gradientLegendTitle)
gradientLegend(valRange=c(min(sortedValues),max(sortedValues)), color=LColor, pos=c(0.05,0.90,0.35,0.93), pos.num=1, dec=3, side=1, inside=TRUE)
pos <- legend('bottomleft', legend=LSelected, col="black", adj=c(0.1,0.1), cex=0.9, pt.cex = 1, title="Selected ligands (*)")

dev.off()
