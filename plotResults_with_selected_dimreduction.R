#
# Script to plot 2D diagram of dissimilarities in a set of ligands using datafile.
# Dissimilarity is defined as (1 - Tanimoto coefficient).
# 
#
# Datafile format of similarities within a set of N compounds C1, ..., CN:
#
# <Name of C1> <Ligand efficiency of C1> <Predicted binding energy of C1> <Similarity with C1> ... <Similarity with CN>
# ...
# <Name of CN> <Ligand efficiency of CN> <Predicted binding energy of CN> <Similarity with C1> ... <Similarity with CN>
#
#
# Authors: Steffen Moeller (steffen_moeller@gmx.de), Natalia Nikitina (nikitina@krc.karelia.ru)
#

library("plotfunctions")

args              <- commandArgs(TRUE)
datafileDistances <- "sample_similarities.dat" #args[1] 
datafileSelected  <- "sample_selected.dat" #args[2] 
datasetname       <- "Sample Dataset" #args[3] 
diagramfname      <- "sample_dimreduction.png" #args[4] 
byEnergy          <- FALSE   # Color plot points by predicted binding energy or by ligand efficiency

data.similarities <- read.table(datafileDistances,row.names=1)
colnames(data.similarities) <- c("Efficiency","Energy",rownames(data.similarities))
similarities <- data.similarities[-c(1,2)]
distances    <- 1-similarities
efficiency   <- data.similarities[,1]
energy       <- data.similarities[,2]

# Classical (metric) multidimensional scaling
fit <- cmdscale(distances, eig=T, k=2)

# Assign point colors depending on the values
colorLow <- '#006400'; colorMedium <- '#ADFF2F'; colorHigh <-'#FFA500'; # Your palette
sortedValues <- if(byEnergy) sort(energy) else sort(efficiency)
ns <- length(sortedValues); qm <- median(sortedValues); nm <- max(which(qm >= sortedValues))[[1]];
q1 <- quantile(sortedValues, 1/3); q2 <- quantile(sortedValues, 2/3); n1 <- max(which(q1 >= sortedValues))[[1]]; n2 <- max(which(q2 >= sortedValues))[[1]]
rbPal  <- colorRampPalette(c(colorLow, colorMedium, colorHigh)); Color <- rbPal(100)[as.numeric(cut(sortedValues, breaks = 100))]; 
LColor <- c(Color[1], Color[nm], Color[ns])
LData  <- c(sprintf("%.1f ... %.1f kcal/mol", sortedValues[1], sortedValues[n1]),
            sprintf("%.1f ... %.1f kcal/mol", sortedValues[n1]+0.1, sortedValues[n2]),
            sprintf("%.1f ... %.1f kcal/mol", sortedValues[n2]+0.1, sortedValues[ns]))

# Jitter the data of ligand similarities so that the points will not overlap
jittereddata <- jitter(fit$points, factor=20)

# Create a data frame of ZINC IDs, cluster numbers and energies
# and fill it with the set of selected ligands
# Join data frames and sort by cluster number, then by energy, then by efficiency
listSelectedCompounds <- scan(datafileSelected, what="", sep="\n")
NSelected <- length(listSelectedCompounds)
selectedCompounds <- matrix(ncol=5, nrow=NSelected)
for (i in 1:NSelected) {
  compoundName <- listSelectedCompounds[i]
  selectedCompounds[i,1] <- compoundName
  selectedCompounds[i,2] <- data.similarities[compoundName,"Energy"]
  selectedCompounds[i,3] <- data.similarities[compoundName,"Efficiency"]
  selectedCompounds[i,4] <- jittereddata[grep(listSelectedCompounds[i], rownames(jittereddata)),1][[1]]
  selectedCompounds[i,5] <- jittereddata[grep(listSelectedCompounds[i], rownames(jittereddata)),2][[1]]}
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
xLimits=range(jittereddata[,1])
yLimits=range(jittereddata[,2])
plot(jittereddata, xlab="", ylab="", xlim=xLimits, ylim=yLimits, col=Color, pch=20, cex=1) 
if(byEnergy) {legendTitle="Predicted binding energy"} else {legendTitle<-"Predicted ligand efficiency"}

# Plot the points of selected ligands from the clusters
points(x=as.numeric(as.character(selectedCompounds$X4)), y=as.numeric(as.character(selectedCompounds$X5)), pch='*', cex=2.5, col=SColor)
text(x=as.numeric(as.character(selectedCompounds$X4)), y=as.numeric(as.character(selectedCompounds$X5)), 
     labels=SLabels, cex=1, pos=sample(c(1:4), replace=T, size=nrow(distances)), col=SColor) 

gradientLegendTitle <- if(byEnergy) "Predicted binding energy (kcal/mol)" else "Ligand efficiency (kcal/mol)  "
text(xLimits[1]+(xLimits[2]-xLimits[1])/5, yLimits[2], gradientLegendTitle)
gradientLegend(valRange=c(min(sortedValues),max(sortedValues)), color=LColor, pos=c(0.05,0.90,0.35,0.93), pos.num=1, dec=3, side=1, inside=TRUE)
pos <- legend('bottomleft', legend=LSelected, col="black", adj=c(0.1,0.1), cex=0.9, pt.cex = 1, title="Selected ligands (*)")

dev.off()
