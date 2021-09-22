#
# Script to generate random sample data describing a set of compounds
# docked against a target and a subset of top compounds by efficiency.
#
# Author: Natalia Nikitina (nikitina@krc.karelia.ru)
#

N <- 500 # Total number of compounds
NSelected <- 5  # Number of selected compounds

dissimilarities<-matrix(ncol=N+3, nrow=N); 
locations<-matrix(ncol=5, nrow=N);

for (i in 1:N) {
  dissimilarities[i,1]=sprintf("Compound%d", i)
  dissimilarities[i,2]=runif(1, min=-500, max=-300) 
  dissimilarities[i,3]=runif(1, min=-15, max=-12)
  for(j in 4:(N+3)) {
    if(i==j-3) {dissimilarities[i,j]=1} else {dissimilarities[i,j] <- runif(1, min=0, max=1); dissimilarities[j-3,i+3] <- dissimilarities[i,j]}
  }
  locations[i,1]=sprintf("Compound%d",i);
  locations[i,2]=dissimilarities[i,2]
  locations[i,3]=dissimilarities[i,3]
  locations[i,4]=runif(1, min=0, max=15);
  locations[i,5]=runif(1, min=0, max=1);
}

selectedCompounds=locations[,1][sort(locations[,2], decreasing=TRUE, index.return=TRUE)$ix[1:NSelected]]
write.table(dissimilarities, file = "sample_dissimilarities.dat", row.names=FALSE, col.names=FALSE, sep = " ", quote = FALSE)
write.table(locations, file = "sample_locations.dat", row.names=FALSE, col.names=FALSE, sep = " ", quote = FALSE)
write.table(selectedCompounds, file = "sample_selected.dat", row.names=FALSE, col.names=FALSE, sep = "\n", quote = FALSE)