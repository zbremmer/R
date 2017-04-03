#
# This file includes random code snippets for various data science procedures. 
#

##################
# One hot encoding
library(dummies)
df <- dummy.data.frame(df, names=c("MyField"), sep="_") # Original MyField col will not be in new df

################################
# Simple heirarchical clustering
# Source : Practical Data Science with R (book)
m <- matrix(nrow = 1, ncol = 1, colnames = c("item")) # Data (can also use df, just need to change labels below)
d <- dist(m, method = "euclidian")
pfit <- hclust(d, method = "ward.D")

plot(pfit, labels = m$item) # draws dendrogram (for df, use df[,"Item"])
rect.hclust(pfit, k=5) # k is number of clusters sugg by dendrogram

  # Gives a vector of the cluster # each row belongs to
groups <- cutree(pfit, k = 5) 

  # Plot first two principal components (can plot any 2, but top 2 most likely to show useful info)
library(ggplot2)
princ <- prcomp(m)
nComp <- 2 
project <- predict(princ, newdata = m)[,1:nComp] # Will rotate data 
project.plus <- cbind(as.data.frame(project), cluster = as.factor(groups), label = m$item) # for df, use df[,"Item"]
ggplot(project.plus, aes(x=PC1, y=PC2)) + geom_point(aes(shape=cluster)) + geom_text(aes(label = label), hjust=0, vjust=1)

############################
# Total Within Sum of Squares (WSS)
# Source : Practical Data Science with R (book)

sqr_edist <- function(x, y) {
  # Calculate squared distance between two vectors
  sum((x - y)^2)
}

wss_cluster <- function(clustermat) {
  # Calculate WSS for a single cluster
  c0 <- apply(clustmat, 2, FUN = mean)
  sum(apply(clustermat, 1, FUN = function(row){sqr_edist(row, c0)}))
}

wss_total <- function(dmatrix, labels){
  wsstot <- 0
  k <- length(unique(labels))
  for(i in 1:k){
    wsstot <- wsstot + wss.cluster(subset(dmatrix, labels==i))
  }
  wsstot
}

## Need iteration to test different values of k
## Should also take method as part of function so we can use kmeans, hclust, etc.
## kmeans returns WSS as one of its outputs (see ch_criterion() below)

############################
# Calinski-Harabasz Index
# Source : Practical Data Science with R (book)
# Works for hclust() and kmeans()

totss <- function(dmatrix){
  # Convenience function to calculate total sum of squares (TSS)
  grandmean <- apply(dmatrix, 2, FUN = mean)
  sum(apply(dmatrix, 1, FUN = function(row){sqr_edist(row, grandmean)}))
}

ch_criterion <- function(dmatrix, kmax, method="kmeans"){
  # Calculate CH index for a number of clusters from 1:kmax
  if(!(method %in% c("kmeans", "hclust"))) {
    stop("Method must be one of c('kmeans', 'hclust')")
  }
  
  npts <- dim(dmatrix)[1] # number of rows
  totss <- totss(dmatrix) # TSS independent of clustering
  wss <- numeric(kamx)
  crit <- numeric(kmax)
  wss[1] <- (npts - 1) * sum(apply(dmatrix, 2, var)) # Calculate WSS for k=1 (which is really just TSS)
  
  for(k in 2:kmax){
    if(method == "kmeans"){
      clustering <- kmeans(dmatrix, k, nstart = 10, iter.max = 100)
      wss[k] <- clustering$tot.withinss #kmeans returns WSS as one of its outputs
    } else {
      d <- dist(dmatrix, method = "euclidian")
      pfit <- hclust(d, method = "ward")
      labels <- cutree(pfit, k = k)
      wss[k] <- wss_total(dmatrix, labels)
    }
  }

  bss <- totss - wss
  crit_num <- bss/(0:(kmax - 1))
  crit_denom <- wss/(npts - 1:kmax)
  list(crit = crit_num/crit_denom, wss = wss, totss = totss)
}
