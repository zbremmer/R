#
# This file includes random code snippets for various data science procedures. 
#

##################
# One hot encoding
library(dummies)
df <- dummy.data.frame(df, names=c("MyField"), sep="_") # Original MyField col will not be in new df

################################
# Simple heirarchical clustering
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
# Bootstrap clustering
library(fpc)
