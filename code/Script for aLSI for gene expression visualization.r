########################################################################################################################
####
#### R-Script to reproduce the results (Figures and Tables) of the paper
#### "Asymmetric latent semantic indexing for gene expression experiments visualization"
#### Submitted to Advances in data analysis and classifciation.
####
#### Authors: Javier Gonzalez, Alberto Munoz and Gabriel Martos
####
#### 15/02/2014
####
########################################################################################################################
# requiered libraries 
library(geneplotter)
library(RColorBrewer) 
library(genefilter)
library(MASS)
library(mclust)
library(gplots)
library(mclust)
library(xtable)

### Load the data
data      <- read.table("nci.data")
labels    <- read.table("labels.txt")[,1]

# Expression estimates on the UN-LOGGED scale
e.mat     <- exp(data)

# look at mean, sd, & cv for each gene across arrays
gene.mean <- apply(e.mat,1,mean) 
gene.sd   <- apply(e.mat,1,sd) 
gene.cv   <- gene.sd/gene.mean

# Plot of the gene.sd over the gene mean and histogram of those with cv. larger that 7.7
blues.ramp <- colorRampPalette(brewer.pal(9,"Blues")[-1]) 
dCol       <- densCols(log(gene.mean),log(gene.sd), colramp=blues.ramp) 

###########
# Figure 3: cut-off for differential expression analysis 
###########
par(mfrow=c(1,2)) 
plot(gene.mean,gene.sd,log='xy',col=dCol,pch=16,cex=0.1,xlim = c(1,2),ylim = c(0.2,5),xlab="Gene mean", ylab="Gene sd.",main="A",cex.main=2) 
abline(v=100,lwd=3,col='red') 
hist(log(gene.cv),xlab = "Gene CV (log. scale)",col="grey",main="B",cex.main=2) 
abline(v=log(.5),lwd=3,col='red')
text(0.8,800,"Theshold: CV>0.5")

# We select those genes with CV larger.5 (genes differentially expressed in the data)
ffun          <- filterfun(cv(0.5,100))
t.fil         <- genefilter(e.mat,ffun)

# Apply filter and transform back on log scale
small.eset    <- e.mat[t.fil,]
small.noeset  <- e.mat[!t.fil,]

# Threshod to consider that a gene is expressed (the maximum value of the not expressed).
thres = max(small.noeset)

###########
# Figure 2: diff-expressed vs. not diff-expressed gene.
###########
par(mfrow=c(1,2))
par(mar=c(5,5,3,3))
plot(1:64,small.eset[99,],ylim=c(0,15),pch =20,xlab="Experiment",ylab="Expression",main="Differentially expressed gene",cex=2)
abline(h = thres ,lwd=2, col = "red",lty=2)
plot(1:64,small.noeset[2,],ylim=c(0,15),pch =20,xlab="Experiment",ylab="Expression",main="Non differentially expressed gene",cex=2)
abline(h = thres ,lwd=2, col = "red",lty=2)

# Select those genes expressed at least one time
ffun2         <- filterfun(pOverA(0.025,thres))
t.fil2        <- genefilter(small.eset,ffun2)
small.eset2   <- small.eset[t.fil2,]       

# Create genes by experiments matrix of differentially expressed genes
colnames(small.eset2) <- labels
X                     <- as.matrix((sign(small.eset2 - thres)+1)/2)  # genes by experiments matrix
colnames(X)           <- labels
n.genes               <- nrow(X)

###########
# Figure 2: Zipf's law, histogram of the norm of the genes
###########
par(mfrow=c(1,1))
norm.genes = apply(X,1,sum)
hist(norm.genes,nclass=25,col = "grey",main = " ",xlab="Genes Norm",ylab="Frequency")
text(16,600,paste(n.genes," differentially expressed genes"),cex=1.7)

###########
# Figure 1 (A and B): Heatmaps pf raw data and differentially expressed genes
###########
heatmapA <- heatmap(log(t(small.eset2)), Rowv=NA, Colv=NA,scale="column",margins = c(0.5, 15))
heatmapB <- heatmap(t(X), Rowv=NA, Colv=NA, scale="none",margins = c(0.5, 15))

# Define the degree of membership of each gene to each type of experiment
n.groups                          <- length(table(labels))
membership.genes.probs            <- matrix(0,n.genes,n.groups )
membership.genes                  <- matrix(0,n.genes,n.groups )
colnames(membership.genes.probs)  <- levels(labels)
colnames(membership.genes)        <- levels(labels)

for (i in 1:n.genes)
  {
  membership.genes.probs[i,] = table(labels[X[i,]==1])/norm.genes[i]
  membership.genes[i,which(membership.genes.probs[i,] == max(membership.genes.probs[i,]))]=1  
  }

### Calculate the (asymmetric) matrix of similarities between genes with respect to the experiments
S0 		            <- X %*% t(X) 
Matrix.norms.inv	<- diag(1/norm.genes) 
S 		            <- Matrix.norms.inv %*% S0 
 
## Calculate the polar decomposition
svdS              <- svd(S)
M1                <- svdS$u%*%diag(svdS$d)%*%t(svdS$u)  
M2                <- svdS$v%*%diag(svdS$d)%*%t(svdS$v) 

### Calculate the (asymmetric) matrix of similarities between genes with respect to the classes
Z                 <- membership.genes%*%t(membership.genes)%*%diag(1/apply(membership.genes,1,sum)) 
svdZ              <- svd(Z)
R1                <- svdZ$u%*%diag(svdZ$d)%*%t(svdZ$u)  
R2                <- svdZ$v%*%diag(svdZ$d)%*%t(svdZ$v) 

## Matrix combination of the two asymmetric similarities
A                 <- (M1+M2)/2 + 0.2*(R1+R2) 

## Extract components
eigA              <- eigen(A)
XX                <- Re(eigA$vectors[,Re(eigA$values)>0.0001])%*%diag(sqrt(Re(eigA$values)[Re(eigA$values)>0.0001]))

# Extract components of the mixture
cluster.lsa       <- Mclust(XX, G =14)

# Cross the clusters with the original classes
cross = matrix(NA,14,14)
for (k in 1:14){cross[k,] = table(membership.genes[,k],cluster.lsa$classification)[2,]}
rownames(cross)   <- levels(labels)
colnames(cross)   <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14") 

###########
# Table 2: clusters vs. orignal classes
###########
xtable(cross)

###########
# Figure 6: MDS of the groups
##########
X.cross = sammon(dist(cross))
plot(X.cross$points,xlab= "Component 1",ylab= "Component 2",type="n",col=1:14,xlim=c(-270,270),ylim=c(-270,270),main = "aLSI mapping of types of Cancer")
text(X.cross$points,rownames(X.cross$points),col="blue",cex=1.1)
                 

### Obtain genes with higest probability in each cluster
mp                <- matrix(NA,5,14)
for (i in 1:14){mp[,i] = order(cluster.lsa$z[,i],decreasing = TRUE)[1:5]} 
colnames(mp)      <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14") 
rownames(mp)      <- c("gene 1","gene 2","gene 3","gene 4","gene 5") 

###########
# Table 1: most probable genes of each cluster
###########
xtable(t(mp))


## Sammon mapping for the combination
mds.S12         <- prcomp(XX)$x[,c(1,2)]  ## save similarity matrix
mds.S13         <- prcomp(XX)$x[,c(1,3)]  ## save s 
cluster.lsa     <- hclust(dist(XX),method="ward")
  
### MDS-Correlation  
C               <- cor(log(t(small.eset2)))
mds.C12         <- cmdscale(1-C,k=3)[,c(1,2)]
mds.C13         <- cmdscale(1-C,k=3)[,c(1,3)]

### MDS-Euclidean 
D               <- dist(log(small.eset2)) 
mds.D12         <- cmdscale(D,k=3)[,c(1,2)]
mds.D13         <- cmdscale(D,k=3)[,c(1,3)]

###########
# Figure 5: genes maps produced by the aLSI, Pearson's correlation and Euclidean distance. 
###########
par(mfrow = c(3,2))
par(mar= c(5,5,1.5,1))
plot(mds.S12,pch=19,cex=0.8,col="grey",type="n",main = "Asymmetric LSI",xlab="1st component",ylab="2nd component",xlim=c(-1,1.4))
for(k in 1:14)
points(mds.S12[which(apply(membership.genes,1,sum)==1 & membership.genes[,k]==1), ],col=k,pch=k,cex=0.9)
legend("topleft",  c("G1","G2","G3","G4","G5","G6","G7"), col=1:7,pch=1:7) 
legend("topright", c("G8","G9","G10","G11","G12","G13","G14"), col=8:14,pch=8:14)  
 
plot(mds.S13,pch=19,cex=0.8,col="grey",type="n",main = "Asymmetric LSI",xlab="1st component",ylab="3th component",xlim=c(-1,1.4))
for(k in 1:14)
points(mds.S13[which(apply(membership.genes,1,sum)==1 & membership.genes[,k]==1), ],col=k,pch=k,cex=0.9,xlab="1st component")
legend("topleft",  c("G1","G2","G3","G4","G5","G6","G7"), col=1:7,pch=1:7) 
legend("topright", c("G8","G9","G10","G11","G12","G13","G14"), col=8:14,pch=8:14)  
  
plot(mds.C12,pch=19,cex=0.8,col="grey",main = "Correlation",xlab="1st component",ylab="2nd component",xlim=c(-1,1.2))  
for(k in 1:14)
points(mds.C12[which(apply(membership.genes,1,sum)==1 & membership.genes[,k]==1), ],col=k,pch=k,cex=0.9)
legend("topleft",  c("G1","G2","G3","G4","G5","G6","G7"), col=1:7,pch=1:7) 
legend("topright", c("G8","G9","G10","G11","G12","G13","G14"), col=8:14,pch=8:14)  

plot(mds.C13,pch=19,cex=0.8,col="grey",main = "Correlation",xlab="1st component",ylab="3th component",xlim=c(-1,1.2))  
for(k in 1:14)
points(mds.C13[which(apply(membership.genes,1,sum)==1 & membership.genes[,k]==1), ],col=k,pch=k,cex=0.9)
legend("topleft",  c("G1","G2","G3","G4","G5","G6","G7"), col=1:7,pch=1:7) 
legend("topright", c("G8","G9","G10","G11","G12","G13","G14"), col=8:14,pch=8:14)  

plot(mds.D12,pch=19,cex=0.8,col="grey",main = "Euclidean distance",xlab="1st component",ylab="2nd component",xlim=c(-20,20))
for(k in 1:14)
points(mds.D12[which(apply(membership.genes,1,sum)==1 & membership.genes[,k]==1), ],col=k,pch=k,cex=0.9)
legend("topleft",  c("G1","G2","G3","G4","G5","G6","G7"), col=1:7,pch=1:7) 
legend("topright", c("G8","G9","G10","G11","G12","G13","G14"), col=8:14,pch=8:14)  

plot(mds.D13,pch=19,cex=0.8,col="grey",main = "Euclidean distance",xlab="1st component",ylab="3th component",xlim=c(-20,20))
for(k in 1:14)
points(mds.D13[which(apply(membership.genes,1,sum)==1 & membership.genes[,k]==1), ],col=k,pch=k,cex=0.9)
legend("topleft",  c("G1","G2","G3","G4","G5","G6","G7"), col=1:7,pch=1:7) 
legend("topright", c("G8","G9","G10","G11","G12","G13","G14"), col=8:14,pch=8:14)  
                                                                                  
  
