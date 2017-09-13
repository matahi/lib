#### Libraries and data
library(hdp) ### This is the package that Nicola Roberts developed to extract mutation signatures - you can find it here https://github.com/nicolaroberts/hdp/tree/master
#### Careful select the old version

#### Preparation of data
set.seed(42)
binomialImpute <- function(x) {x[is.na(x)] <- rbinom(sum(is.na(x)), 1, mean(x, na.rm=TRUE)); return(x)}
d <- as.matrix(mut.df) 
d[!d %in% c(0,1)] <- NA
genotypesImputed <- apply(d,2, binomialImpute)
# genotypesImputed <- round(ImputeMissing(d))
# Replace by zero


# clin.features <- c("hyperdiploidy", "t.11.14.", "t.14.16.", "t.14.20.", "t.4.14", "t.6.14.", "t.8.14.")
# clin.data <- genotypesImputed[,which(colnames(genotypesImputed) %in% clin.features)]
# 
# genotypesImputed <- genotypesImputed[, -which(colnames(genotypesImputed) %in% clin.features)]

#### Setting up hdp
#DPsetup
n <- ncol(genotypesImputed)
shape <- 1
invscale <- 1
hdp <- hdp_init(ppindex=0, #index of the parent DP for initial DP
cpindex=1, #index of alphaa and alphab for initial DP
hh=rep(1/n,n), #params for base distn (uniform Dirichlet)
alphaa=shape,
alphab=invscale)

hdp <- hdp_adddp(hdp,
numdp=nrow(genotypesImputed), # one DP for every sample in that cancer type
pp=1, # parent DP for group i is the i-th+1 overall DP because of the grandparent at position 1
cp=1) # index of alphaa and alphab for each DP

# Assign the data from each patient to a child DP
hdp <- hdp_setdata(hdp = hdp, dpindex=1:nrow(genotypesImputed)+1, data=genotypesImputed)

# Activate the DPs with specified number of classes (signatures)
# hdp <- dp_activate(hdp, 1:(nrow(genotypesImputed)+1), 5)
hdp <- dp_activate(hdp, 1:(nrow(genotypesImputed)+1), 7)
# hdp <- dp_activate(hdp, 1:(nrow(genotypesImputed)+1), 3) ## mm tries

#DP parameters
# DPparam
burnin <- 5000
#postsamples <- 1000
postsamples <- 10000
spacebw <- 20
cpsamples <- 10
# cpsamples <- 5 ## mm tries

#DPrun, cache=TRUE, cache.lazy=FALSE
set.seed(42)
output <- hdp_posterior(hdp, #activated hdp structure
numburnin=burnin,
numsample=postsamples,
numspace=spacebw,
doconparam=cpsamples)

# DP checks

save(output, file="results/clustering_update_full/output_hdp_posterior.RData")

plot(output$lik, type='l'); abline(v=burnin, lty=3)

pdf('results/clustering_update_full/classes.pdf')
plot(output$numclass, type='l')
dev.off()

library(ggplot2)
# class.info <- data.frame(num_class= output$numclass, iteration= 1:length(output$numclass))

# pdf('results/clustering_update_full/hist_classes.pdf')
# ggplot(class.info) + geom_histogram(aes(x=num_class)) + theme_bw()
# dev.off()

output <- get(load("results/clustering_update_full/output_hdp_posterior.RData"))

# AML classes
#posteriorMerged, cache=TRUE
#posteriorMerged <- hdp_extract_signatures(output, prop.explained=0.95, cos.merge=0.99)
posteriorMerged <- hdp_extract_signatures(output, prop.explained=0.99, cos.merge=0.99)
### varying prop.explained change the number of clusters

## plot number of class as a function of prop.explained

## prop.vect <- seq(0.9,1,0.01)
## 
## clust.size <- sapply(1:length(prop.vect), function(x)
##                      {
##                              print(x)
##                              posteriorMerged <- hdp_extract_signatures(output,prop.explained=prop.vect[x],cos.merge=0.99)
##                              return(ncol(posteriorMerged$sigs_qq[[1]]))
##                      })
## 
## save(clust.size, file="../../big_data/clust_size_hdp.RData")
## 
## pdf('../../results/HDP/prop_size_09_1.pdf')
## plot(prop.vect, clust.size)
## dev.off()

#classes, results='asis', cache=TRUE
#posteriorMeans <- Reduce("+",posteriorMerged$sigs_qq)/length(posteriorMerged$sigs_qq)
posteriorSamples <- array(unlist(posteriorMerged$sigs_qq), dim=c(dim(posteriorMerged$sigs_qq[[1]]), length(posteriorMerged$sigs_qq)))
rownames(posteriorSamples) <- colnames(genotypesImputed)
colnames(posteriorSamples) <- 1:ncol(posteriorSamples) -1
posteriorMeans <- rowMeans(posteriorSamples, dim=2)
posteriorQuantiles <- apply(posteriorSamples, 1:2, quantile, c(0.025,.5,0.975))
posteriorMode <- apply(posteriorSamples, 1:2, function(x) {t <- table(x); as.numeric(names(t)[which.max(t)])})

library(knitr)
#kable(posteriorQuantiles[2,,], "html", table.attr = 'id="posteriorMedian"') # Posterior median
# kable(posteriorQuantiles[2,,], "html", table.attr = 'id="posteriorMedian"') # Posterior median

# the five most prevalent lesions contributing to class assignment - 
genes <- apply(posteriorMeans, 2, function(x) {paste(ifelse(x>10,rownames(posteriorMeans),"")[order(x, decreasing = TRUE)[1:5]], collapse=";")})
genes <- gsub(";+$","",genes)
genes

write(genes, file="results/clustering_update_full/genes.txt")

#Assignment from posterior samples

library(RColorBrewer)
col <- c(brewer.pal(9,"Set1")[c(9,1:8)], brewer.pal(8,"Dark2"))
posteriorProbability <- apply(sapply(posteriorMerged$sigs_nd_by_dp, colMeans)[,-1],2,function(x) (x+.Machine$double.eps)/sum(x+.Machine$double.eps))
o <- order(apply(posteriorProbability,2,which.max))

pdf('results/clustering_update_full/probability_assign.pdf')
barplot(posteriorProbability[,o], col=col, border=NA, ylab="Probability", xlab="Patient")
#barplot(posteriorProbability[,o[10:100]], col=col, border=NA, ylab="Probability", xlab="Patient")
dev.off()

Sum.info <- data.frame(Size= rowMeans(posteriorProbability)*ncol(posteriorProbability), Prob=rowMeans(posteriorProbability), genes)

### Classes
#DPclass
dpClass <- factor(apply(posteriorProbability, 2, which.max)-1)
table(dpClass)

pdf('results/clustering_update_full/ordered_probability_assign.pdf')
plot(seq(0,1,l=ncol(posteriorProbability)),sort(apply(posteriorProbability,2,max)), type='l', ylim=c(0,1) , xlab="Fraction of patients", ylab="Assignment probability")
dev.off()

pdf('results/clustering_update_full/distribution_assign.pdf')
boxplot(apply(posteriorProbability,2,max) ~ dpClass, col=col, ylab="Probability", xlab="Class")
dev.off()

rotatedLabel <- function(x0 = seq_along(labels), y0 = rep(par("usr")[3], length(labels)), labels, pos = 1, cex=1, srt=45, ...) {
w <- strwidth(labels, units="user", cex=cex)
h <- strheight(labels, units="user",cex=cex)
u <- par('usr')
p <- par('plt')
f <- par("fin")
xpd <- par("xpd")
par(xpd=NA)
text(x=x0 + ifelse(pos==1, -1,1) * w/2*cos(srt/360*2*base::pi), y = y0 + ifelse(pos==1, -1,1) * w/2 *sin(srt/360*2*base::pi) * (u[4]-u[3])/(u[2]-u[1]) / (p[4]-p[3]) * (p[2]-p[1])* f[1]/f[2] , labels, las=2, cex=cex, pos=pos, adj=1, srt=srt,...)
par(xpd=xpd)
}

pdf('results/clustering_update_full/gene_signature.pdf')
par(mar=c(6,3,1,1)+.1, cex=.8)
o <- order(colSums(genotypesImputed), decreasing=TRUE)
driverPrevalence <- t(sapply(split(as.data.frame(as.matrix(genotypesImputed)), dpClass), colSums)[o,])
b <- barplot(driverPrevalence, col=col, las=2, legend=TRUE, border=NA, args.legend=list(border=NA), names.arg=rep("", ncol(genotypesImputed)))
abline(h=seq(100,500,100), col="white")
rotatedLabel(b, labels=colnames(genotypesImputed)[o])

dev.off()

#Driver signatures

par(mar=c(6,3,1,1)+.1, cex=.8)
t <- table(dpClass)
i <- 0; for(c in levels(dpClass)){i <- 1+i 
pdf(paste0('results/clustering_update_full/class_signature_',i,'.pdf'),width=16, height=9)
b <- barplot(posteriorQuantiles[2,o,c]/t[i], col=col[i], las=2, legend=FALSE, border=NA,  names.arg=rep("", ncol(genotypesImputed)), main=paste("Class", i-1,t[i], "Patients","f =",round(t[i]/1540,2)), ylim=c(0,1))
segments(b, posteriorQuantiles[1,o,c]/t[i], b, posteriorQuantiles[2,o,c]/t[i], col="white")
segments(b, posteriorQuantiles[2,o,c]/t[i], b, posteriorQuantiles[3,o,c]/t[i], col=col[i])
rotatedLabel(b, labels=colnames(genotypesImputed)[o])
dev.off()
}


