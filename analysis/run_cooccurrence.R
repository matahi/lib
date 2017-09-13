run_cooccurrence <- function(genomicData, out.dir="", cutoff.occurrence=10, col.axis=NULL)
{
        # genomicData <- out.merge$dat
        # col.axis <- colors.feature_type[grouping_Cardiff]

        genomicData <- genomicData[, colSums(genomicData,na.rm=T) >= cutoff.occurrence]

        genomicGroups <- factor(colnames(genomicData),levels=colnames(genomicData))

        ##
        logPInt <- sapply(1:ncol(genomicData), function(i) sapply(1:ncol(genomicData), function(j) {f<- try(fisher.test(genomicData[,i], genomicData[,j]), silent=TRUE); if(class(f)=="try-error") 0 else ifelse(f$estimate>1, -log10(f$p.val),log10(f$p.val))} ))
        odds <- sapply(1:ncol(genomicData), function(i) sapply(1:ncol(genomicData), function(j) {f<- try(fisher.test(table(genomicData[,i], genomicData[,j])), silent=TRUE); if(class(f)=="try-error") f=NA else f$estimate} ))
        pairs <- sapply(1:ncol(genomicData), function(i) colMeans(genomicData * genomicData[,i], na.rm=TRUE))
        diag(logPInt) <- 0
        diag(odds) <- 1
        colnames(odds) <- rownames(odds) <- colnames(logPInt) <- rownames(logPInt) <- colnames(genomicData)
        odds[odds<1e-3] = 1e-4
        odds[odds>1e3] = 1e4
        odds[10^-abs(logPInt) > 0.05] = 1
        logOdds=log10(odds)
        diag(logPInt) <- NA

        ########### HEATMAP ALL 
        library(RColorBrewer)
        library(gplots)

        pdf(paste0(out.dir,"/cooccurrence.pdf"), width=16, height=10, pointsize=8) 
        par(bty="n", mgp = c(2,.5,0), mar=c(4,4,4,4)+.1, las=2, tcl=-.33)

        ix = TRUE#colnames(interactions) %in% colnames(all_genotypes)
        d <- dist(t(genomicData[,ix]) + 10*as.numeric(genomicGroups))
        h = hclust(d, method="average")
        o = order(genomicGroups,-colSums(genomicData, na.rm=TRUE))#order(cmdscale(d, k=1))#h$order #c(h$order,(length(h$order) +1):ncol(interactions))
        M <-  matrix( NA, ncol=ncol(odds), nrow=nrow(odds))
        M[lower.tri(M)] <- cut(logOdds[o,o][lower.tri(M)], breaks = c(-4:0-.Machine$double.eps,0:4), include.lowest=TRUE)
        M[upper.tri(M, diag=TRUE)] <- as.numeric(cut(pairs[o,o][upper.tri(M, diag=TRUE)]*nrow(genomicData), breaks=c(-1,0,5,10,20,50,100,200,600))) + 9 
        ## Need to print max(M) then adjust colours accordingly
        # print(max(M))
        image(x=1:ncol(logPInt), y=1:nrow(logPInt), M, col=c(brewer.pal(9,"BrBG"), c("white",brewer.pal(7,"Blues"))), breaks=0:max(M,na.rm=TRUE), xaxt="n", yaxt="n", xlab="",ylab="", xlim=c(0, ncol(logPInt)+3), ylim=c(0, ncol(logPInt)+3))
        #image(x=1:ncol(logPInt), y=1:nrow(logPInt), M, col=c(brewer.pal(8,"BrBG"), c("white",brewer.pal(7,"Blues"))), breaks=0:max(M,na.rm=TRUE), xaxt="n", yaxt="n", xlab="",ylab="", xlim=c(0, ncol(logPInt)+3), ylim=c(0, ncol(logPInt)+3))
        l <- colnames(logPInt)[o]

        mtext(side=1, at=1:ncol(logPInt), l, cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1))
        mtext(side=2, at=1:ncol(logPInt), l, cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1))

        # colors of axes 
        # if (is.null(col.axis))
        # {
        #         mtext(side=1, at=1:ncol(logPInt), l, cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1))
        #         mtext(side=2, at=1:ncol(logPInt), l, cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1))
        # } else
        # {
        #         iy <- 1:ncol(logPInt)
        #         for (i in 1:length(unique(col.axis)))
        #         {
        #                 index.i <- which(col.axis==unique(col.axis)[i])
        #                 axis(1, iy[index.i], labels = l[index.i], las = 2, line = -0.5, tick = 0,
        #                      col.axis=unique(col.axis)[i], # added line
        #                      cex.axis = .66)

        #                 axis(2, iy[index.i], labels = l[index.i], las = 2, line = -0.5, tick = 0,
        #                      col.axis=unique(col.axis)[i], # added line
        #                      cex.axis = .66)
        #         }
        # }

        # beautify
        abline(h=0:ncol(logPInt)+.5, col="white", lwd=.5)
        abline(v=0:ncol(logPInt)+.5, col="white", lwd=.5)

        ## rectangle around features
        if (!is.null(col.axis))
        {
                changing_name <- sapply(2:length(col.axis), function(n){col.axis[n]!=col.axis[n-1]})
                pos <- c(0.5, which(changing_name)+0.5,length(col.axis)+0.5)
                rect( xleft=pos[-length(pos)], ybottom=pos[-length(pos)], xright=pos[-1], ytop=pos[-1], border=unique(col.axis), lwd=2  )
        }

        #
        P <- 10^-abs(logPInt[o,o])
        P[upper.tri(P)] <- NA
        w = arrayInd(which(p.adjust(P, method="BH") < .1), rep(nrow(logPInt),2))
        points(w, pch=".", col="salmon3")
        w = arrayInd(which(p.adjust(P) < .05), rep(nrow(logPInt),2))
        points(w, pch="*", col="salmon3")

        # Frequency position
        # image(y = 1:9 +18, x=rep(ncol(logPInt),2)+c(2,3), z=matrix(c(1:8), nrow=1), col=c("white",brewer.pal(7,"Blues")), add=TRUE)
        # axis(side = 4, at = seq(1,7) + 19, cex.axis=.66, tcl=-.15, label=c(1,5,10,20,50,100,200), las=1, lwd=.5)
        # mtext(side=4, at=28, "Frequency", las=2, line=-1,cex=.66)

        image(y = 1:9 +12, x=rep(ncol(logPInt),2)+c(2,3), z=matrix(c(1:8), nrow=1), col=c("white",brewer.pal(7,"Blues")), add=TRUE)
        axis(side = 4, at = seq(1,7) + 13, cex.axis=.66, tcl=-.15, label=c(1,5,10,20,50,100,200), las=1, lwd=.5)
        mtext(side=4, at=22, "Frequency", las=2, line=-1,cex=.66)

        # Odd-ratio position
        # image(y = 1:8 +5, x=rep(ncol(logPInt),2)+c(2,3), z=matrix(c(1:8), nrow=1), col=brewer.pal(8,"BrBG"), add=TRUE)
        # axis(side = 4, at = seq(1,7) + 5.5, cex.axis=.66, tcl=-.15, label=10^seq(-3,3), las=1, lwd=.5)
        # mtext(side=4, at=14, "Odds ratio", las=2, line=-1,cex=.66)

        image(y = 1:8 +4, x=rep(ncol(logPInt),2)+c(2,3), z=matrix(c(1:8), nrow=1), col=brewer.pal(8,"BrBG"), add=TRUE)
        axis(side = 4, at = seq(1,7) + 4.5, cex.axis=.66, tcl=-.15, label=10^seq(-3,3), las=1, lwd=.5)
        mtext(side=4, at=13, "Odds ratio", las=2, line=-1,cex=.66)


        # Significance position
        # mtext(side=4, at=4, "Significance", las=2, line=-1,cex=.66)
        # points(x=rep(ncol(logPInt),2)+2.5, y=1:2, pch=c("*","."))
        # image(x=rep(ncol(logPInt),2)+c(2,3), y=(2:3) +0.5, z=matrix(1), col=brewer.pal(3,"BrBG"), add=TRUE)
        # mtext(side=4, at=3:1, c("P > 0.05", "FDR < 0.1", "FWER < 0.05"), cex=.66, line=0.2)

        mtext(side=4, at=4, "Significance", las=2, line=-1,cex=.66)
        points(x=rep(ncol(logPInt),2)+2.5, y=1:2, pch=c("*","."))
        image(x=rep(ncol(logPInt),2)+c(2,3), y=(2:3) +0.5, z=matrix(1), col=brewer.pal(3,"BrBG"), add=TRUE)
        mtext(side=4, at=3:1, c("P > 0.05", "FDR < 0.1", "FWER < 0.05"), cex=.66, line=0.2)

        dev.off()
        #####barplot

        return(odds)

}
