plot_precedence <- function(chronology.df, gene.list, Reordering=NULL)
{
        cutoff.occurrence <- 20
        # TODO: PUT RANDOM IF total < 5

        gene.list <- targets
        # Ratio wins/losses
        Ratio.win <- sapply(1:nrow(chronology.df), function(k)
                            {
                                    if (chronology.df$Win[k] + chronology.df$Loss[k] > 0)
                                            return(chronology.df$Win[k]/(chronology.df$Win[k] + chronology.df$Loss[k]))

                                    return(NA)
                            })
        tmp.df <- cbind(wins=Ratio.win, chronology.df)

        # Ratio Conclusive/Inconclusives
        Ratio.conclusive <- sapply(1:nrow(chronology.df), function(k)
                            {
                                    if (chronology.df$Win[k] + chronology.df$Loss[k] > 0)
                                            return((chronology.df$Win[k] + chronology.df$Loss[k])/(chronology.df$total[k]))

                                    return(NA)
                            })
        tmp.df <- cbind(wins=Ratio.win, chronology.df)


        library(reshape2)

        ##
        # i) Ratio win loss
        wins.df <- dcast(tmp.df, Gene1 ~ Gene2, value.var="wins")
        rownames(wins.df) <- wins.df[,1]
        wins.df <- wins.df[,-1]
        wins.df <- cbind(NA,wins.df)
        colnames(wins.df)[1] <- gene.list[1]
        #TODO: check for missing
        wins.df <- wins.df[gene.list, match(gene.list,colnames(wins.df))]
        wins.df <- as.matrix(wins.df)

        ##
        # ii) Ratio win loss
        inconclusive.df <- dcast(tmp.df, Gene1 ~ Gene2, value.var="Inconclusive")
        rownames(inconclusive.df) <- inconclusive.df[,1]
        inconclusive.df <- inconclusive.df[,-1]
        inconclusive.df <- cbind(NA,inconclusive.df)
        colnames(inconclusive.df)[1] <- gene.list[1]
        #TODO: check for missing
        inconclusive.df <- inconclusive.df[gene.list, match(gene.list,colnames(inconclusive.df))]
        inconclusive.df <- t(inconclusive.df)

        ##
        # iii) Occurrence
        occ.df <- dcast(tmp.df, Gene1 ~ Gene2, value.var="total")
        rownames(occ.df) <- occ.df[,1]
        occ.df <- occ.df[,-1]
        occ.df <- cbind(NA,occ.df)
        colnames(occ.df)[1] <- gene.list[1]
        #TODO: check for missing
        occ.df <- occ.df[gene.list, match(gene.list,colnames(occ.df))]
        occ.df <- as.matrix(occ.df)


        #############################################
        #############################################
        # What is pairs??????
        ########### HEATMAP ALL 
        library(RColorBrewer)
        library(gplots)

        par(bty="n", mgp = c(2,.5,0), mar=c(4,4,4,4)+.5, las=2, tcl=-.33)
        ix = TRUE#colnames(interactions) %in% colnames(all_genotypes)
        M <-  matrix( NA, ncol=ncol(wins.df), nrow=nrow(wins.df))

        ####
        M[lower.tri(M,diag=T)] <- cut(inconclusive.df[lower.tri(inconclusive.df,diag=T)], breaks=c(-1,0,5,10,20,50,100), include.lowest=T)

        # breaks = c(-4:0-.Machine$double.eps,0:4)  
        # Breaks is size 10

        M[upper.tri(M)] <- as.numeric(cut(wins.df[upper.tri(wins.df)], breaks= c(-.Machine$double.eps, seq(0,1,0.2))   )) + 8 # TODO: Look at how many points

        # toto <- M[upper.tri(M, diag=TRUE)]
        # tmp <- pairs[o,o][upper.tri(M, diag=TRUE)]*nrow(dat.integrate)

        # Size is 7

        ## Need to print max(M) then adjust colours accordingly

        colorset <- c(
                      brewer.pal(7, "Greens"), #BrBG
                      "grey90",
                      brewer.pal((max(na.omit(M[upper.tri(M,diag=T)]))-8), "RdBu")) #Blues

        M[is.na(M)] <- 8

        # print(max(M))
        image(x=1:ncol(wins.df), y=1:nrow(wins.df), t(M)[,nrow(M):1], col=colorset, breaks=0:max(M,na.rm=TRUE), xaxt="n", yaxt="n", xlab="",ylab="", xlim=c(0, ncol(wins.df)+3), ylim=c(0, ncol(wins.df)+3),
              main="Pairwise Precedences")
        l <- colnames(wins.df)

        if (is.null(Reordering))
        {
                mtext(side=1, at=1:ncol(wins.df), l, cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1))
                mtext(side=2, at=1:ncol(wins.df), rev(l), cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1))
        } else
        {
                mtext(side=1, at=1:ncol(wins.df), l, cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1), col=Reordering$color.features)
                kmtext(side=2, at=1:ncol(wins.df), rev(l), cex=.66, font=ifelse(grepl("^[A-Z]",l),3,1), col=Reordering$color.features)
        }

        ######################

        for (k in 0:ncol(wins.df))
        {
                # horizontal
                segments(x0=k+.5, x1=k+.5, y0=.5, y1=ncol(wins.df)+.5, col="grey", lwd=.5)
                # vertical
                segments(x0=.5, x1=ncol(wins.df)+.5, y0=k+.5, y1=k+.5, col="grey", lwd=.5)
        }


        w = arrayInd(which( occ.df > cutoff.occurrence), rep(nrow(occ.df),2))
        inv_w <- w
        inv_w[,2] <- ncol(occ.df) - w[,1] + 1
        inv_w[,1] <- w[,2]
        points(inv_w, pch="*", col="salmon3")

        # w = arrayInd(which(p.adjust(P) < .05), rep(nrow(wins.df),2))
        # points(w, pch="*", col="salmon3")

        ######################


        offset.val <- 8.5

        image(y = 1:7 + offset.val, x=ncol(wins.df) +c(2,3), z=matrix(c(1:7), nrow=1), col=c("white",brewer.pal(max(na.omit(M[upper.tri(M,diag=T)]))-10, "Greens")), add=TRUE)
        axis(side = 4, at = seq(1,7) + offset.val + .5, cex.axis=.66, tcl=-.15, label=c(1,5,10,20,50,100,200), las=1, lwd=.5)
        mtext(side=4, at=8 + offset.val, "Inconclusive cases", las=2, line=-1,cex=.66)

        offset.precedence <- 2.5

        image(y = 1:5 + offset.precedence, x=ncol(wins.df)+c(2,3), z=matrix(c(1:5), nrow=1), col=brewer.pal(11,"RdBu")[seq(1,11,2)], add=TRUE)
        axis(side = 4, at = seq(0,5) + offset.precedence + .5, cex.axis=.66, tcl=-.15, label=seq(0,1,0.2)*100, las=1, lwd=.5)
        mtext(side=4, at=6 + offset.precedence, "Precedence (%)", las=2, line=-1,cex=.66)

        offset.text <- -1.5

        mtext(side=4, at=3 + offset.text, "Significance", las=2, line=-1,cex=.66)
        points(x=ncol(wins.df)+2.5, y=2 + offset.text, pch="*")
        mtext(side=4, at=2 + offset.text, paste0("Occurrence > ",cutoff.occurrence), cex=.66, line=0.2)

}
