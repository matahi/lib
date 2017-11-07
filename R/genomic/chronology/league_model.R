league_model <- function (chronology.df, n.repeat=1000, win=2, tie=1, loss=0, cutoff.occurrence=5, add.bars=T, list.genes=NULL) {

        # tmp <- chronology.df

        # chronology.df <- chronology.df[-which(chronology.df$Gene1 %in% c("MLL_PTD", "FLT3_ITD", "t_15_17")),]
        # chronology.df <- chronology.df[-which(chronology.df$Gene2 %in% c("MLL_PTD", "FLT3_ITD", "t_15_17")),]

        ###
        chronology.df.reverse <- chronology.df
        chronology.df.reverse$Win <- chronology.df$Loss
        chronology.df.reverse$Loss <- chronology.df$Win
        chronology.df.reverse$Gene1 <- chronology.df$Gene2
        chronology.df.reverse$Gene2 <- chronology.df$Gene1
        chronology.df.full <- rbind(chronology.df, chronology.df.reverse)

        ###
        genes <- as.character(unique(chronology.df.full$Gene1))

        ###
        score.gene <- sapply(1:length(genes), function(k)
                             {
                                     chronology.gene <- chronology.df.full[chronology.df.full$Gene1==genes[k],]

                                     if (any(chronology.gene$total <= cutoff.occurrence))
                                             chronology.gene[ chronology.gene$total <= cutoff.occurrence,c("total", "Win","Loss","Inconclusive")] <- c(9,3,3,3 )


                                     score.matrix <- sapply(1:nrow(chronology.gene), function(j)
                                                            {
                                                                    prob.win <- chronology.gene$Win[j]/chronology.gene$total[j]
                                                                    prob.loss <- chronology.gene$Loss[j]/chronology.gene$total[j]
                                                                    prob.tie <- chronology.gene$Inconclusive[j]/chronology.gene$total[j]

                                                                    return(sample(c(win,tie,loss), size=n.repeat, replace=TRUE, prob=c(prob.win,prob.tie, prob.loss)))
                                                            })

                                     score.repeat <- rowSums(score.matrix)

                                     return(score.repeat)
                             })

        # dummy variable
        dummy.score <- rowSums(matrix(sample(c(win,tie,loss), size=(length(genes)-1)*n.repeat, replace=TRUE, prob=c(1/3,1/3, 1/3)),nrow=n.repeat, ncol=length(genes)-1))

        dummy.mean <- -mean(dummy.score)
        dummy.sd <- sd(dummy.score)


        # invert rank
        # score.rank <- t(apply(score.gene,1, function(x){length(genes)-rank(x, ties="random")+1}))
        # colnames(score.rank) <- genes
        # score.rank <- data.frame(score.rank)
        #mean.rank <- apply(score.rank,2,mean)

        colnames(score.gene) <- genes
        score.gene <- -data.frame(score.gene) ## Opposite of score
        mean.score <- apply(score.gene,2,mean)

        library(reshape2)

        gene.m <- melt(score.gene)
        gene.m$variable <- factor(gene.m$variable, levels=genes[order(mean.score,decreasing=T)])
        source("./src/lib/R/plot/misc/data_summary.R")

        pp <- ggplot(gene.m, aes(x=variable,y=value)) + stat_summary(fun.data=data_summary) + xlab("")+ylab("relative order") + coord_flip() +
              theme_bw() + theme(text= element_text(size=25), axis.text.x=element_blank(), axis.ticks=element_blank(), legend.text= element_text(size=20), panel.grid=element_blank()) 

      if (add.bars)
              pp <- pp +geom_hline(aes(yintercept=dummy.mean - dummy.sd), linetype="longdash") + geom_hline(aes(yintercept=dummy.mean + dummy.sd), linetype="longdash")

      ###
      # add colors to x axis
      ###
      ### text.info <- element_text(face = c("bold","plain"), color="black", size=16)

      # text.info <- ifelse(levels(gene.m$variable) %in% list.genes, "bold", "plain")
      # pp <- pp + theme(axis.text.y=element_text(face=text.info))

      if (!is.null(list.genes))
      {
              text.color <- ifelse(levels(gene.m$variable) %in% list.genes, "#2c7bb6", "black")
              pp <- pp + theme(axis.text.y=element_text(color=text.color))
      }


      ###
      return(pp)
}
