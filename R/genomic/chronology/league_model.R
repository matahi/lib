league_model <- function (chronology.df, n.repeat=1000, type="score", win=2, tie=1, loss=0, cutoff.occurrence=5, add.bars=T, group.genes=NULL) {

        # cutoff.occurrence <- 5
        # n.repeat <- 1000
        # win <- 2
        # tie <- 1
        # loss <- 0
        # add.bars <- T
        # group.genes <- Reordering.features

        require(tidyverse)

        ###
        genes <- unique(chronology.df$Gene1)

        ###
        precedence.df <- chronology.df %>% filter(Gene1 != Gene2)

        ###
        score.gene <- sapply(1:length(genes), function(k)
                             {
                                     chronology.gene <- precedence.df[precedence.df$Gene1==genes[k],]

                                     if (any(chronology.gene$Total <= cutoff.occurrence))
                                             chronology.gene[ chronology.gene$Total <= cutoff.occurrence,c("Total", "Win","Loss","Inconclusive")] <- c(9,3,3,3)

                                     score.matrix <- sapply(1:nrow(chronology.gene), function(j)
                                                            {
                                                                    prob.win <- chronology.gene$Win[j]/chronology.gene$Total[j]
                                                                    prob.loss <- chronology.gene$Loss[j]/chronology.gene$Total[j]
                                                                    prob.tie <- chronology.gene$Inconclusive[j]/chronology.gene$Total[j]

                                                                    return(sample(c(win,tie,loss), size=n.repeat, replace=TRUE, prob=c(prob.win,prob.tie, prob.loss)))
                                                            })

                                     score.repeat <- rowSums(score.matrix)

                                     return(score.repeat)
                             })

        # dummy variable
        dummy.score <- rowSums(matrix(sample(c(win,tie,loss), size=(length(genes)-1)*n.repeat, replace=TRUE, prob=c(1/3,1/3, 1/3)),nrow=n.repeat, ncol=length(genes)-1))

        dummy.mean <- -mean(dummy.score)
        dummy.sd <- sd(dummy.score)

        # # invert rank
        # if (type == "rank")
        # {
        # score.rank <- t(apply(score.gene,1, function(x){length(genes)-rank(x, ties="random")+1}))
        # colnames(score.rank) <- genes
        # }
        # score.rank <- data.frame(score.rank)
        #mean.rank <- apply(score.rank,2,mean)

        colnames(score.gene) <- genes
        score.gene <- -data.frame(score.gene) %>% tbl_df() ## Opposite of score
        mean.score <- apply(score.gene,2,mean)

        score.df <- gather(score.gene,
                           key = "alteration",
                           value = "score")

        gene.m <- melt(score.gene)
        gene.m$variable <- factor(gene.m$variable, levels=genes[order(mean.score,decreasing=T)])

        score.df <- score.df %>%
                mutate(alteration = factor(alteration, levels = names(mean.score)[order(mean.score, decreasing=T)]))

        source("./src/lib/R/plot/misc/data_summary.R")

        pp <- ggplot(score.df, aes(x=alteration,y=score)) +
                stat_summary(fun.data=data_summary) +
                xlab("") +
                ylab("relative order") +
                coord_flip() +
                theme(axis.text.x = element_blank(),
                      axis.ticks = element_blank(),
                      legend.text = element_text(size=20),
                      panel.grid = element_blank())

      if (add.bars)
              pp <- pp +
                      geom_hline(aes(yintercept=dummy.mean - dummy.sd), linetype="longdash") + 
                      geom_hline(aes(yintercept=dummy.mean + dummy.sd), linetype="longdash")

      if (!is.null(group.genes))
      {
              cols.axis <- group.genes$colours[ match(levels(score.df$alteration), group.genes$alterations )] 
              pp <- pp + theme(axis.text.y=element_text(color=cols.axis))
      }

      ###
      league.out <- list(ranking = score.df, p = pp)

      return(league.out)
}
