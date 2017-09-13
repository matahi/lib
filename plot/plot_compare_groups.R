plot_compare_groups <- function(dat.list,features=NULL, add.significance=F, threshold=0.005, max_num=100)
{
        # dat.list <- cyto.list 
        # dat.list <- genetics.list
        # dat.list <- CNA.list
        # dat.list <- dat.list.1
        # dat.list <- dat.binary.list
        # features <- NULL
        # threshold <- 0.01
        # add.significance <- F

        
        require(ggplot2)
        require(scales)
        require(dplyr)
        require(reshape2)

        if (!is.null(features))
        {
                dat.common <- lapply(1:length(dat.list), function(n)
                                     {
                                             dat.feature <- data.frame(matrix(0, nrow=nrow(dat.list[[n]]), ncol=length(features)))
                                             colnames(dat.feature) <- features

                                             rownames(dat.feature) <- rownames(dat.list[[n]])

                                             dat.feature[,colnames(dat.list[[n]])] <- dat.list[[n]]

                                             return(dat.feature)
                                             #return(dat.list[[n]][,features])
                                     })
        } else
        {
                #features <- Reduce("intersect",lapply(dat.list,colnames)) 
                features <- Reduce("union",lapply(dat.list,colnames)) 

                # features <- Reduce("intersect", lapply(1:length(dat.list), function(n){names(which(colSums(dat.list[[n]],na.rm=T)>0))})) 
                dat.common <- lapply(1:length(dat.list), function(n)
                                     {
                                             dat.df <- data.frame(matrix(0,nrow=nrow(dat.list[[n]]), ncol=length(features)))
                                             rownames(dat.df) <- rownames(dat.list[[n]])
                                             colnames(dat.df) <- features

                                             dat.df[, colnames(dat.list[[n]])] <- dat.list[[n]]

                                             return(dat.df)
                                     })

        }

        dat.processed <- lapply(1:length(dat.common), function(n)
                                {
                                        n_event <- colSums(dat.common[[n]],na.rm=T)
                                        n_NA <- colSums(is.na(dat.common[[n]]))
                                        n_sample <- nrow(dat.common[[n]])

                                        tmp.df <- data.frame(features=names(n_event), n_event,n_sample, frequency=n_event/n_sample, n_NA, dataset=names(dat.list)[n])

                                        # order features
                                        tmp.df$features <- factor(tmp.df$features, levels=tmp.df$features)

                                        return(tmp.df)
                                })

        dat.processed.opp <- lapply(1:length(dat.common), function(n)
                                {
                                        n_event <- colSums(dat.common[[n]], na.rm=T)
                                        n_NA <- colSums(is.na(dat.common[[n]]))
                                        n_sample <- nrow(dat.common[[n]])

                                        if (n==1)
                                        {
                                                tmp.df <- data.frame(features=names(n_event), n_event,n_sample, frequency=n_event/n_sample, n_NA, dataset=names(dat.list)[n])
                                        } else
                                        {
                                                tmp.df <- data.frame(features=names(n_event), n_event,n_sample, frequency=-n_event/n_sample, n_NA, dataset=names(dat.list)[n])
                                        }

                                        # order features
                                        tmp.df$features <- factor(tmp.df$features, levels=tmp.df$features)

                                        return(tmp.df)
                                })

        # features with significantly different ratios
        sig.features <- sapply(1:length(features), function(n)
                               {
                                       Dat1.pos <- dat.processed[[1]]$n_event[n]
                                       Dat2.pos <- dat.processed[[2]]$n_event[n]

                                       Dat1.tot <- dat.processed[[1]]$n_sample[n]  - dat.processed[[1]]$n_NA[n]
                                       Dat2.tot <- dat.processed[[2]]$n_sample[n] - dat.processed[[2]]$n_NA[n]

                                       return(fisher.test(rbind(c(Dat1.pos, Dat1.tot - Dat1.pos),
                                                                c(Dat2.pos, Dat2.tot - Dat2.pos)))$p.value)
                               })
        sig.features <- p.adjust(sig.features, method="fdr")
        names(sig.features) <- features

        # ordering variables according to frequency
        dat.final <- Reduce('rbind', dat.processed)

        features.info <- sapply(1:length(features), function(n)
                                {
                                        freqs <- dat.final$frequency[dat.final$features==features[n]]

                                        if (any(freqs>threshold))
                                        {
                                                return("frequent")
                                        } else
                                        {
                                                return("rare")
                                        }
                                })
        names(features.info) <- features

        # Data.frame
        dat.frequent <- dat.final[dat.final$features %in% names(features.info)[features.info=="frequent"],]
        dat.rare <- dat.final[dat.final$features %in% names(features.info)[features.info=="rare"],]
        dat.opp <- Reduce('rbind', dat.processed.opp) 

        # Order
        order_var <- order(dat.processed.opp[[1]]$frequency,decreasing=T)
        dat.final$features <- factor(dat.final$features,levels=dat.processed.opp[[1]]$features[order_var])
        dat.frequent$features <- factor(dat.frequent$features,levels=dat.processed.opp[[1]]$features[order_var])
        dat.rare$features <- factor(dat.rare$features,levels=dat.processed.opp[[1]]$features[order_var])
        dat.opp$features <- factor(dat.opp$features,levels=dat.processed.opp[[1]]$features[order_var])

        # Plots
        p <- ggplot(dat.final) + geom_bar(aes(x=features,y=frequency,fill=dataset), stat="identity", position="dodge") + scale_y_continuous(labels=percent) + ylab("Frequency") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="top") 

        p.frequent <- ggplot(dat.frequent) + geom_bar(aes(x=features,y=frequency,fill=dataset), stat="identity", position="dodge") + scale_y_continuous(labels=percent) + ylab("Frequency") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top") 

        p.rare <- ggplot(dat.rare) + geom_bar(aes(x=features,y=frequency,fill=dataset), stat="identity", position="dodge") + scale_y_continuous(labels=percent) + ylab("Frequency") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top") 

        # dat.1 <- subset(dat.opp, frequency>0)
        # dat.2 <- subset(dat.opp, frequency<0)

        dat.1 <- dat.opp[dat.opp$dataset==unique(dat.opp$dataset)[1],]
        dat.2 <- dat.opp[dat.opp$dataset==unique(dat.opp$dataset)[2],]

        ymax <- max(dat.1$frequency)
        ymin <- min(dat.2$frequency)

        yfinal <- max(ymax,-ymin) + 0.05

        p.opp <- ggplot() + geom_bar(dat=dat.1,aes(x=features,y=frequency,fill=dataset), stat="identity", position="identity") + 
        geom_bar(dat=dat.2,aes(x=features,y=frequency,fill=dataset), stat="identity", position="identity") + 
        scale_y_continuous(labels=percent,limits=c(-yfinal,yfinal)) + ylab("Frequency") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top") 

        if (add.significance)
        {
                significant.features <- sig.features[which(sig.features < 0.05)]

                if (!length(significant.features)==0)
                {
                        for (k in 1:length(significant.features))

                        {
                                #print(k)
                                # Recover position
                                x.pos <- which(levels(dat.final$features)==names(significant.features)[k])

                                # Star
                                if (significant.features[k] <0.001)
                                {
                                        star <- "***"
                                } else if (significant.features[k] < 0.01)
                                {
                                        star <- "**"
                                } else 
                                {
                                        star <- "*"
                                } 

                                # y level
                                y.pos <- max(dat.final$frequency[dat.final$features==names(significant.features)[k]])

                                p <- p + annotate("text", x=x.pos, y=y.pos+0.02, label=star) +
                                geom_segment(x=x.pos-0.25 , xend=x.pos+0.25, y=y.pos+0.01, yend=y.pos+0.01)

                                y.pos.opp <- dat.1$frequency[dat.1$features==names(significant.features)[k]]
                                p.opp <- p.opp + annotate("text", x=x.pos, y=y.pos.opp +0.02, label=star)
                        }

                }

        }

        return(list(p=p, p.opp=p.opp, p.frequent=p.frequent, p.rare=p.rare, sig.features=sig.features,dat.final=dat.final))
}
