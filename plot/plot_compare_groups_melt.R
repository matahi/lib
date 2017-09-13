plot_compare_groups_melt <- function(dat.list, threshold=0.01)
{
        # dat.list <- mut.onc
        # dat.list <-  full.df.common
        # dat.list <- karyotypes 
        # features <- NULL
        # threshold <- 0.01

        require(ggplot2)
        require(scales)
        require(dplyr)
        require(reshape2)


        # transform
        features <- Reduce("union", lapply(dat.list,function(x){unique(x$Gene)})) 

        dat.common <- lapply(1:length(dat.list), function(n)
                            {
                                    tmp <- dcast(dat.list[[n]], Sample ~ Gene)
                                    tmp <- tmp[,-1] # Remove sample names
                                    tmp[tmp>0] <- 1

                                    final.df <- data.frame(matrix(NA, nrow=nrow(tmp), ncol=length(features)))
                                    colnames(final.df) <- features

                                    for (i in 1:ncol(tmp))
                                    {
                                            final.df[,colnames(tmp)[i]] <- tmp[,i]
                                    }

                                    return(final.df)
                            })



        dat.processed <- lapply(1:length(dat.common), function(n)
                                {
                                        #n_event <- colSums(dat.common[[n]])
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


        dat.frequent <- dat.final[dat.final$features %in% names(features.info)[features.info=="frequent"],]
        dat.rare <- dat.final[dat.final$features %in% names(features.info)[features.info=="rare"],]

                # 
        p <- ggplot(dat.final) + geom_bar(aes(x=features,y=frequency,fill=dataset), stat="identity", position="dodge") + scale_y_continuous(labels=percent) + ylab("Frequency") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top") 

        p.frequent <- ggplot(dat.frequent) + geom_bar(aes(x=features,y=frequency,fill=dataset), stat="identity", position="dodge") + scale_y_continuous(labels=percent) + ylab("Frequency") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top") 
        p.rare <- ggplot(dat.rare) + geom_bar(aes(x=features,y=frequency,fill=dataset), stat="identity", position="dodge") + scale_y_continuous(labels=percent) + ylab("Frequency") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top") 

        # Order
        dat.opp <- Reduce('rbind', dat.processed.opp) 
        order_var <- order(dat.processed.opp[[1]]$frequency,decreasing=T)
        dat.opp$features <- factor(dat.opp$features,levels=dat.processed.opp[[1]]$features[order_var])

        dat.1 <- subset(dat.opp, frequency>=0)
        dat.2 <- subset(dat.opp, frequency<=0)

        p.opp <- ggplot() + geom_bar(dat=dat.1,aes(x=features,y=frequency,fill=dataset), stat="identity") + 
                            geom_bar(dat=dat.2,aes(x=features,y=frequency,fill=dataset), stat="identity") + 
                            scale_y_continuous(labels=percent) + ylab("Frequency") + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top") 

        return(list(p=p, p.opp=p.opp, p.frequent=p.frequent, p.rare=p.rare, sig.features=sig.features,dat.final=dat.final))
}
