plot_bar_groups <- function(dat,features=NULL,type="frequency",group=T)
{
        # dat <- karyotype.noXY
        # features <- NULL

        require(ggplot2)
        require(scales)
        require(dplyr)
        require(reshape2)

        if (!is.null(features))
        {
                dat.features <- dat[,features]
        } else
        {
                dat.features <- dat
        }


        n_event <- rowSums(dat)
        names(n_event) <- NULL
        n_sample <- nrow(dat)

        # create the dataframe
        tmp <- melt(dat.features)
        tmp.df <- data.frame(tmp, n_event,n_sample)
        tmp.ter <- tmp.df[which(tmp.df$value==1),]

        # ordering variables according to frequency
        freq <- colSums(dat.features)
        order_var <- order(freq,decreasing=T)
        tmp.ter$variable <- factor(tmp.ter$variable,levels=names(freq)[order_var])

        # Processing of the data
        event.summary <- sapply(1:nrow(tmp.ter), function(n)
                                {
                                        if (tmp.ter$n_event[n]==1)
                                        {
                                                return("Isolated")
                                        } else if (tmp.ter$n_event[n]==2)
                                        {
                                                return("+1")
                                        } else
                                        {
                                                return("Complex")
                                        }
                                })
        tmp.ter <- cbind(tmp.ter, summary=event.summary)
        tmp.ter$summary <- factor(tmp.ter$summary, levels=c("Isolated","+1","Complex"))

        if (type=="frequency")
        {
                if (group)
                {
                        p <- ggplot(tmp.ter) + geom_bar(aes(x=variable,n_sample=n_sample,y=..count../n_sample,fill=factor(summary))) + scale_y_continuous(labels=percent) + ylab("Frequency")
                } else
                {
                        p <- ggplot(tmp.ter) + geom_bar(aes(x=variable,n_sample=n_sample,y=..count../n_sample)) + scale_y_continuous(labels=percent) + ylab("Frequency")
                }
        } else if (type=="count")
        {
                if (group)
                {
                        p <- ggplot(tmp.ter) + geom_bar(aes(x=variable,fill=factor(summary))) 
                } else
                {
                        p <- ggplot(tmp.ter) + geom_bar(aes(x=variable)) 
                }

        }

        return(p)
}
