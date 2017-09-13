### Reverse the legend
# + guides(fill = guide_legend(reverse=TRUE))

require(ggplot2)

# blank_theme
blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )

# function pie_chart
pie_chart <-  function(dat.df, feature.freq, feature.facet, num_col=NA)
{
        # dat.df <-  Nmut.df
        # feature.freq <- "NMut"
        # feature.facet <- "panel"
        if (is.na(num_col))
                num_col <- length(unique(dat.df[,feature.facet]))

        # summary.dat.df <- table(dat.df[,feature.freq], dat.df[,feature.facet])
        summary.dat.df <- table(dat.df[,feature.freq], dat.df[,feature.facet],useNA="ifany")

        # Significant difference in annotation
        pval <- chisq.test(summary.dat.df)

        n.df <- data.frame(as.matrix(summary.dat.df))
        rownames(n.df) <- NULL

        # Normalize columns
        for (i in 1:ncol(summary.dat.df))
        {
                summary.dat.df[,i]  <- summary.dat.df[,i] /colSums(summary.dat.df)[i]
        }

        #

        summary.df <- data.frame(summary.dat.df)
        if (any(is.na(summary.df$Var1))) # Replace NA by "NA"
        {
                levels.info <-  c(na.omit(levels(summary.df$Var1)), "NA")
                summary.df$Var1 <-  as.character(summary.df$Var1)
                summary.df$Var1[is.na(summary.df$Var1)] <- "NA" 
                summary.df$Var1 <- factor(summary.df$Var1, levels=levels.info)
        }

        colnames(summary.df) <- c("freq","facet","value")

        library(plyr)
        library(dplyr)
        library(scales)
        mydf = ddply(summary.df, .(facet), transform, position = cumsum(value) - 0.5*value) 
        mydf$percent_value <- percent(mydf$value)

        mydf$percent_value[mydf$value <=0.01] <- ""

        mydf <- cbind(mydf, n_total=n.df$Freq)

        # Add the NAs

        # for some reason we have to do that
        mydf$freq <- factor(mydf$freq, levels=rev(levels(mydf$freq)))

        #p <- ggplot(mydf, aes(x="",y=value,fill=freq)) + geom_bar(stat="identity",width=1)  + coord_polar("y",direction=-1) + facet_wrap(~facet, ncol=4) + 
        p <- ggplot(mydf, aes(x="",y=value,fill=freq)) + geom_bar(stat="identity",width=1)  + coord_polar("y",direction=1) + facet_wrap(~facet, ncol=num_col) + 
                                                          blank_theme  + theme(axis.text.x=element_blank())

        if (sum(mydf$value>0.05,na.rm=T)>0)
        {
                p <- p + geom_text(data=mydf[which(mydf$value>0.05),,drop=F],aes(label = percent_value, x="", y=position),color="black", size=5)
        }

        if (sum((mydf$value!=0)&(mydf$value<=0.05), na.rm=T)>0)
        {
                p <- p + geom_text(data=mydf[which((mydf$value!=0)&(mydf$value<=0.05)),,drop=F],aes(label = percent_value, x=1.73, y=position),color="black", size=5)
        }

        if (sum((mydf$value>=0.01)&(mydf$value<=0.05), na.rm=T)>0)
        {
                p <- p + geom_segment(data=mydf[which((mydf$value>=0.01)&(mydf$value<=0.05)),,drop=F],aes(x=1.52, xend=1.65, y=position, yend=position),color="black")
                #geom_text(data=mydf[mydf$value!=0,],aes(label = percent_value, x="", y=position),color="white", size=5)
        }



        return(list(p=p, pval=pval, mydf=mydf))
}


# function pie_chart
pie_chart_old <-  function(dat.df, feature.freq, feature.facet)
{
        # dat.df <- dataset.info
        # feature.freq <- "WHO"
        # feature.facet <- "HDP"

        summary.dat.df <- table(dat.df[,feature.freq], dat.df[,feature.facet])

        # Significant difference in annotation
        pval <- chisq.test(summary.dat.df)

        summary.dat.df[,1]  <- summary.dat.df[,1] /colSums(summary.dat.df)[1]
        summary.dat.df[,2]  <- summary.dat.df[,2] /colSums(summary.dat.df)[2]

        summary.df <- data.frame(summary.dat.df)
        colnames(summary.df) <- c("freq","facet","value")

        library(plyr)
        library(dplyr)
        library(scales)
        mydf = ddply(summary.df, .(facet), transform, position = cumsum(value) - 0.5*value) 

        p <- ggplot(mydf, aes(x="",y=value,fill=freq)) + geom_bar(stat="identity",width=1)  + coord_polar("y") + facet_grid(~facet) + 
                                                          blank_theme  + theme(axis.text.x=element_blank()) + 
                                                          geom_text(data=mydf[mydf$value!=0,],aes(label = percent(value), y=position), size=5)

        return(list(p=p, pval=pval))
}


