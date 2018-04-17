plot_comutation_diff  <- function (dat.hotspot.binary, 
                                      dat.binary, 
                                      OR.range = c(0.01,100), 
                                      disease.type="disease"
                                      ) {

        ############################
        dat.hotspot.binary <- dat.hotspot.binary %>% 
                                        select(4:ncol(dat.hotspot.binary))
        dat.binary <- dat.binary.frequent %>% 
                select(4:ncol(dat.binary.frequent))
        OR.range <- c(0.01, 100)
        ############################

        dat1 <- as.matrix(dat.hotspot.binary)
        dat2 <- as.matrix(dat.binary)

        pp <- try(compare_comutation(dat1=as.matrix(dat.hotspot.binary),
                                     dat2=as.matrix(dat.binary),
                                     OR.range=OR.range),
                  silent=T)

        
        label.OR <- as.character(10^seq(log10(OR.range[1]), log10(OR.range[2])))
        
        # Data.frame
        pp <- data.frame(pp)
        pp.m <- melt(pp) 
        pp.m <- cbind(pp.m, hotspot = rownames(pp)) %>% tbl_df()
        pp.m <- pp.m %>% 
                mutate(hotspot= as.character(hotspot))
        # pp.m <- left_join(pp.m, hotspot.summary, by="hotspot")
        # pp.m <- pp.m %>% mutate(info = paste0(hotspot, 
        #                                       " (n=",
        #                                       n_hotspot,
        #                                       ")"))

        colnames(pp.m) <- c("genes", "OR", "hotspot")
        
        if (length(unique(pp.m$OR))==1)
                pp.m$OR[1] <- pp.m$OR[1]+0.01
        
        
        p1 <- ggplot(pp.m) + 
                geom_tile(aes_string(x="genes", 
                                     y=1, 
                                     fill="OR"),
                          colour="grey20") + 
                theme_minimal() + 
                theme(axis.title.x=element_blank(), 
                      panel.border=element_blank(), 
                      panel.grid=element_blank(), 
                      axis.ticks=element_blank(), 
                      plot.title=element_text(size=14,face="bold") ,  
                      axis.text.x=element_text(angle=90,hjust=1, vjust=0.5)) + 
                scale_fill_gradientn(colours=brewer.pal(length(label.OR),"RdBu"), 
                             name="Odd Ratio", 
                             labels=label.OR, 
                             breaks=1:length(label.OR), 
                             limits=c(1,length(label.OR))) +
                ylab("IDH1")
        

        pp <- p1
        
                # Size
        pp <- pp + theme(text = element_text(size=20))
        
        #return(list(pp.m, p.sign))
        return(pp)
        
}


