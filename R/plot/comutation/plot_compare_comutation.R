compare_comutation <- function (dat1, dat2, OR.range=c(0.01,100)) {

                # Odds ratio
                odds <- sapply(1:ncol(dat2), 
                               function(i) sapply(1:ncol(dat1), 
                                                  function(j) {f<- try(fisher.test(table(dat2[,i], 
                                                                                         dat1[,j])), 
                                                                       silent=TRUE); 
                                                  if(class(f)=="try-error") 
                                                          f=NA 
                                                  else f$estimate} ))

                # if dat1 has only 1 row
                odds <- matrix(odds, nrow = ncol(dat1), ncol = ncol(dat2))

                # Pvals
                logPInt <- sapply(1:ncol(dat2), 
                                  function(i) sapply(1:ncol(dat1), 
                                                     function(j) {f<- try(fisher.test(dat2[,i], 
                                                                                      dat1[,j]), 
                                                                          silent=TRUE); 
                                                     if(class(f)=="try-error") 
                                                             0 
                                                     else ifelse(f$estimate>1, 
                                                                 -log10(f$p.val),
                                                                 log10(f$p.val))} ))

                # if dat1 has only 1 row
                logPInt <- matrix(logPInt, nrow = ncol(dat1), ncol = ncol(dat2))

                # 
                rownames(logPInt) <- rownames(odds) <- colnames(dat1)
                colnames(logPInt) <- colnames(odds) <- colnames(dat2)

                # Thresholding the odd ratios to stand between the OR.ranges
                odds[is.infinite(odds)] <- NA

                odds[odds<10^(log10(OR.range)[1]+1)] <- OR.range[1]
                odds[odds>10^(log10(OR.range)[2]-1)] <- OR.range[2]

                odds[10^-abs(logPInt) > 0.05] = 1

                ## 
                logOdds=log10(odds)
                M <-  matrix( NA, ncol=ncol(odds), nrow=nrow(odds))
                M <- cut(logOdds, 
                         breaks = c(log10(OR.range)[1]:0-.Machine$double.eps,0:log10(OR.range)[2]), 
                         include.lowest=TRUE)  
                M <- matrix(as.numeric(M), ncol=ncol(odds), nrow=nrow(odds))

                rownames(M) <- rownames(odds)
                colnames(M) <- colnames(odds)

                # Output
                return(M)
}


odds_significance <- function (dat.hotspot.binary, dat.binary) {

        library(DescTools)
        library(epiR)

        odds <- sapply(1:ncol(dat.binary), 
                       function(i) sapply(1:ncol(dat.hotspot.binary), 
                                          function(j) {f<- try(fisher.test(table(dat.binary[,i], 
                                                                                 dat.hotspot.binary[,j])), 
                                                               silent=TRUE); 
                                          if(class(f)=="try-error") 
                                                  f=NA 
                                          else 
                                                  f$estimate} ))
        
        
        
        # Look at 2 by 2 differences in odd ratios
        ncomb <- combn(1:ncol(dat.hotspot.binary),2)
        
        p.sign <- Reduce("rbind", lapply(1:ncol(ncomb), function(combi)
                                       {
                                               # BDT
                                               pvals <- sapply(1:ncol(dat.binary), function(k)
                                                               {
                                                                       dat.table <- as.table(array(c(table(dat.binary[,k], 
                                                                                                           dat.hotspot.binary[,ncomb[1,combi]]),
                                                                                                     table(dat.binary[,k], 
                                                                                                           dat.hotspot.binary[,ncomb[2,combi]])), 
                                                                                                   dim=c(2,2,2)
                                                                                                   ))
        
                                                                       tt <- try(BreslowDayTest(dat.table), 
                                                                                 silent=T)
        
                                                                       if (class(tt)=="try-error")
                                                                               return(NA)
        
                                                                       return(BreslowDayTest(dat.table)$p.val)
                                                               })
                                               pvalues.bdt <- p.adjust(pvals)
        
                                               return(data.frame(hotspot1=ncomb[1,combi],
                                                                 hotspot2=ncomb[2,combi],
                                                                 gene = 1:ncol(dat.binary),
                                                                 gene_name = colnames(dat.binary)[1:ncol(dat.binary)],
                                                                 pval = pvalues.bdt))
                                       }))
        
        # p.sign %>% tbl_df() %>% filter(pval < 0.05)

        return(p.sign)
}


plot_compare_comutation  <- function (dat.hotspot.binary, dat.binary, OR.range = c(0.01,100), disease.type="disease", hotspot.orders=NULL) {

        pp <- try(compare_comutation(dat1=as.matrix(dat.hotspot.binary),
                                     dat2=as.matrix(dat.binary),
                                     OR.range=OR.range),
                  silent=T)
        
        label.OR <- as.character(10^seq(log10(OR.range[1]), log10(OR.range[2])))
        
        hotspot.summary <- data.frame(hotspot = colnames(dat.hotspot.binary),
                                      n_hotspot = dat.hotspot.binary %>% 
                                              colSums()) %>% 
                                         tbl_df() %>%
                                         mutate(hotspot = as.character(hotspot),
                                                n_hotspot = as.integer(n_hotspot) )
        
        # Calculate significance

        if (!is.null(hotspot.orders))
        {
                p.sign <- odds_significance(as.matrix(dat.hotspot.binary)[,hotspot.orders], 
                                            as.matrix(dat.binary))
        } else
        {
                p.sign <- odds_significance(as.matrix(dat.hotspot.binary), 
                                            as.matrix(dat.binary))
        }

        # Data.frame
        pp <- data.frame(pp)
        pp.m <- melt(pp) 
        pp.m <- cbind(pp.m, hotspot = rownames(pp)) %>% tbl_df()
        pp.m <- pp.m %>% 
                mutate(hotspot= as.character(hotspot))
        pp.m <- left_join(pp.m, hotspot.summary, by="hotspot")
        pp.m <- pp.m %>% mutate(info = paste0(hotspot, 
                                              " (n=",
                                              n_hotspot,
                                              ")"))


        if (!is.null(hotspot.orders))
        {
                pp.m <- pp.m %>%
                        mutate(hotspot=factor(hotspot, levels=hotspot.orders))

                pp.m <- pp.m %>%
                        mutate(info = factor(info, levels= pp.m %>%
                                             group_by(hotspot) %>%
                                                     slice(1) %>%
                                                     .$info))

        }

        colnames(pp.m) <- c("genes", "OR", "hotspot", "n_hotspot", disease.type)
        
        if (length(unique(pp.m$OR))==1)
                pp.m$OR[1] <- pp.m$OR[1]+0.01
        
        
        p1 <- ggplot(pp.m) + 
                geom_tile(aes_string(x="genes", 
                                     y=disease.type, 
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
                ylab(paste0(disease.type, 
                               " (n=",
                               dat.hotspot.binary %>%
                                       nrow(),
                               ")"))
        
        
        pp <- p1 + 
                geom_rect(data=p.sign %>% filter(pval < 0.05),
                       aes(xmin=gene-0.5, 
                           xmax=gene+0.5,
                           ymin=hotspot1-0.5,
                           ymax=hotspot2+0.5
                           ),
                       fill = NA,
                       colour="black",
                       size=1) +
                geom_point(data=p.sign %>% filter(pval < 0.05),
                           aes(x=gene,
                               y=hotspot1+0.5),
                           shape=8,
                           size=5,
                           colour="black")

        # Size
        pp <- pp + theme(text = element_text(size=20))
        
        #return(list(pp.m, p.sign))
        return(pp)
        
}


