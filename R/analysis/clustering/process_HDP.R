process_HDP <- function (dat, HDP.posterior, HDP.output ) {

        # Todo: update
        require(hdp)
        distn <- comp_dp_distn(HDP.posterior)$mean[-1,]

        HDP.class <- sapply(1:nrow(distn), function(n)
                            {
                                    if (all(is.na(distn[n,])))
                                    {
                                            return(1)
                                    } else
                                    {
                                            return(which.max(distn[n,])+1)
                                    }
                            })
        n.HDP <- length(unique(HDP.class))

        HDP.post <- sapply(1:nrow(distn), function(n)
                           {
                                   if (all(is.na(distn[n,])))
                                   {
                                           return(1)
                                   } else
                                   {
                                           return(max(distn[n,]))
                                   }
                           })

        #### Reorder sample in classes 
        HDP.classes <- sort(unique(HDP.class))

        Index.bis <- lapply(1:length(HDP.classes), function(n)
                            {
                                    dat.bis <- data.frame(dat[which(HDP.class==HDP.classes[n]),,drop=F])

                                    order.gene <- order(colSums(dat.bis),decreasing=T)

                                    gene.order <- paste(paste0("-",colnames(dat.bis)[order.gene[1:10]]), collapse=",")

                                    Index.samples <- do.call(order,-dat.bis[,order.gene])
                            })

        final.order <- Reduce("c", lapply(1:length(HDP.classes), function(n)
                                          {
                                                  class.index <- which(HDP.class==HDP.classes[n])

                                                  index.final <- class.index[Index.bis[[n]]]
                                          }))

        return(list(class=HDP.class, order=final.order, prob=distn))
        
}
