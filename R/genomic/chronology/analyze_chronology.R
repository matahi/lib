analyze_chronology <- function(dat.genetics, mutation.df, cutoff.occurrence, out.dir="")
{
        require(ggplot2)

        M <- as.matrix(dat.genetics[,colSums(dat.genetics)>=cutoff.occurrence])
        A <- t(M) %*% M

        co_occ.genes <- which(A>= cutoff.occurrence,arr.ind=T)
        co_occ.once <- sapply(1:nrow(co_occ.genes), function(n)
                              {
                                      co_occ.genes[n,1] < co_occ.genes[n,2]
                              })
        co_occ.list <- co_occ.genes[co_occ.once,]

        genes.order <- lapply(1:nrow(co_occ.list), function(i)
                              {
                                      Gene1 <- colnames(M)[co_occ.list[i,1]]
                                      Gene2 <- colnames(M)[co_occ.list[i,2]]
                                      print(paste0( "Gene1=", Gene1, "; Gene2=", Gene2))

                                      # Samples 
                                      samples.Gene1_2 <- as.character(rownames(dat.genetics)[intersect(which(dat.genetics[,Gene1]==1),which(dat.genetics[,Gene2]==1))])

                                      # What happens if 2 mutations?
                                      order.info <- sapply(1:length(samples.Gene1_2), function(k)
                                                           {

                                                                   dat.Gene1 <- mutation.df[(mutation.df$SampleName==samples.Gene1_2[k])&(mutation.df$Gene==Gene1),]
                                                                   dat.Gene2 <- mutation.df[(mutation.df$SampleName==samples.Gene1_2[k])&(mutation.df$Gene==Gene2),]

                                                                   if (nrow(dat.Gene1)>1|nrow(dat.Gene2)>1)
                                                                   {
                                                                           dat.sample <- mutation.df[mutation.df$SampleName==samples.Gene1_2[k],]

                                                                           out <- "inconclusive"

                                                                           for (i in 1:nrow(dat.Gene1))
                                                                           {
                                                                                   if (all(dat.Gene1$LCI.corr[i]>dat.Gene2$UCI.corr,na.rm=T))
                                                                                   {
                                                                                           out <- paste0(Gene1,"_first")
                                                                                   } 
                                                                           }

                                                                           for (i in 1:nrow(dat.Gene2))
                                                                           {
                                                                                   if (all(dat.Gene2$LCI.corr[i]>dat.Gene1$UCI.corr,na.rm=T))
                                                                                   {
                                                                                           out <- paste0(Gene2,"_first")
                                                                                   } 
                                                                           }

                                                                   } else
                                                                   {
                                                                           if (any(is.na(c(dat.Gene1$LCI.corr, dat.Gene1$UCI.corr, 
                                                                                           dat.Gene2$LCI.corr, dat.Gene2$UCI.corr))))
                                                                           {
                                                                                   out <- NA
                                                                           } else if (dat.Gene1$UCI.corr < dat.Gene2$LCI.corr)
                                                                           {
                                                                                   out <- paste0(Gene2,"_first")
                                                                           } else if (dat.Gene2$UCI.corr < dat.Gene1$LCI.corr)
                                                                           {
                                                                                   out <- paste0(Gene1,"_first")
                                                                           } else
                                                                           {
                                                                                   out <- "inconclusive"
                                                                           }

                                                                   }

                                                                   mutation.df.sample <- mutation.df[(mutation.df$SampleName==samples.Gene1_2[k]),]

                                                                   # pdf(paste0("results/", out.dir, "/", Gene1,"_",Gene2, "_",out,"_",samples.Gene1_2[k],".pdf"))
                                                                   # print(
                                                                   #       ggplot(mutation.df.sample) + geom_point(aes(x=Gene, y=VAF.corr, colour=(Gene %in% c(Gene1,Gene2))),size=5) + ylim(0,1) + 
                                                                   #       geom_segment(aes(x=Gene, xend=Gene, y=LCI.corr, yend=UCI.corr, colour=(Gene %in% c(Gene1,Gene2)))) + theme_bw() +
                                                                   #       theme(legend.position="none")
                                                                   #       #scale_colour_manual(values=col.gene) + theme(legend.position="none")
                                                                   #       )
                                                                   # dev.off()

                                                                   return(out)
                                                           })
                              })

        names(genes.order) <- sapply(1:nrow(co_occ.list), function(i)
                                     {
                                             Gene1 <- colnames(M)[co_occ.list[i,1]]
                                             Gene2 <- colnames(M)[co_occ.list[i,2]]
                                             return(paste0(Gene1,"_",Gene2))
                                     })

        return(genes.order)


}
