analyze_chronology <- function(dat.binary, dat.genetics, Gene1, Gene2)
{
        # system("mkdir -p results/supp/Ross/chronology")
        # dat.genetics <- dat.genetics.bis
        # Gene1 <- "KRAS"
        # Gene2 <- "NPM1"

        dat.binary <- dat
        # dat.genetics <- dat.genetics.bis 
        Gene1 <- "NRAS"
        Gene2 <- "SF3B1"

        # Samples 
        samples.Gene1_2 <- as.character(rownames(dat.binary)[intersect(which(dat.binary[,Gene1]==1),which(dat.binary[,Gene2]==1))])

        # What happens if 2 mutations?
        order.info <- sapply(1:length(samples.Gene1_2), function(k)
                             {
                                     dat.Gene1 <- dat.genetics[(dat.genetics$ID==samples.Gene1_2[k])&(dat.genetics$Gene==Gene1),]
                                     dat.Gene2 <- dat.genetics[(dat.genetics$ID==samples.Gene1_2[k])&(dat.genetics$Gene==Gene2),]

                                     if (nrow(dat.Gene1)==0|nrow(dat.Gene2)==0)
                                     {
                                             print("Careful binary matrix and VCF files do not correspond")
                                             out <- NA
                                     } else if (nrow(dat.Gene1)>1|nrow(dat.Gene2)>1)
                                     {
                                             dat.sample <- dat.genetics[dat.genetics$ID==samples.Gene1_2[k],]

                                             out <- "inconclusive"

                                             for (i in 1:nrow(dat.Gene1))
                                             {
                                                     if (all(dat.Gene1$LCI[i]>dat.Gene2$UCI,na.rm=T))
                                                     {
                                                             out <- paste0(Gene1," first")
                                                     } 
                                             }

                                             for (i in 1:nrow(dat.Gene2))
                                             {
                                                     if (all(dat.Gene2$LCI[i]>dat.Gene1$UCI,na.rm=T))
                                                     {
                                                             out <- paste0(Gene2," first")
                                                     } 
                                             }

                                     } else
                                     {
                                             if (any(is.na(c(dat.Gene1$LCI, dat.Gene1$UCI, 
                                                             dat.Gene2$LCI, dat.Gene2$UCI))))
                                             {
                                                     out <- NA
                                             } else if (dat.Gene1$UCI < dat.Gene2$LCI)
                                             {
                                                     out <- paste0(Gene2," first")
                                             } else if (dat.Gene2$UCI < dat.Gene1$LCI)
                                             {
                                                     out <- paste0(Gene1," first")
                                             } else
                                             {
                                                     out <- "inconclusive"
                                             }

                                     }

                                     dat.genetics.sample <- dat.genetics[(dat.genetics$ID==samples.Gene1_2[k]),]

                                     # pdf(paste0("./results/supp/Ross/chronology/", Gene1,"_",Gene2, "_",out,"_",samples.Gene1_2[k],".pdf"))
                                     # print(
                                     #       ggplot(dat.genetics.sample) + geom_point(aes(x=Gene, y=VAF, colour=(Gene %in% c(Gene1,Gene2))),size=5) + ylim(0,1) + 
                                     #               geom_segment(aes(x=Gene, xend=Gene, y=LCI, yend=UCI, colour=(Gene %in% c(Gene1,Gene2)))) + theme_bw() +
                                     #                       theme(legend.position="none")
                                     #               #scale_colour_manual(values=col.gene) + theme(legend.position="none")
                                     #               )
                                     # dev.off()

                                     return(out)
                             })

        names(order.info) <- samples.Gene1_2
        return(order.info)


}
