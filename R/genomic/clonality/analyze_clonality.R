analyze_clonality <- function(dat.binary, dat.genetics, gene)
{

        Gene.samples <- rownames(dat.binary)[which(dat.binary[,gene]==1)]

        # dat.binary <- dat.integrate
        # dat.genetics <- dat.genetics.adj
        # gene <- "FLT3_ITD"


        # dat.genetics.FLT3_ITD <- dat.genetics[dat.genetics$PROTEIN_VARIANT=="FLT3_ITD",]

        # pos.info <- sapply(1:nrow(dat.genetics), function(n)
        #                    {
        #                            PROTEIN_NUM <- strsplit(dat.genetics$PROTEIN_CHANGE[n],"[^0-9]+")[[1]][2]
        #                            as.numeric(PROTEIN_NUM)
        #                    })

        # dat.genetics.FLT3_TKD <- dat.genetics[which((dat.genetics$Gene=="FLT3")&(pos.info > 800)&(pos.info <860)),]

        Gene.clonality <- sapply(1:length(Gene.samples), function(k)
                                 {
                                         # subset the genetics df to the sample
                                         dat.sample <- dat.genetics[(dat.genetics$ID == Gene.samples[k]),]
                                         dat.sample.gene <- dat.genetics[(dat.genetics$ID == Gene.samples[k])&(dat.genetics$variable==gene),]
                                         dat.sample.other <- dat.genetics[(dat.genetics$ID == Gene.samples[k])&(dat.genetics$variable!=gene),]

                                         ###
                                         # Assess clonality
                                         ###

                                         # If no other mutations we assume that the gene is clonal
                                         if (nrow(dat.sample.other)==0)
                                                 return("clonal")

                                         #
                                         Gene.max <- dat.sample.gene[which.max(dat.sample.gene$VAF.corr),] #

                                         comparison <- sapply(1:nrow(dat.sample.other), function(j)
                                                              {
                                                                      m <- round(matrix(c(
                                                                                          Gene.max$VAF.corr*Gene.max$TotRead,
                                                                                          Gene.max$TotRead - Gene.max$VAF.corr*Gene.max$TotRead,
                                                                                          dat.sample.other$VAF.corr[j]*dat.sample.other$TotRead[j],
                                                                                          dat.sample.other$TotRead[j] - dat.sample.other$VAF.corr[j]*dat.sample.other$TotRead[j]),
                                                                                        ncol=2))

                                                                      f2 <- try(fisher.test(m, alternative="less")$p.value < 0.01, silent=T)

                                                                      if (class(f2) !="try-error")
                                                                              if(f2)
                                                                              {
                                                                                      if (Gene.max$VAF.corr + dat.sample.other$VAF.corr[j] > max(dat.sample$VAF.corr,na.rm=T))
                                                                                              return("subclonal")

                                                                                      ## TODO: 1. CAN BE IMPROVED TO DO NESTED PIGEONHOLE
                                                                                      ## TODO: 2. CAN BE IMPROVED TO DO CI instead of hardthreshold


                                                                                      # Some code below:

                                                                                      ## NMut.1 <- Gene.max$VAF.corr*Gene.max$TotRead
                                                                                      ## NMut.2 <- dat.sample.other$VAF.corr[j]*dat.sample.other$TotRead[j]
                                                                                      ## Ntot <- Gene.max$TotRead + dat.sample.other$TotRead[j]

                                                                                      ## # Idx.max <- which.max(dat.sample$VAF.corr)

                                                                                      ## m.clonality <- round(matrix(c(
                                                                                      ##                               NMut.1 + NMut.2,
                                                                                      ##                               Gene.max$TotRead + dat.sample.other$TotRead[j] - (NMut.1 + NMut.2),
                                                                                      ##                               (Gene.max$TotRead + dat.sample.other$TotRead[j])/2, 
                                                                                      ##                               (Gene.max$TotRead + dat.sample.other$TotRead[j])/2,
                                                                                      ##                               ncol=2)))

                                                                                      ## f.clonality <- try(fisher.test(m.clonality, alternative="greater")$p.value < 0.01, silent=T) # Verify pigeonhole principle

                                                                                      ## if (class(f.clonality) !="try-error")
                                                                                      ##         if (f.clonality)
                                                                                      ##                 return("subclonal")


                                                                                      #####################

                                                                                      # Lower CI
                                                                                      # lCI <- qbeta(0.05, NMut.1 + NMut.2 +1 , Ntot - (NMut.1 + NMut.2) + 1)
                                                                                      # if (!is.na(lCI))
                                                                                      #         if (lCI>max(dat.sample$VAF.corr,na.rm=T)) # We suppose that the maxVAF is a clonal mutation
                                                                                      #                 return("subclonal")

                                                                                      # if (!is.na(lCI))
                                                                                      #         if (lCI>0.5) 
                                                                                      #                 return("subclonal")


                                                                              }

                                                                      return("clonal")
                                                              })

                                         if (any(comparison=="subclonal"))
                                                 return("subclonal")

                                         return("clonal")
                                 })
        return(Gene.clonality)
}
