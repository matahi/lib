analyze_clonal_state <- function(dat.binary, dat.genetics, gene) {

        Gene.samples <- rownames(dat.binary)[which(dat.binary[,gene]==1)]

        Gene.clonality <- sapply(1:length(Gene.samples), function(k)
                                 {
                                         # subset the genetics df to the sample
                                         dat.sample <- dat.genetics[(dat.genetics$ID == Gene.samples[k]),]

                                         ###
                                         # Assess clonality
                                         ###

                                         # If only one mutation then clonal
                                         if (nrow(dat.sample)==1)
                                                 return("unique clone")

                                         #
                                         # Gene.max <- dat.sample.gene[which.max(dat.sample.gene$VAF),]

                                         pair.indexes <- t(combn(1:nrow(dat.sample),2))

                                         pair.comparison <- sapply(1:nrow(pair.indexes), function(j)
                                                              {
                                                                      idx.1 <- pair.indexes[j,1]
                                                                      idx.2 <- pair.indexes[j,2]

                                                                      m <- round(matrix(c(
                                                                                          dat.sample$VAF.corr[idx.1]*dat.sample$TotRead[idx.1],
                                                                                          dat.sample$TotRead[idx.1] - dat.sample$VAF.corr[idx.1]*dat.sample$TotRead[idx.1],
                                                                                          dat.sample$VAF.corr[idx.2]*dat.sample$TotRead[idx.2],
                                                                                          dat.sample$TotRead[idx.2] - dat.sample$VAF.corr[idx.2]*dat.sample$TotRead[idx.2]),
                                                                                        ncol=2))

                                                                      f1 <- try(fisher.test(m, alternative="less")$p.value < 0.01, silent=T)
                                                                      f2 <- try(fisher.test(m, alternative="greater")$p.value < 0.01, silent=T)

                                                                      if (class(f1) !="try-error")
                                                                              if((f1) & dat.sample$VAF.corr[idx.1] + dat.sample$VAF.corr[idx.2] > max(dat.sample$VAF.corr,na.rm=T))
                                                                                              return("clonal heterogeneity")

                                                                      if (class(f2) !="try-error")
                                                                              if((f2) & dat.sample$VAF.corr[idx.1] + dat.sample$VAF.corr[idx.2] > max(dat.sample$VAF.corr,na.rm=T))
                                                                                              return("clonal heterogeneity")

                                                                      return("unique clone")
                                                              })

                                         if (any(pair.comparison=="clonal heterogeneity"))
                                                 return("clonal heterogeneity")

                                         return("unique clone")
                                 })
        return(Gene.clonality)
}
