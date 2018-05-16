analyze_chronology <- function(dat.binary, dat.genetics, Gene1, Gene2)
{
        require(tidyverse)
        # Samples 
        samples.Gene1_2 <- dat.binary %>% filter(dat.binary[[Gene1]]==1,dat.binary[[Gene2]]==1) %>% .$ID

        if (length(samples.Gene1_2)==0)
                return(NULL)

        # What happens if 2 mutations?
        order.info <- sapply(1:length(samples.Gene1_2), function(k)
                             {
                                     dat.Gene1 <- dat.genetics %>% filter(ID == samples.Gene1_2[k], variable == Gene1) 
                                     dat.Gene2 <- dat.genetics %>% filter(ID == samples.Gene1_2[k], variable == Gene2) 

                                     ######
                                     if (nrow(dat.Gene1)==0|nrow(dat.Gene2)==0)
                                     {
                                             print("Careful binary matrix and VCF files do not correspond")
                                             print(Gene1)
                                             print(Gene2)
                                             return("inconclusive")
                                     } else 
                                     {
                                             dat.sample <- dat.genetics %>% filter(ID == samples.Gene1_2[k])

                                             Gene1.max <- dat.Gene1[which.max(dat.Gene1$VAF),]
                                             Gene2.max <- dat.Gene2[which.max(dat.Gene2$VAF),]

                                             m <- round(matrix(c(
                                                                 Gene1.max$VAF.corr*Gene1.max$TotRead,
                                                                 Gene1.max$TotRead - Gene1.max$VAF.corr*Gene1.max$TotRead,
                                                                 Gene2.max$VAF.corr*Gene2.max$TotRead,
                                                                 Gene2.max$TotRead - Gene2.max$VAF.corr*Gene2.max$TotRead),
                                                               ncol=2))

                                             f1 <- try(fisher.test(m, alternative="greater")$p.value < 0.01, silent=T)
                                             if (class(f1) !="try-error")
                                             {
                                                     if (f1 & (Gene1.max$VAF.corr + Gene2.max$VAF.corr > max(dat.sample$VAF.corr,na.rm=T)))
                                                     {
                                                             return(paste0(Gene1, " first"))
                                                     } 
                                             }

                                             f2 <- try(fisher.test(m, alternative="less")$p.value < 0.01, silent=T)
                                             if (class(f2) !="try-error")
                                             {
                                                     if (f2 & (Gene1.max$VAF.corr + Gene2.max$VAF.corr > max(dat.sample$VAF.corr,na.rm=T)))
                                                     {
                                                             return(paste0(Gene2, " first"))
                                                     } 
                                             }

                                             return("inconclusive")
                                     } 
                             })

        names(order.info) <- samples.Gene1_2
        return(order.info)

}





#         for(s in clinicalData$PDID)
#         {
#                 l <- list()
#                 for(i in which(mutationData$SAMPLE_NAME==s & ix))
#                 {
#                         for(j in which(mutationData$SAMPLE_NAME==s & ix))
#                         { 
#                                 if(!is.na(cn[i]) & !is.na(cn[j]) & i!=j){  
#                                         m <- round(matrix(c(
#                                                             mcf[i]*depth[i],
#                                                             depth[i]-mcf[i]*depth[i],
#                                                             mcf[j]*depth[j],
#                                                             depth[j]-mcf[j]*depth[j]),
#                                                           ncol=2))
#                                         f <- try(fisher.test(m, alternative="greater")$p.value < 0.01, silent=T)
#                                         if (class(f) !="try-error")
#                                         {
#                                                 if (f & mcf[i] >= 1 - mcf[j]) ## TODO ADVANCED PIDGEONHOLE
#                                                 {
# 
#                                                 }
#                                         }
#                                 }
#                         }
#                 }
#         }


