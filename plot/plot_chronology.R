plot_chronology <- function(mutation.df, gene, out.dir="results/", col.gene)
{
        require(ggplot2)

        gene.samplelist <- unique(as.character(mutation.df[mutation.df$Gene==gene,"ID"]))

        if (length(gene.samplelist)==0)
        {
                return(NA)
        } else
        {
                chronology.list <- sapply(1:length(gene.samplelist), function(n)
                              {
                                      gene.sample <- gene.samplelist[n]
                                      

                                      dat.Gene <- mutation.df[(mutation.df$ID==gene.sample)&(mutation.df$Gene==gene),,drop=F]
                                      dat.Others <- mutation.df[(mutation.df$ID==gene.sample)&(mutation.df$Gene!=gene),,drop=F]

                                      if (nrow(dat.Others)==0)
                                      {
                                              chronology <- paste0(gene,"_alone")

                                              
                                      } else if (nrow(dat.Gene)>1|nrow(dat.Others)>1)
                                      {
                                              dat.sample <- mutation.df[mutation.df$ID==gene.sample,]

                                              chronology <- "inconclusive"

                                              for (i in 1:nrow(dat.Gene))
                                              {
                                                      if (all(dat.Gene$LCI[i]>dat.Others$UCI,na.rm=T))
                                                      {
                                                              chronology <- paste0(gene,"_first")
                                                      } 
                                              }

                                              for (i in 1:nrow(dat.Others))
                                              {
                                                      if (all(dat.Others$LCI[i]>dat.Gene$UCI,na.rm=T))
                                                      {
                                                              chronology <- paste0(gene,"_late")
                                                      } 
                                              }

                                      } else
                                      {
                                              if (any(is.na(c(dat.Gene$LCI, dat.Gene$UCI, 
                                                              dat.Others$LCI, dat.Others$UCI))))
                                              {
                                                      chronology <- NA
                                              } else if (dat.Gene$UCI < dat.Others$LCI)
                                              {
                                                      chronology <- paste0(gene,"_late")
                                              } else if (dat.Others$UCI < dat.Gene$LCI)
                                              {
                                                      chronology <- paste0(gene,"_first")
                                              } else
                                              {
                                                      chronology <- "inconclusive"
                                              }

                                      }

                                      mutation.df.sample <- mutation.df[mutation.df$ID==gene.sample,,drop=F]
                                      mutation.df.sample$Gene.info <- sapply(1:nrow(mutation.df.sample), function(n)
                                                                             {
                                                                                     if (mutation.df.sample$Gene[n]==gene)
                                                                                     {
                                                                                             return("Gene")
                                                                                     } else
                                                                                     {
                                                                                             return("Others")
                                                                                     }
                                                                             })

                                      pdf(paste0(out.dir, "/", gene.sample,"_",chronology,".pdf"))
                                      print(
                                            ggplot(mutation.df.sample) + geom_point(aes(x=Gene, y=VAF, colour=Gene.info),size=5) + ylim(0,1) + 
                                            geom_segment(aes(x=Gene, xend=Gene, y=LCI, yend=UCI, colour=Gene.info)) + theme_bw() +
                                            scale_colour_manual(values=col.gene) + theme(legend.position="none")
                                            )
                                      dev.off()
                                      
                                      return(chronology)


                              })
                names(chronology.list) <-gene.samplelist

                return(chronology.list)
        }
}

