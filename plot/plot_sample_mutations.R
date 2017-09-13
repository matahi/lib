plot_sample_mutations <- function () {
        
        Gene.samples <- unique(dat[dat$Gene==Gene_study,"Sample"])

        Gene.clonality <- sapply(1:length(Gene.samples), function(k)
                                 {
                                         if (sum(dat$Sample %in% Gene.samples[k])==1)
                                         {
                                                 clonality <- "clonal"
                                         } else 
                                         {
                                                 dat.df.Gene <- dat[dat$Sample %in% Gene.samples[k],]

                                                 Gene.UCI <- max(dat.df.Gene[dat.df.Gene$Gene==Gene_study,"UCI.corr"])
                                                 Others.LCI <- na.omit(dat.df.Gene[dat.df.Gene$Gene!=Gene_study,"LCI.corr"])

                                                 if (is.na(Gene.UCI))
                                                 {
                                                         return(NA)
                                                 } else if (any(Gene.UCI < Others.LCI))
                                                 {
                                                         clonality <- "subclonal"
                                                 } else
                                                 {
                                                         clonality <- "clonal"
                                                 }

                                                 pdf(paste0("results/5_clonality/", Gene_study , "/",clonality,"_", Gene.samples[k],".pdf"))
                                                 print(
                                                       ggplot(dat.df.Gene) + geom_point(aes(x=Gene, y=VAF.corr, colour=(Gene==Gene_study)),size=5) + ylim(0,1) + 
                                                       geom_segment(aes(x=Gene, xend=Gene, y=LCI.corr, yend=UCI.corr, colour=(Gene==Gene_study))) + theme_bw() +
                                                       scale_colour_manual(values=cols) + theme(legend.position="none")
                                                       )
                                                 dev.off()

                                         }
                                         return(clonality)
                                 })
        names(Gene.clonality) <- Gene.samples

}
