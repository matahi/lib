run_clonality <- function(dat, Gene_study, dat.survival, cols=colors.gene_study, cols.survival=colors.survival_binary)
{

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
        save(Gene.clonality, file=paste0("results/5_clonality/", Gene_study , "/clonality.RData"))

        ## Plot VAF distribution as a function of clonality
        VAF.gene <- Reduce("rbind",lapply(1:length(Gene.samples), function(k)
                                          {
                                                  dat.df.Gene <- dat[dat$Sample %in% Gene.samples[k],]

                                                  Gene.VAF <- max(dat.df.Gene[dat.df.Gene$Gene==Gene_study,"VAF.corr"])
                                                  Gene.Others <- dat.df.Gene[dat.df.Gene$Gene!=Gene_study,"VAF.corr"]

                                                  return(data.frame(VAF=c(Gene.VAF, Gene.Others), Gene.info=c("Gene", rep("Others", length(Gene.Others))), clonality= Gene.clonality[k] ))
                                          }))

        VAF.gene <- VAF.gene[!is.na(VAF.gene$clonality),]

        pdf(paste0("results/5_clonality/", Gene_study , "/VAF_distribution_clonality.pdf"))
        print(ggplot(VAF.gene) + geom_boxplot(aes(x=clonality, y=VAF, fill=Gene.info), notch=T) + theme_bw() + scale_fill_manual(values=cols) + theme(legend.position="none") )
        dev.off()

        pval <- t.test(VAF.gene$VAF[(VAF.gene$clonality=="subclonal")&(VAF.gene$Gene.info=="Gene")], VAF.gene$VAF[VAF.gene$clonality=="clonal"&(VAF.gene$Gene.info=="Gene")])

        ## Plot VAF distribution as a function of clonality
        max.VAF <- Reduce("rbind",lapply(1:length(Gene.samples), function(k)
                                         {
                                                 dat.df.Gene <- dat[dat$Sample %in% Gene.samples[k],]

                                                 max.VAF <- max(dat.df.Gene[,"VAF.corr"],na.rm=T)

                                                 return(data.frame(VAF=max.VAF , clonality= Gene.clonality[k] ))
                                         }))

        max.VAF <- max.VAF[!is.na(max.VAF$clonality),]

        pdf(paste0("results/5_clonality/", Gene_study , "/max_VAF_vs_clonality.pdf"))
        print(ggplot(max.VAF) + geom_boxplot(aes(x=clonality, y=VAF)) + theme_bw() + scale_fill_manual(values=cols) + theme(legend.position="none") )
        dev.off()

        # survival clonal vs subclonal
        # EFS
        dat.survival.gene <- dat.survival[which(dat.survival[,Gene_study]==1),]
        dat.survival.gene <- cbind(dat.survival.gene, clonality=Gene.clonality)

        dat.survival.gene$clonality <- factor(as.character(dat.survival.gene$clonality), levels=c("subclonal","clonal"))
        

        source("src/lib/run_univariate_survival.R")

        cox.fit.OS <- univariate_survival_clonality(dat.survival=dat.survival.gene, surv.type="OS", Gene=Gene_study, stat=F, cols=cols.survival)

        cox.fit.EFS <- univariate_survival_clonality(dat.survival=dat.survival.gene, surv.type="EFS", Gene=Gene_study, stat=F, cols=cols.survival)

        return(list(clonality=Gene.clonality, pval=pval, pval.EFS=summary(cox.fit.EFS)$coefficients[5], pval.OS=summary(cox.fit.OS)$coefficients[5]))
}
