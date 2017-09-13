analyze_HDP_MDS <- function(out.HDP, out.merge, out.clinical=NULL, grouping=NULL, out.dir="", display_legend=F)
{
        # n.HDP <- length(unique(out.HDP$class))
        # out.dir <- "clustering/ALL"
        # grouping <- grouping_AML_MDS
        # display_legend <- F

        # out.clinical <- dat.survival.impute

        # grouping <- NULL
        # out.merge <- dat.integrate.impute.bis 
        # out.clinical <- dat.clinical
        # out.HDP 

        n.HDP <- length(unique(out.HDP$class))

        ## Lets process HDP

        prob <- apply(out.HDP$prob,1,max)

        library(knitr)

        ### Use the clusters.processed

        # Survival by HDP class
        HDP.df <- data.frame(HDP=out.HDP$class, out.clinical)

        # Survival
        source("src/lib/analysis/run_univariate_survival.R")
        pp <- univariate_survival_HDP(HDP.df, surv.type="OS", cols=cols, stats=F, out.dir=out.dir, cluster.names=NULL)

        # pp <- univariate_survival_HDP(HDP.df, surv.type="EFS", cols=cols, stats=T, out.dir=out.dir)

        # Other clinical variables
        pdf(paste0('results/', out.dir,'/HDP_Age.pdf'))
        print(ggplot(HDP.df) + geom_boxplot(aes(x=HDP,y=Age, fill=HDP)) + scale_fill_manual(values=cols) + theme_bw()
                             + theme(axis.text.x=element_text(angle=90)))
        dev.off()

        table(clusters.processed, useNA="always")

}
