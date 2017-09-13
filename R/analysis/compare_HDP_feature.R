compare_HDP_feature <- function (HDP.df, out.merge, feature.compare, out.dir="") {

        features.comparison <- unique(HDP.df[,feature.compare])

        for (class.num in levels(HDP.df$HDP))
        {
                dat.list <- lapply(1:length(features.comparison), function(k)
                                   {
                                           return(data.frame(out.merge$dat[ HDP.df$HDP== class.num & HDP.df[,feature.compare]==features.comparison[k],, drop=F]))
                                   })

                names(dat.list) <- features.comparison

                if (all(sapply(dat.list, nrow) > 5))
                {
                        source("src/lib/plot/plot_compare_groups.R")
                        comparison <- plot_compare_groups(dat.list, add.significance=T)

                        # cols.clusters <- brewer.pal(9,"Set3")[c(3,5)]

                        pdf(paste0(out.dir, "/HDP_", class.num, "_", feature.compare,".pdf"), width=16, height=9)
                        print(comparison$p + scale_fill_manual(values=colors.myeloid) +  theme(text= element_text(size=20)))
                        dev.off()
                }
        }


        return(NULL)

}
