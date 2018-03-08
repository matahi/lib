surv.strat <- function (dat, filter.var, filter.value, analysis.var, analysis.name="", cols.surv = NULL) {

        require(survival)
        require(tidyverse)
        require(rlang)
        require(glue)

        # subset
        dat.current <- dat %>%
                filter((!!sym(filter.var)) %in% filter.value)

        # # variable repartition
        # dat.clinical.analysis %>%
        #         filter((!!sym(filter.var)) %in% filter.value) %>%
        #         count((!!sym(analysis.var)))

        if (all(is.na(dat.current$OS_Status)))
        {
                p.surv <- ggplot() + 
                        geom_point(x=1, y=NA) +
                        ylim(0,1) + xlim(0,5000) +
                        ggtitle(analysis.name) +
                        theme(legend.position="top",
                              plot.title = element_text(hjust=0.5)) +
                        xlab("Time") + ylab("Survival")
                return(p.surv)
        } else
        {
                dat.surv <- survfit(as.formula(glue("Surv(OS, OS_Status) ~ {analysis.var}")), data= dat.current)
        }

        if (is.null(cols.surv))
        {
                cols.final <- brewer.pal(length(unique(dat.current[[analysis.var]])), "Set1")
                names(cols.final) <- unique(dat.current[[analysis.var]])
        } else
        {
                cols.final <- cols.surv[names(cols.surv) %in% as.character(dat.current %>% filter(!is.na(OS_Status)) %>% .[[analysis.var]])]
        }

        p.surv <- ggsurv(dat.surv, CI=F, surv.col=cols.final, order.legend=F) + 
                ylim(0,1) +
                ggtitle(analysis.name) +
                theme(legend.position="top",
                plot.title = element_text(hjust=0.5)) +
                xlim(0,5000)
}


