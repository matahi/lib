chisq <- function(n_occ, n_tot,  pop.occ, pop.tot)
{

        # n <- 1726
        # n_occ <- ID_VAR.df$n_occ[n]
        # n_tot <-  ID_VAR.df$n_tot[n]
        # pop.occ <- ID_VAR.df$pop.occ[n]
        # pop.tot <- ID_VAR.df$pop.tot[n]

        if (any(is.na(c(n_occ, n_tot,pop.occ, pop.tot))))
        {
                return(NA)
        } else
        {
                return(chisq.test(rbind(c(n_occ,n_tot-n_occ),
                                        c(pop.occ, pop.tot-pop.occ)))$p.value)
        }
}
