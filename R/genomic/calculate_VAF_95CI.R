calculate_VAF_95CI <- function(Nmut, Ntot)
{
        ## calculate the 95 CI for the VAF given the number of mutated reads Nmut and the total number of reads Ntot
        return(qbeta(c(0.025,0.975), Nmut +1, Ntot - Nmut + 1))
}
