impute_continuous <- function (dat, method=c("mean"), lambda=1, seed=42) {

        # Set seed
        set.seed(seed)

        # pre-processing dat
        dat <- as.matrix(dat)

        # Initialize output
        out <- NULL 
        out$dat <- dat

        if (method=="mean")
        {
                out$impute <- dat 
                out$impute[is.na(dat)] <- mean(dat, na.rm=T)
        }

        if (method=="gaussian")
        {
                out$impute <- dat 
                mean_info <- mean(dat, na.rm=T)
                sd_info <- sd(dat, na.rm=T)
                out$impute[is.na(dat)] <- rnorm(sum(is.na(dat)), mean_info, sd_info)
        }

        if (method=="beta")
        {
                out$impute <- dat 
                out$impute[is.na(dat)] <- mean(dat, na.rm=T)
        }

        return(out)
}
