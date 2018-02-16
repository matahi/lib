impute_binary <- function (dat, method=c("zero","one","binomial","majority","svd", "svd-frequency"), lambda=1, seed=42) {


        # Set seed
        set.seed(seed)

        # pre-processing dat
        dat <- as.matrix(dat)

        # Verify that the matrix is binary
        if (!all(dat %in% c(0,1,NA)))
                stop("Matrix is not binary...")

        # Initialize output
        out <- NULL 
        out$dat <- dat

        if (method=="zero")
        {
                out$impute <- dat 
                out$impute[is.na(dat)] <- 0
        }

        if (method=="one")
        {
                out$impute <- dat 
                out$impute[is.na(dat)] <- 1
        }

        # i. Binomial imputation
        if (method=="binomial")
        {
                binomialImpute <- function(x) {x[is.na(x)] <- rbinom(sum(is.na(x)), 1, mean(x, na.rm=TRUE)); return(x)}
                out$impute <- apply(dat, 2, binomialImpute)
        }

        # ii. Majority imputation
        if (method=="majority")
        {
                majorityImpute <- function(x) {x[is.na(x)] <- ifelse(mean(x, na.rm=TRUE)>0.5,1,0); return(x)}
                out$impute <- apply(dat, 2, majorityImpute)
        }

        # iii. Svd imputation
        if (method=="svd")
        {
                require(softImpute)
                lambda0(dat)

                lambda.grid <- c(0,1,5,10,20)

                # Several value of lambda
                dat.soft.impute.complete <- lapply(1:length(lambda.grid), function(n)
                                                   {
                                                           fits <- softImpute(dat, trace=F, rank.max=10, lambda=lambda.grid[n], type="svd")

                                                           # In this case we re-compute everything i.e de-noising + imputation
                                                           dat.soft.impute.complete <- softImpute::complete(dat, fits)
                                                           dat.soft.impute.binary <- dat.soft.impute.complete

                                                           # In this case we only do imputation
                                                           # dat.soft.impute.full <- fits$u %*% (fits$d * diag(length(fits$d))) %*% t(fits$v)
                                                           # dat.soft.impute.binary <- dat.soft.impute.full

                                                           # Should we binarize the output
                                                           dat.soft.impute.binary[dat.soft.impute.binary <=0.5] <- 0
                                                           dat.soft.impute.binary[dat.soft.impute.binary >0.5] <- 1

                                                           colnames(dat.soft.impute.binary) <- colnames(dat)

                                                           return(as.data.frame(dat.soft.impute.binary))
                                                   })

                # Best lambda value

        }

        # iv. Svd imputation with frequency? or is it just that svd imputation is not well done?
        if (method=="svd-frequency")
        {
                require(softImpute)
                lambda0(dat)

                lambda.grid <- c(0,1,5,10,20)

                # Several value of lambda
                dat.soft.impute.complete <- lapply(1:length(lambda.grid), function(n)
                                                   {
                                                           fits <- softImpute(dat, trace=F, rank.max=10, lambda=lambda.grid[n], type="svd")

                                                           # In this case we re-compute everything i.e de-noising + imputation
                                                           dat.soft.impute.complete <- softImpute::complete(dat, fits)
                                                           dat.soft.impute.binary <- dat.soft.impute.complete

                                                           # In this case we only do imputation
                                                           # dat.soft.impute.full <- fits$u %*% (fits$d * diag(length(fits$d))) %*% t(fits$v)
                                                           # dat.soft.impute.binary <- dat.soft.impute.full

                                                           # Should we binarize the output
                                                           dat.soft.impute.binary[dat.soft.impute.binary <=0.5] <- 0
                                                           dat.soft.impute.binary[dat.soft.impute.binary >0.5] <- 1

                                                           colnames(dat.soft.impute.binary) <- colnames(dat)

                                                           return(as.data.frame(dat.soft.impute.binary))
                                                   })

                # Best lambda value

        }

        # Return output
        return(out)
}
