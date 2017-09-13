run_predict_survival_RFX <- function (dat,surv, groups, prop=0.8) 
{
        # dat <- dat.scale
        # surv <- surv.format

        idx.train <- sample(1:nrow(dat), floor(prop*nrow(dat)))


        # Cutting
        dat.train <- dat[idx.train,]
        surv.train <- surv[idx.train,]

        # Training
        cox.RFX <- CoxRFX(dat, surv, groups=groups, which.mu=unique(groups[(groups != "gene_gene")]))
 
        # Testing
        dat.test <- as.data.frame(dat[-idx.train,])
        surv.test <- surv[-idx.train,]
       
        surv.predict <- predict(cox.RFX, dat.test)
        surv.predict.final <- surv.predict[-idx.train]

        # surv.test.formatted <- Surv(surv.test[,1], surv.test[,2])

        return(survConcordance(surv.test ~ surv.predict.final)$concordance)

}



