run_predict_survival_RFX <- function (dat,surv, groups, prop=0.8) 
{
        idx.train <- sample(1:nrow(dat), floor(prop*nrow(dat)))

        # Cutting
        dat.train <- dat[idx.train,]
        surv.train <- surv[idx.train]

        # Training
        cox.RFX <- CoxRFX(dat.train, surv.train, groups=groups, which.mu=setdiff(names(groups), "Genes_Genes"))
 
        # Testing
        dat.test <- as.data.frame(dat[-idx.train,])
        surv.test <- surv[-idx.train]
       
        surv.predict <- predict(cox.RFX, dat.test)

        return(survConcordance(surv.test ~ surv.predict)$concordance)

}



