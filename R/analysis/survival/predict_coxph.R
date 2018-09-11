run_predict_coxph <- function (dat,surv, prop=0.8) 
{
        # dat <- dat.filter
        # surv <- surv.filter


        idx.train <- sample(1:nrow(dat), floor(prop*nrow(dat)))

        # Cutting
        dat.train <- dat[idx.train,,drop=F]
        surv.train <- surv[idx.train,,drop=F]

        # Training
        dat.training <- data.frame(dat.train, surv=surv.train) 
        cox.fit <- coxph(Surv(surv.time, surv.status) ~ ., dat.training)
 
        # Testing
        dat.test <- data.frame(dat[-idx.train,,drop=F])
        surv.test <- surv[-idx.train,]
       
        surv.predict <- predict(cox.fit, dat.test)

        surv.test.formatted <- Surv(surv.test[,1], surv.test[,2])

        return(survConcordance(surv.test.formatted ~ surv.predict)$concordance)
}



