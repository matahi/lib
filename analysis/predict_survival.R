predict_survival <- function (dat, surv) 
{
        # Status?
        library(glmnet)
        tmp <- cv.glmnet(as.matrix(dat), as.matrix(surv), family="cox")

        # beta.mat <- tmp$beta
        lambda.1se <- tmp$lambda.1se

        tmp.bis <- glmnet(as.matrix(dat), as.matrix(surv), lambda=lambda.1se, family="cox")

        return(tmp.bis)

}

run_predict_survival <- function (dat,surv, prop=0.8) 
{
        idx.train <- sample(1:nrow(dat), floor(prop*nrow(dat)))

        # Cutting
        dat.train <- dat[idx.train,]
        surv.train <- surv[idx.train,]

        # Training
        cox.glm <- predict_survival(dat.train, surv.train)
 
        # Testing
        dat.test <- dat[-idx.train,]
        surv.test <- surv[-idx.train,]
       
        surv.predict <- predict(cox.glm, as.matrix(dat.test))

        surv.test.formatted <- Surv(surv.test[,1], surv.test[,2])

        return(survConcordance(surv.test.formatted ~ surv.predict)$concordance)

}



