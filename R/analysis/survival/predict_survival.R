predict_survival <- function (dat, surv, lambda.value= "1se") 
{
        # Status?
        library(glmnet)
        tmp <- cv.glmnet(as.matrix(dat), as.matrix(surv), family="cox")

        # beta.mat <- tmp$beta
        lambda.model <- tmp$lambda.1se

        if (lambda.value == "min")
                lambda.model <- tmp$lambda.min

        # or take the lambda optimal?

        tmp.bis <- glmnet(as.matrix(dat), as.matrix(surv), lambda=lambda.model, family="cox")

        return(tmp.bis)

}

run_predict_survival <- function (dat,surv, prop=0.8, lambda.value = "1se") 
{
        idx.train <- sample(1:nrow(dat), floor(prop*nrow(dat)))

        # Cutting
        dat.train <- dat[idx.train,]
        surv.train <- surv[idx.train,]

        # Training
        cox.glm <- predict_survival(dat.train, surv.train, lambda.value = lambda.value)
 
        # Testing
        dat.test <- dat[-idx.train,]
        surv.test <- surv[-idx.train,]
       
        surv.predict <- predict(cox.glm, as.matrix(dat.test))

        surv.test.formatted <- Surv(surv.test[,1], surv.test[,2])

        return(survConcordance(surv.test.formatted ~ surv.predict)$concordance)

}



