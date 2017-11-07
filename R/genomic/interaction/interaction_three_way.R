interaction_three_way <- function (n_tot, n1, n2, n3, observed.val=NULL, bootstraps=10000) {

        # n_tot <- 2330
        # n1 <- 320
        # n2 <- 220 
        # n3 <- 300 
        # observed.val <- 30

        three_way.bootstrap <- sapply(1:bootstraps, function(n)
                                      {
                                              vect1 <- rbinom( n_tot, 1 , prob= n1/n_tot)
                                              vect2 <- rbinom( n_tot, 1 , prob= n2/n_tot)
                                              vect3 <- rbinom( n_tot, 1 , prob= n3/n_tot)

                                              return(sum((vect1==1)&(vect2==1)&(vect3==1)))
                                      })

        ###
        CI.95 <- quantile(three_way.bootstrap, c(0.025, 0.975))
        hist(three_way.bootstrap)
        abline(v=CI.95[1], col="blue") # 95% lower CI
        abline(v=CI.95[2], col="blue") # 95% higher CI

        if (!is.null(observed.val))
        {
                abline(v=observed.val, col="red")
        }

        if is.null(observed.val)
        {
                return(list(test= NULL     ,   infos=list(mean=mean(three_way.bootstrap), median=median(three_way.bootstrap), low95= CI.95[1], upp95=CI.95[2])))
        } else
        {
                test <- (observed.val > CI.95[2])|(observed.val < CI.95[1])
                return(list(test= test     ,   infos=list(mean=mean(three_way.bootstrap), median=median(three_way.bootstrap), low95= CI.95[1], upp95=CI.95[2])))
        }


        
}
