glm_clinical_mutation <- function(dat.survival, non_zeros, clinical.feature, n_rep=100, prop=0.8)
{
        
        require(glmnet)
        require(lars)

        # clinical.feature <- clinical.list[n]

        # dat.survival <- dat.survival[-which(dat.survival$Age==0),]

        #n_rep <- 10

        cor.list <- lapply(1:n_rep, function(n)
                           {
                                   # print(n)
                                   #######
                                   n_size <- nrow(dat.survival)

                                   # Binary classification make sure to have samples from both classes in train and test
                                   if (length(unique(dat.survival[,clinical.feature]))==2)
                                   {
                                           classes <- unique(dat.survival[,clinical.feature])

                                           idx.1 <- which(dat.survival[,clinical.feature]==classes[1])
                                           idx.2 <- which(dat.survival[,clinical.feature]==classes[2])

                                           idx.train.1 <- sample(n_size, floor(prop * length(idx.1)))
                                           idx.train.2 <- sample(n_size, floor(prop * length(idx.2)))

                                           dat.train <- dat.survival[c(idx.train.1,idx.train.2),]
                                           dat.test <- dat.survival[-c(idx.train.1,idx.train.2),]

                                   } else # Regression
                                   {

                                           n.train <- floor(prop * n_size)
                                           idx.train <- sample(n_size,n.train)

                                           dat.train <- dat.survival[idx.train,]
                                           dat.test <- dat.survival[-idx.train,]

                                   }

                                   # filter genes not present in train
                                   non_zeros.table <- sapply(1:length(non_zeros), function(k)
                                                             {
                                                                     length(unique(dat.train[,non_zeros[k]]))
                                                             })

                                   non_zeros.filtered <- non_zeros[non_zeros.table!=1]

                                   #
                                   gene_list <- paste(paste0("factor(",non_zeros.filtered,")"),collapse="+")
                                   f <- paste0(clinical.feature," ~", gene_list)

                                   # Do a glm
                                   fit2 <- glm(as.formula(f)  , data=dat.train)
                                   var.predict <- predict(fit2, dat.test, type="response")

                                   R_test <- cor(var.predict, dat.test[,clinical.feature], use="pairwise.complete.obs")

                                   return(R_test)

                                   ## # Do a stability selection 1000 times then rank the variables
                                   ## # fit.lars <- lars(dat.train[,non_zeros.filtered], dat.train[,clinical.feature])
                                   ## x <- as.matrix(dat.train[,non_zeros.filtered])
                                   ## y <- dat.train[,clinical.feature]

                                   ## fit.lars <- lars(x, y)
                                   ## features.rank <- order(fit.lars$entry)

                                   ## #
                                   ## x.test <-  as.matrix(dat.test[,non_zeros.filtered])
                                   ## y.test <-  dat.test[,clinical.feature]

                                   ## #
                                   ## tmp <- predict(fit.lars, x.test)

                                   ## MSE <- sapply(1:ncol(tmp$fit), function(n)
                                   ##               {
                                   ##                       sum((tmp$fit[,n] - y.test)*(tmp$fit[,n] - y.test))
                                   ##               })

                                   ## k.opt <- which.min(MSE)

                                   #
                                   # return(list(R=R_test,features.rank= features.rank, k.opt=k.opt))
                           })

        # feat.rank <- Reduce("rbind", lapply(1:length(cor.list), function(n)
        #                                     {
        #                                             matrix(cor.list[[n]]$features.rank,nrow=1)
        #                                     }))
        # feat.rank.df <- data.frame(feat.rank)
        # colnames(feat.rank.df) <- non_zeros.filtered

        # k.opt <- sapply(1:length(cor.list), function(n){cor.list[[n]]$k.opt})
        # # R <- sapply(1:length(cor.list), function(n){cor.list[[n]]$R})

        # #k.final <- floor(mean(k.opt))
        # k.final <- 10


        # ############

        # first.k <- lapply(1:nrow(feat.rank), function(n)
        #                   {
        #                           return(which(feat.rank[n,] <= k.final))

        #                   })

        # tmp <- factor(unlist(first.k),levels=1:ncol(feat.rank))
        # 
        # toto <- table(tmp)

        # non_zeros.filtered[order(toto,decreasing=T)]

        # cor(dat.survival$Age, dat.survival$EZH2)

        # cor(dat.survival$Age, dat.survival$NPM1)

        # cor.info <- sapply(1:length(non_zeros.filtered), function(n)
        #                    {
        #                            cor(dat.survival$Age, dat.survival[,non_zeros.filtered[n]])
        #                    })

        # pval <- sapply(1:length(non_zeros.filtered), function(n)
        #                    {
        #                            t.test(dat.survival$Age[dat.survival[,non_zeros.filtered[n]]==1], dat.survival$Age[dat.survival[,non_zeros.filtered[n]]==0])$p.value
        #                    })
        # p.adj <- p.adjust(pval, method="fdr")

        # variables <- which(p.adj < 0.05)

        # 

        # print(ggplot(dat.survival.impute) + geom_boxplot(aes(x=factor(NPM1),y=Age)) + theme_bw())
        # print(ggplot(dat.survival.impute) + geom_boxplot(aes(x=factor(FLT3_ITD),y=Age)) + theme_bw())

        # # library(reshape2)

        # # feat.m <- melt(feat.rank.df)
        # # feat.m$variable <- factor(feat.m$variable, levels=non_zeros.filtered)

        # # ggplot(feat.m) + geom_boxplot(aes(x=variable, y=value)) + theme(axis.text.x=element_text(angle=90, hjust=1))

        # # avg.rank <- apply(feat.rank,2,sum)

        # # avg.df <- data.frame(feature=factor(non_zeros.filtered,levels=non_zeros.filtered), rank=avg.rank)

        # # quartz()
        # # ggplot(avg.df) + geom_boxplot(aes(x=feature, y=rank)) + theme(axis.text.x=element_text(angle=90, hjust=1))

        return(cor.list)

        # #####
        # qplot(var.predict, dat.test$Age)
}

run_glm_clinical_mutation <- function(dat.survival, non_zeros, clinical.list, n_rep=1000, prop=0.8)
{
        # #
        # dat.survival <- dat.survival.impute
        # non_zeros <- sparse.OS$non_zeros
        # n_rep <- 1000
        # n <- 2

        # dat.survival <- dat.survival.impute
        # non_zeros <- sparse.OS$non_zeros
        # prop <- 0.8


        # prop <- .8


        #

        cor.df <- Reduce("rbind", lapply(1:length(clinical.list), function(n)
                                         {
                                                 #print(n)
                                                 cor.out <- glm_clinical_mutation(dat.survival, non_zeros, clinical.list[n])

                                                 return(data.frame(R=unlist(cor.out), feature=clinical.list[n]))
                                         }))

        # Plot the correlation
        require(ggplot2)

        p.R <- ggplot(cor.df) + geom_boxplot(aes(x=feature, y=R))

        # Calculate the correlation with each clinical variables
        # p.corr <- ggplot(cor.df) + geom_boxplot(aes(x=feature, y=R))


        return(p.R)
}


# ##### Beta_2
# beta_2 <- dat$beta_2
# 
# gaussianImpute <- function(x) 
# {
#         mean_var <- mean(x,na.rm=T)
#         sd_var <- sd(x,na.rm=T)
#         x[is.na(x)] <- rnorm(sum(is.na(x)),  mean_var, sd_var); 
#         return(x)
# }
# 
# set.seed(42)
# beta_2 <- gaussianImpute(beta_2)


