## First lets run lasso

sparse_multivariate_survival <- function(dat.survival, surv.type, feature_list, cols, out.dir="results/")
{
        # dat.survival <- dat.cluster
        # surv.type <- "OS"
        # feature_list <- gene.list.cluster
        # cols <- colors.survival_binary
        # out.dir <- paste0("results/datasets/",dataset,"/",levels(clusters)[n])



        require(survival)
        require(glmnet)

        final.dir <- paste0(out.dir, "/survival/", surv.type, "/sparse/")
        system(paste("mkdir -p", final.dir))

        # dat.survival <- dat.survival.impute
        # surv.type <- "OS"
        # feature_list <- gene_frequent
        # cols <- colors.survival_binary

        if (surv.type=="OS")
        {
                ylab <- "Overall survival"
        } else if (surv.type=="EFS")
        {
                ylab <- "Disease-free survival"
        }

        # We have to take the imputed matrix?
        if (sum(is.na(dat.survival$OS))>0)
                dat.survival <- dat.survival[ -which(is.na(dat.survival$OS)),]
        if (sum(dat.survival$OS==0)>0)
                dat.survival <- dat.survival[ -which(dat.survival$OS==0),]

        # Remove samples with NAs now
        feature.filter <- union(which(colSums(dat.survival[, feature_list])==0), which(colSums(dat.survival[, feature_list])==nrow(dat.survival)))
        if (length(feature.filter)>0)
                feature_list <- feature_list[-feature.filter]



        #


        x <- as.matrix(dat.survival[,feature_list])
        y <- as.numeric(dat.survival[ ,surv.type])
        status <- as.numeric(dat.survival[,paste0(surv.type,"_Status")])

        cv.fit <- cv.glmnet(x, Surv(y,status), family = "cox", alpha = 1)
        fit.final <- glmnet(x, Surv(y,status), family = "cox", alpha = 1, lambda=cv.fit$lambda.min)

        non_zeros <- feature_list[which(fit.final$beta!=0)]

        if (length(non_zeros)==0)
                return(NULL)

        #gene_list <- paste(paste0("factor(",non_zeros,")"),collapse="+")
        gene_list <- paste(paste0("factor(",feature_list,")"),collapse="+")

        f <- paste0("Surv(",surv.type,",",surv.type,"_Status) ~", gene_list)
        cox.fit <- coxph(as.formula(f), data= dat.survival)

        lasso.coef <- sapply(1:length(summary(cox.fit)$coefficients[,5]), function(n)
                             {
                                     if (summary(cox.fit)$coefficients[n,5]>0.05)
                                     {
                                             return("NS")
                                     } else if (summary(cox.fit)$coefficients[n,5]<0.0001)
                                     {
                                             return("<.0001")
                                     } else if (summary(cox.fit)$coefficients[n,5]<0.05)
                                     {
                                             return("<.05")
                                     }
                             })

        surv.info <- cbind(round(summary(cox.fit)$conf.int[,c(1,3:4)],2),lasso.coef)

        write.table(surv.info, file=paste0(final.dir,"/features.txt"), col.names=T, row.names=T, quote=F, sep="\t")

        # i) do forest plot
        # require(forestplot)
        # pdf(paste0("results/3_survival/",surv.type,"_forestplot.pdf"))
        # forestplot(non_zeros , summary(cox.fit)$conf.int[,1],  summary(cox.fit)$conf.int[,3], summary(cox.fit)$conf.int[,4])
        # dev.off()

        # ii) plot univariate of the non_zero variables 
        for (i in 1:length(non_zeros))
        {
                uni_feature <- non_zeros[i]
                form <- paste0("Surv(",surv.type,",", surv.type,"_Status) ~ factor(", uni_feature,")")

                KM.fit <- survfit(as.formula(form), data= dat.survival)
                cox.fit <- coxph(as.formula(form), data= dat.survival)

                legend.names <- paste0(uni_feature,"=",0:1)

                pdf(paste0(final.dir,"/",uni_feature ,".pdf"))
                plot(KM.fit, col=cols, mark=NA, lwd =2, ,xlab="days", ylab =ylab)
                legend("topright", bty="", legend.names, col=cols,cex=0.7, lwd=2, lty=1)
                dev.off()

                # write.table(summary(cox.fit)$coefficients, file=paste0("results/3_survival/lasso_",uni_feature ,".txt"),col.names=T, row.names=T, quote=F, sep="\t")
        }


        # return surv.info
        return(list(surv.info=surv.info,non_zeros=non_zeros))

}

