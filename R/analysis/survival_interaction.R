survival_interaction <- function(dat.survival, surv.type, feature1, feature2, cols=c("#abd9e9","#2c7bb6","#fdae61","#d7191c"), plot=T, out.dir="results/")
{
        require(survival)
        # #
        # dat.survival <- dat.survival.TP53
        # surv.type <- "OS"
        # feature1 <- "TP53"
        # feature2 <- "del"
        # cols <- colors.survival_interaction


        # dat.survival <- dat.final
        # surv.type <- "OS"
        # feature1 <- "TP53"
        # feature2 <- "Complex"


        #

        final.dir <- paste0(out.dir, "/survival/",surv.type,"/interaction/")
        system(paste("mkdir -p", final.dir))

        if (surv.type=="OS")
        {
                ylab <- "Overall survival"
        } else if (surv.type=="EFS")
        {
                ylab <- "Disease-free survival"
        }
        form <- paste0("Surv(", surv.type,", ", surv.type,"_Status) ~ factor(", feature1,") + factor(", feature2, ")" )

        KM.fit <- survfit(as.formula(form), data= dat.survival)
        cox.fit <- coxph( as.formula(gsub("\\+","*",form)), data= dat.survival)

        legend.names <- paste0(feature1,"=",c(0,0,1,1),", ", feature2,"=", c(0,1,0,1))

        if (plot)
        {
                pdf(paste0(final.dir, feature1, "_", feature2, ".pdf"))
                plot(KM.fit, col=cols, mark=NA, lwd =2, ,xlab="days", ylab=ylab)
                legend("topright", bty="", legend.names, col=cols,cex=0.7, lwd=2, lty=1)
                dev.off()
        }

        pval.interaction <- summary(cox.fit)$coefficients[3,5]

        out <- list(pval.interaction=pval.interaction, cox.fit=cox.fit, table=table(dat.survival[,feature1], dat.survival[,feature2]))

}


survival_interaction_treatment <- function(dat.survival, surv.type, feature, cols, plot=T)
{
        n_treatments <- length(unique(dat.survival$Treatment))

        if (surv.type=="OS")
        {
                ylab <- "Overall survival"
        } else if (surv.type=="EFS")
        {
                ylab <- "Disease-free survival"
        }
        form <- paste0("Surv(", surv.type,"_followUp, ", surv.type,"_status) ~ factor(", feature,") + factor(Treatment)" )

        KM.fit <- survfit(as.formula(form), data= dat.survival)
        cox.fit <- coxph( as.formula(gsub("\\+","*",form)), data= dat.survival)

        legend.treatment <- paste0("Treatment= ", unique(dat.survival$Treatment))
        legend.feature <- paste0(feature,"= ", 0:1)

        if (plot)
        {
                pdf(paste0("results/3_survival/interaction_",surv.type,"_treatment_", feature, ".pdf"))
                plot(KM.fit, col=cols, mark=NA, lwd =2, lty=rep(1:2,each=4),xlab="days", ylab=ylab)
                legend("topright", bty="", c(legend.treatment,"", legend.feature), col=c(cols,"white", "black", "black") ,cex=0.7, lwd=2, lty=c(rep(1,n_treatments),0,1,2))
                dev.off()
        }

        # out <- list(cox.fit=cox.fit, table=table(dat.survival[,"Treatment"], dat.survival[,feature]))

        summary.cox <- coef(summary(cox.fit))
        Index.interaction <- (n_treatments+1):(2*n_treatments -1)

        pval.interaction <- summary.cox[Index.interaction,5]

        return(pval.interaction)
}

survival_interaction_three_way <- function(dat.survival, surv.type, features, cols, plotting=T, out.dir="results/")
{
        # dat.survival <- dat.clinical.bis
        # surv.type <- "OS"
        # features <- c("RUNX1", "SRSF2", "TET2")
        # cols <- colors.survival_binary
        # plotting <-T
        # out.dir=paste0("results/datasets/",dataset)        

        final.dir <- paste0(out.dir, "/survival/",surv.type,"/three_way/")
        system(paste("mkdir -p", final.dir))

        if (surv.type=="OS")
        {
                ylab <- "Overall survival"
        } else if (surv.type=="EFS")
        {
                ylab <- "Disease-free survival"
        }

        iter <- combn(1:3, 2) # Enumerate all 

        for ( k in 1:3 )
        {
                pdf(paste0(final.dir,paste(features,collapse="_"), "_", k, ".pdf"), width=16, height=8)
                pos <- cbind(c(0,0,1,1), c(0,1,0,1))

                par(mfrow=c(1,4))
                for (j in 1:nrow(pos))
                {
                        dat.surv <- dat.survival[ dat.survival[  , features[iter[1,k]]]==pos[j,1] & dat.survival[  , features[iter[2,k]]]==pos[j,2],] 
                        last.feature <- features[setdiff(1:3, iter[,k])]

                        form <- paste0("Surv(", surv.type,", ", surv.type,"_Status) ~ factor(", last.feature ,")"  )

                        KM.fit <- survfit(as.formula(form), data= dat.surv)
                        cox.fit <- coxph(as.formula(form), data= dat.surv)

                        legend.names <- paste0(last.feature,"=",0:1)

                        if (plotting)
                        {
                                plot(KM.fit, col=cols, mark=NA, lwd =2, ,xlab="days", ylab=ylab, xlim=c(0,5000))
                                title(paste0(features[iter[1,k]], "=", pos[j,1], ",", features[iter[2,k]], "=", pos[j,2]  ))
                                text(x=4000,y=0.8,labels=paste0("P=",signif(summary(cox.fit)$coefficients[1,5],2) ))
                                legend("topright", bty="", legend.names, col=cols,cex=0.7, lwd=2, lty=1)
                        }
                }
                dev.off()
        }

}


survival_interaction_treatment <- function(dat.survival, surv.type, feature, cols, plot=T)
{
        n_treatments <- length(unique(dat.survival$Treatment))

        if (surv.type=="OS")
        {
                ylab <- "Overall survival"
        } else if (surv.type=="EFS")
        {
                ylab <- "Disease-free survival"
        }
        form <- paste0("Surv(", surv.type,"_followUp, ", surv.type,"_status) ~ factor(", feature,") + factor(Treatment)" )

        KM.fit <- survfit(as.formula(form), data= dat.survival)
        cox.fit <- coxph( as.formula(gsub("\\+","*",form)), data= dat.survival)

        legend.treatment <- paste0("Treatment= ", unique(dat.survival$Treatment))
        legend.feature <- paste0(feature,"= ", 0:1)

        if (plot)
        {
                pdf(paste0("results/3_survival/interaction_",surv.type,"_treatment_", feature, ".pdf"))
                plot(KM.fit, col=cols, mark=NA, lwd =2, lty=rep(1:2,each=4),xlab="days", ylab=ylab)
                legend("topright", bty="", c(legend.treatment,"", legend.feature), col=c(cols,"white", "black", "black") ,cex=0.7, lwd=2, lty=c(rep(1,n_treatments),0,1,2))
                dev.off()
        }

        # out <- list(cox.fit=cox.fit, table=table(dat.survival[,"Treatment"], dat.survival[,feature]))

        summary.cox <- coef(summary(cox.fit))
        Index.interaction <- (n_treatments+1):(2*n_treatments -1)

        pval.interaction <- summary.cox[Index.interaction,5]

        return(pval.interaction)
}



