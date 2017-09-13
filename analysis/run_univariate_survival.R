require(survival)


univariate_survival <- function(dat.survival, feature, surv.type="OS",  cols=NULL, stats=T, cens=NULL,  out.dir="results/")
{
        # dat.survival <- dat.clinical.ter
        # surv.type <- "OS" 
        # feature <- "chronology"
        # cols <- cols.chronology
        # stats <- F
        # cens <- NULL
        # out.dir <- "./results/supp/Ross"

        final.dir <- paste0(out.dir, "/survival/", surv.type,"/univariate/")
        system(paste0("mkdir -p ",  final.dir))

        if (surv.type=="OS")
        {
                ylab <- "Overall survival"
        } else if (surv.type=="EFS")
        {
                ylab <- "Disease-free survival"
        }

        if (!is.null(cens))
        {
                surv.status <- paste0(surv.type,"_Status")
                dat.survival[which(dat.survival[,surv.type] >= cens), surv.type] <- cens
                dat.survival[which(dat.survival[,surv.type] >= cens), surv.status] <- F
        }

        form <- paste0("Surv(",surv.type,",", surv.type,"_Status) ~ factor(", feature,")")

        KM.fit <- survfit(as.formula(form), data= dat.survival)
        cox.fit <- coxph(as.formula(form), data= dat.survival)

        # legend.names
        legend.names <- levels(factor(dat.survival[,feature]))

        # quartz()
        # Plot
        pdf(paste0(final.dir,feature ,".pdf"))
        plot(KM.fit, col=cols, mark=NA, lwd =2, ,xlab="days", ylab =ylab)
        #plot(KM.fit, col=cols, mark=NA, lwd =2, ,xlab="years", ylab =ylab)
        legend("topright", bty="", legend.names, col=cols,cex=0.7, lwd=2, lty=1)
        dev.off()

        # # median info
        # if (any(is.na(dat.survival[,feature])))
        # {
        #         surv.out <- read.table(textConnection(capture.output(KM.fit)),skip=3,header=TRUE)
        # } else
        # {
        #         surv.out <- read.table(textConnection(capture.output(KM.fit)),skip=3,header=TRUE)
        #         # KM.fit
        # }

        # if (stats)
        # {
        #         Absent <- round(surv.out$median[1]/30,1)
        #         Present <- round(surv.out$median[2]/30,1)
        #         freq <- round(100*sum(dat.survival[,feature],na.rm=T)/(nrow(dat.survival)-  sum(is.na(dat.survival[,feature]))),1)
        #         pval <- summary(cox.fit)$coefficients[5]

        #         return(data.frame(feature=feature, Absent=Absent, Present=Present, freq, pval=pval, surv.type=surv.type))
        # } else
        # {
        #         # return(cox.fit)
        #         return(summary(cox.fit)$conf.int)
        # }

        return(list(cox.fit,KM.fit,form))
}

univariate_survival_binary <- function(dat.survival, surv.type, feature, cols=c("#d7191c","#2c7bb6"), stats=T, cens=NULL,  out.dir="results/")
{
        final.dir <- paste0(out.dir, "/survival/", surv.type,"/univariate/")
        system(paste0("mkdir -p ",  final.dir))

        if (surv.type=="OS")
        {
                ylab <- "Overall survival"
        } else if (surv.type=="EFS")
        {
                ylab <- "Disease-free survival"
        }

        if (!is.null(cens))
        {
                surv.status <- paste0(surv.type,"_Status")
                dat.survival[which(dat.survival[,surv.type] >= cens), surv.type] <- cens
                dat.survival[which(dat.survival[,surv.type] >= cens), surv.status] <- F
        }

        form <- paste0("Surv(",surv.type,",", surv.type,"_Status) ~ factor(", feature,")")

        KM.fit <- survfit(as.formula(form), data= dat.survival)
        cox.fit <- coxph(as.formula(form), data= dat.survival)

        # legend.names
        legend.names <- paste0(feature,"=",0:1)

        # Plot
        pdf(paste0(final.dir,feature ,".pdf"))
        plot(KM.fit, col=cols, mark=NA, lwd =2, ,xlab="days", ylab =ylab)
        legend("topright", bty="", legend.names, col=cols,cex=0.7, lwd=2, lty=1)
        dev.off()

        # median info
        if (any(is.na(dat.survival[,feature])))
        {
                surv.out <- read.table(textConnection(capture.output(KM.fit)),skip=3,header=TRUE)
        } else
        {
                surv.out <- read.table(textConnection(capture.output(KM.fit)),skip=3,header=TRUE)
                # KM.fit
        }

        if (stats)
        {
                Absent <- round(surv.out$median[1]/30,1)
                Present <- round(surv.out$median[2]/30,1)
                freq <- round(100*sum(dat.survival[,feature],na.rm=T)/(nrow(dat.survival)-  sum(is.na(dat.survival[,feature]))),1)
                pval <- summary(cox.fit)$coefficients[5]

                return(data.frame(feature=feature, Absent=Absent, Present=Present, freq, pval=pval, surv.type=surv.type))
        } else
        {
                # return(cox.fit)
                return(summary(cox.fit)$conf.int)
        }

}

run_univariate_survival <- function(dat.survival, surv.type, feature_list, cols,stats=T, out.dir="results/")
{
        # dat.survival <- dat.survival.impute
        # surv.type <- "OS"
        # feature_list <- gene_frequent
        # cols <- colors.survival_binary

        final.dir <- paste0(out.dir, "/survival/", surv.type,"/univariate/")
        system(paste0("mkdir -p ",  final.dir))
        # i <- 1

        univariate.out <- Reduce('rbind', lapply(1:length(feature_list), function(i)
                                                 {
                                                         out <- univariate_survival(dat.survival, surv.type=surv.type,feature=feature_list[i] , cols=cols, stats=stats, out.dir=final.dir )
                                                 }))

        if (!stats)
                return(univariate.out)

        univariate.out$qval <- p.adjust(univariate.out$pval, method="fdr")

        univariate.out$pval.final <- sapply(1:nrow(univariate.out), function(n)
                                            {
                                                    if (is.na(univariate.out$qval[n]))
                                                    {
                                                            return(NA)
                                                    } else if (univariate.out$qval[n]>0.1)
                                                    {
                                                            return("NS")
                                                    } else if (univariate.out$qval[n]<0.05)
                                                    {
                                                            return("< .05")
                                                    } else
                                                    {
                                                            return("< .1")
                                                    }
                                            })
        rownames(univariate.out) <- univariate.out$feature

        return(univariate.out)
}

univariate_survival_clonality <- function(dat.survival, surv.type, Gene, cols=colors.survival_binary, stats=F)
{

        if (surv.type=="OS")
        {
                ylab <- "Overall survival"
        } else if (surv.type=="EFS")
        {
                ylab <- "Disease-free survival"
        }

        form <- paste0("Surv(",surv.type,",", surv.type,"_Status) ~ clonality")

        KM.fit <- survfit(as.formula(form), data= dat.survival)
        cox.fit <- coxph(as.formula(form), data= dat.survival)

        # legend.names
        legend.names <- paste0(Gene," ",c("subclonal", "clonal"))

        # Plot
        pdf(paste0("results/5_clonality/",Gene,"/", surv.type,"_clonality.pdf"))
        plot(KM.fit, col=cols, mark=NA, lwd =2, ,xlab="days", ylab =ylab)
        legend("topright", bty="", legend.names, col=cols,cex=0.7, lwd=2, lty=1)
        dev.off()

        # median info
        if (any(is.na(dat.survival[,"clonality"])))
        {
                surv.out <- read.table(textConnection(capture.output(KM.fit)),skip=3,header=TRUE)
        } else
        {
                surv.out <- read.table(textConnection(capture.output(KM.fit)),skip=2,header=TRUE)
        }

        if (stats)
        {
                Absent <- round(surv.out$median[1]/30,1)
                Present <- round(surv.out$median[2]/30,1)
                freq <- round(100*sum(dat.survival[,"clonality"],na.rm=T)/(nrow(dat.survival)-  sum(is.na(dat.survival[,"clonality"]))),1)
                pval <- summary(cox.fit)$coefficients[5]

                return(data.frame(Gene=Gene, Absent=Absent, Present=Present, freq, pval=pval, surv.type=surv.type))
        } else
        {
                return(cox.fit)
        }

}


univariate_survival_HDP <- function(dat.survival, surv.type="OS", cols=NULL, stats=T, out.dir="3_survival", cluster.names=NULL)
{
        require(survival)
        # dat.survival <- HDP.df
        # surv.type <- 'OS'
        # cols <- cols
        # names(cols) <- cluster.names

        if (is.numeric(dat.survival$HDP))
        {
                HDP.bis <- paste0("HDP", dat.survival$HDP)
                cluster.names <- paste0("HDP", sort(unique(dat.survival$HDP)))
                dat.survival$HDP <- factor(HDP.bis, levels=cluster.names)
        } else
        {
                cluster.names <- levels(dat.survival$HDP)
        }

        if (any(dat.survival[,surv.type]<0,na.rm=T))
        {
                print("some negative OS times")
                dat.survival[which(dat.survival[,surv.type]<0) ,surv.type] <- - dat.survival[which(dat.survival[,surv.type]<0) ,surv.type]
        }

        if (surv.type=="OS")
        {
                ylab <- "Overall survival"
        } else if (surv.type=="EFS")
        {
                ylab <- "Disease-free survival"
        }

        #form <- paste0("Surv(",surv.type,",", surv.type,"_Status) ~ factor(HDP)")
        form <- paste0("Surv(",surv.type,",", surv.type,"_Status) ~ HDP")

        KM.fit <- survfit(as.formula(form), data= dat.survival)
        cox.fit <- coxph(as.formula(form), data= dat.survival)

        # legend.names
        #legend.names <- paste0("cluster =",1:length(unique(dat.survival$HDP)))
        if (is.null(cluster.names))
        {
                legend.names <- paste0("cluster= ",levels(dat.survival$HDP))
        } else if (length(cluster.names)!= length(levels(dat.survival$HDP)))
        {
                stop("cluster names is different than number of HDP clusters")
        } else
        {
                legend.names <- cluster.names
        }

        if (is.null(cols))
        {
                library(RColorBrewer)
                source("./src/lib/plot/gg_color_hue.R")
                cols <- gg_color_hue(length(levels(dat.survival$HDP)))
                names(cols) <- levels(dat.survival$HDP)
        }

        cols <- brewer.pal(length(cluster.names),"Set3")
        names(cols) <- cluster.names

        # Plot
        pdf(paste0(out.dir,"/",surv.type,"_HDP.pdf"))
        plot(KM.fit, col=cols[levels(dat.survival$HDP)], mark=NA, lwd =2, ,xlab="days", ylab =ylab)
        # plot(KM.fit, col=cols[unique(dat.survival$HDP)], mark=NA, lwd =2, ,xlab="days", ylab =ylab)
        legend("topright", bty="", legend=legend.names, col=cols[levels(dat.survival$HDP)],cex=0.7, lwd=2, lty=1)
        dev.off()

        # # median info
        # if (any(is.na(dat.survival[,"HDP"])))
        # {
        #         surv.out <- read.table(textConnection(capture.output(KM.fit)),skip=3,header=TRUE)
        # } else
        # {
        #         # surv.out <- read.table(textConnection(capture.output(KM.fit)),skip=2,header=TRUE)
        #         surv.out <- read.table(textConnection(capture.output(KM.fit)),skip=3,header=TRUE)
        # }

        # if (stats)
        # {
        #         Absent <- round(surv.out$median[1]/30,1)
        #         Present <- round(surv.out$median[2]/30,1)
        #         freq <- round(100*sum(dat.survival[,"HDP"],na.rm=T)/(nrow(dat.survival)-  sum(is.na(dat.survival[,"HDP"]))),1)
        #         pval <- summary(cox.fit)$coefficients[5]

        #         return(data.frame(feature="HDP", Absent=Absent, Present=Present, freq, pval=pval, surv.type=surv.type))
        # } else
        # {
        #         # return(cox.fit)
        #         return(summary(cox.fit)$conf.int)
        # }

}


