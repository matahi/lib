analyze_HDP <- function(out.HDP, out.merge, out.clinical=NULL, grouping=NULL, out.dir="", display_legend=F)
{
        # n.HDP <- length(unique(out.HDP$class))
        # out.dir <- "clustering/ALL"
        # grouping <- grouping_AML_MDS
        # display_legend <- F

        # out.clinical <- dat.survival.impute

        # grouping <- NULL
        # out.merge <- dat.integrate.impute.bis 
        # out.clinical <- dat.clinical

        out.clinical <- out.merge$dat.clinical
        out.merge <- out.merge$dat

        # out.HDP 

        n.HDP <- length(unique(out.HDP$class))

        ## Lets process HDP

        prob <- apply(out.HDP$prob,1,max)

        manual.clusters <- c("No drivers",
                             "CEBPA_bi",
                             "Chromatin-spliceosome",
                             "NPM1",
                             "TP53 Complex",
                             "t_8_21",
                             "inv_16",
                             "IDH2_172",
                             "t_15_17")

        clusters.processed <- rep(NA, length(out.HDP$class))

        samples.manual <- vector("list", n.HDP)

        for (k in 1:length(manual.clusters))
        {
                 if (manual.clusters[k]=="No drivers")
                 {
                         # print("No drivers")
                         tmp <- apply(out.merge[out.HDP$class==k,],1,sum)

                         samples.manual[[k]] <- which(out.HDP$class==k)[which(tmp != 0)]

                         tmp.df <- data.frame(prob=prob[which(out.HDP$class==k)], class="No drivers", stringsAsFactors=F)
                         tmp.df[which(tmp != 0),"class"] <- "No class"

                         # pdf(paste0("results/",out.dir,"/no_drivers_manual.pdf"))
                         # print(ggplot(tmp.df, aes(x=class,y=prob)) + geom_boxplot() + geom_jitter())
                         # dev.off()

                         # Reassign classes
                         clusters.processed[out.HDP$class==k] <- "No drivers"
                         clusters.processed[which(out.HDP$class==k)[which(tmp != 0)]] <- "No class"

                } else if (manual.clusters[k]=="CEBPA_bi")
                {
                         # print("CEBPA_bi")

                        clusters.processed[out.HDP$class==k] <- "CEBPA_bi"

                        if (any(out.merge[out.HDP$class==k,"CEBPA_bi"]==0))
                        {
                                samples.manual[[k]] <- which(out.HDP$class==k)[which(out.merge[out.HDP$class==k,"CEBPA_bi"]==0)]

                                tmp.df <- data.frame(prob=prob[which(out.HDP$class==k)], class="CEBPA_bi", stringsAsFactors=F)
                                tmp.df[which(out.merge[out.HDP$class==k,"CEBPA_bi"]==0),"class"] <- "No-class"


                                pdf(paste0("results/",out.dir,"/CEBPA_bi_manual.pdf"))
                                print(ggplot(tmp.df, aes(x=class,y=prob)) + geom_boxplot() + geom_jitter())
                                dev.off()

                                # Reassign classes
                                clusters.processed[samples.manual[[k]]] <- "No class"
                        }

                } else if (manual.clusters[k]=="Chromatin-spliceosome")
                {

                         # print("Chromatin-spliceosome")
                        clusters.processed[out.HDP$class==k] <- "Chromatin-spliceosome"

                        gene.list <- c("ASXL1", "RUNX1", "MLL", "STAG2", "SF3B1", "U2AF1", "EZH2","SRSF2")
                        gene.present <- gene.list[gene.list %in% colnames(out.merge)]
                        tmp <- apply(out.merge[out.HDP$class==k,gene.present],1,sum)

                        samples.manual[[k]] <- which(out.HDP$class==k)[which(tmp < 2)]

                        tmp.df <- data.frame(prob=prob[which(out.HDP$class==k)], class="Chromatin-spliceosome", stringsAsFactors=F)
                        tmp.df[which(tmp < 2),"class"] <- "No-class"


                        pdf(paste0("results/",out.dir,"/Chromatin_spliceosome_manual.pdf"))
                        print(ggplot(tmp.df, aes(x=class,y=prob)) + geom_boxplot() + geom_jitter())
                        dev.off()

                        clusters.processed[samples.manual[[k]]] <- "No class"

                } else if (manual.clusters[k]=="NPM1")
                {

                        clusters.processed[out.HDP$class==k] <- "NPM1"

                        if (any(out.merge[out.HDP$class==k,"NPM1"]==0))
                        {
                                samples.manual[[k]] <- which(out.HDP$class==k)[which(out.merge[out.HDP$class==k,"NPM1"]==0)]

                                tmp.df <- data.frame(prob=prob[which(out.HDP$class==k)], class="NPM1", stringsAsFactors=F)
                                tmp.df[which(out.merge[out.HDP$class==k,"NPM1"]==0),"class"] <- "No-class"


                                pdf(paste0("results/",out.dir,"/NPM1_manual.pdf"))
                                print(ggplot(tmp.df, aes(x=class,y=prob)) + geom_boxplot() + geom_jitter())
                                dev.off()
                        }

                        clusters.processed[samples.manual[[k]]] <- "No class"

                } else if (manual.clusters[k]=="TP53 Complex")
                {

                        clusters.processed[out.HDP$class==k] <- "TP53 Complex"

                        # if (any(out.merge[out.HDP$class==k,"NPM1"]==0))
                        # {
                        #         samples.manual[[k]] <- which(out.HDP$class==k)[which(out.merge[out.HDP$class==k,"NPM1"]==0)]

                        #         tmp.df <- data.frame(prob=prob[which(out.HDP$class==k)], class="NPM1", stringsAsFactors=F)
                        #         tmp.df[which(out.merge[out.HDP$class==k,"NPM1"]==0),"class"] <- "No-class"


                        #         pdf(paste0("results/",out.dir,"/NPM1_manual.pdf"))
                        #         print(ggplot(tmp.df, aes(x=class,y=prob)) + geom_boxplot() + geom_jitter())
                        #         dev.off()
                        # }

                        # clusters.processed[samples.manual[[k]]] <- "No class"

                } else if (manual.clusters[k]=="t_8_21")
                {
                        clusters.processed[out.HDP$class==k] <- "t_8_21"

                        if (any(out.merge[out.HDP$class==k,"t_8_21"]==0))
                        {
                                samples.manual[[k]] <- which(out.HDP$class==k)[which(out.merge[out.HDP$class==k,"t_8_21"]==0)]

                                tmp.df <- data.frame(prob=prob[which(out.HDP$class==k)], class="t_8_21", stringsAsFactors=F)
                                tmp.df[which(out.merge[out.HDP$class==k,"t_8_21"]==0),"class"] <- "No-class"


                                pdf(paste0("results/",out.dir,"/t8_21_manual.pdf"))
                                print(ggplot(tmp.df, aes(x=class,y=prob)) + geom_boxplot() + geom_jitter())
                                dev.off()

                                clusters.processed[samples.manual[[k]]] <- "No class"
                        }
                } else if (manual.clusters[k]=="inv_16")
                {
                        clusters.processed[out.HDP$class==k] <- "inv_16"

                        if (any(out.merge[out.HDP$class==k,"inv_16"]==0))
                        {
                                samples.manual[[k]] <- which(out.HDP$class==k)[which(out.merge[out.HDP$class==k,"inv_16"]==0)]

                                tmp.df <- data.frame(prob=prob[which(out.HDP$class==k)], class="inv_16", stringsAsFactors=F)
                                tmp.df[which(out.merge[out.HDP$class==k,"inv_16"]==0),"class"] <- "No-class"


                                pdf(paste0("results/",out.dir,"/inv16_manual.pdf"))
                                print(ggplot(tmp.df, aes(x=class, y=prob)) + geom_boxplot() + geom_jitter())
                                dev.off()

                                clusters.processed[samples.manual[[k]]] <- "No class"
                        }
                } else if (manual.clusters[k]=="IDH2_172")
                {
                        clusters.processed[out.HDP$class==k] <- "IDH2_172"
                } else if (manual.clusters[k]=="t_15_17")
                {
                        clusters.processed[out.HDP$class==k] <- "t_15_17"

                        if (any(out.merge[out.HDP$class==k,"t_15_17"]==0))
                        {
                                samples.manual[[k]] <- which(out.HDP$class==k)[which(out.merge[out.HDP$class==k,"t_15_17"]==0)]

                                tmp.df <- data.frame(prob=prob[which(out.HDP$class==k)], class="t_15_17", stringsAsFactors=F)
                                tmp.df[which(out.merge[out.HDP$class==k,"t_15_17"]==0),"class"] <- "No-class"


                                pdf(paste0("results/",out.dir,"/t15_17_manual.pdf"))
                                print(ggplot(tmp.df, aes(x=class, y=prob)) + geom_boxplot() + geom_jitter())
                                dev.off()

                                clusters.processed[samples.manual[[k]]] <- "No class"
                        }
                }         
        }

        table(clusters.processed)

        #### Do the ambiguous classification
        ambiguous.cutoff <- 0.2 

        out.diff <- Reduce("rbind",lapply(1:nrow(out.HDP$prob), function(n)
                            {
                                    if (any(is.na(out.HDP$prob[n,])))
                                    {
                                            return(data.frame(diff=NA,class2=NA))
                                    } else
                                    {
                                            max1 <- max(out.HDP$prob[n,])
                                            max2 <- max(out.HDP$prob[n,-which(out.HDP$prob[n,]==max1)])

                                            class2 <- order(out.HDP$prob[n,], decreasing=T)[2]

                                            return(data.frame(diff=max1-max2,
                                                              class2=class2))
                                    }
                            }))


        diff.df <- data.frame(prob=out.diff$diff, class=out.HDP$class,class2=out.diff$class2+1)

        pdf(paste0("results/", out.dir,"/prob_diff.pdf"))
        print(ggplot(diff.df) + geom_boxplot(aes(x=factor(class), y= prob)))
        dev.off()

        ambiguous.df <- diff.df[diff.df$prob< ambiguous.cutoff,]

        table(ambiguous.df$class)

        library(knitr)
        # table(ambiguous.df$class, ambiguous.df$class2+10)
        kable(table(ambiguous.df$class, ambiguous.df$class2), format="latex")

        clusters.processed[ diff.df$prob < ambiguous.cutoff] <- "Ambiguous"

        kable(table(clusters.processed, out.HDP$class), format="latex")
        table(clusters.processed)

        ### Use the clusters.processed

        clusters.processed <- factor(clusters.processed, levels=c("t_8_21",
                                                                  "inv_16",
                                                                  "t_15_17",
                                                                  "NPM1",
                                                                  "TP53 Complex",
                                                                  "Chromatin-spliceosome",
                                                                  "CEBPA_bi",
                                                                  "IDH2_172",
                                                                  "Ambiguous",
                                                                  "No class",
                                                                  "No drivers"))

        Index.bis <- lapply(1:length(unique(clusters.processed)), function(n)
                            {
                                    dat.bis <- data.frame(out.merge[which(clusters.processed==levels(clusters.processed)[n]),,drop=F])

                                    order.gene <- order(colSums(dat.bis),decreasing=T)

                                    gene.order <- paste(paste0("-",colnames(dat.bis)[order.gene[1:10]]), collapse=",")

                                    Index.samples <- do.call(order,-dat.bis[,order.gene])
                            })

        final.order <- Reduce("c", lapply(1:length(unique(clusters.processed)), function(n)
                                          {
                                                  class.index <- which(clusters.processed==levels(clusters.processed)[n])

                                                  index.final <- class.index[Index.bis[[n]]]
                                          }))



        # AML-MDS
        require(RColorBrewer)
        colmap <- colorRampPalette(c('grey','blue'))(100)
        clust.size <- as.vector(table(clusters.processed))
        source("src/lib/plot/gg_color_hue.R")
        cols <- brewer.pal(length(levels(clusters.processed)),"Set3")
        colors.HDP <- rbind(cols[as.numeric(clusters.processed)],cols[as.numeric(clusters.processed)])

        # colors.dataset.heatmap <- colors.datasets[as.character(out.merge$dat.clinical$dataset)]
        # colors.WHO.heatmap <- colors.MDS_type[out.merge$dat.clinical$WHO]

        # TODO: Add IPSS

        # Plot heatmap clustering
        # source("src/lib/plot/heatmap.3.R")
        if (is.null(grouping))
        {
                # source("src/lib/plot/heatmap.4.R")
                library(gplots)
                #pdf(paste0('results/', out.dir,'/cluster_post_processed.pdf'), width=16, height=9)
                #pdf(paste0('results/', out.dir,'/cluster_post_processed.pdf'), width=16, height=30)
                pdf(paste0('results/', out.dir,'/cluster_post_processed.pdf'))
                heatmap.2(t(out.merge[final.order,]), scale="none", Rowv=NULL, Colv=NULL, dendrogram="none",
                          trace="none", labCol=NA, key=F,
                          colsep=cumsum(clust.size)[-length(clust.size)], sepcolor="black",
                          ColSideColors=colors.HDP[1,final.order],
                          # ColSideColors=t(colors.HDP[,final.order]),
                          col=colmap)
                dev.off()

        } else
        {
                source("src/lib/plot/heatmap.4.R")
                library(gplots)
                pdf(paste0('results/', out.dir,'/cluster.pdf'), width=16, height=9)
                heatmap.4(t(out.merge[out.HDP$order,]), scale="none", Rowv=NULL, Colv=NULL, dendrogram="none",
                          trace="none", labCol=NA, key=F,
                          colsep=cumsum(clust.size)[-length(clust.size)], sepcolor="black",
                          #ColSideColors=t(colors.heatmap[,out.HDP$order]),
                          #col.axis= colors.feature_type[factor(grouping)],
                          col=colmap)
                if (display_legend)
                        legend("bottomleft", legend=unique(out.merge$dat.clinical$WHO), fill=unique(colors.WHO.heatmap), bty="n")

                dev.off()

        }

        # Survival by HDP class
        HDP.df <- data.frame(HDP=clusters.processed, out.clinical)

        # Survival
        source("src/lib/analysis/run_univariate_survival.R")
        pp <- univariate_survival_HDP(HDP.df, surv.type="OS", cols=cols, stats=F, out.dir=out.dir, cluster.names=levels(clusters.processed))

        # pp <- univariate_survival_HDP(HDP.df, surv.type="EFS", cols=cols, stats=T, out.dir=out.dir)

        # Other clinical variables
        pdf(paste0('results/', out.dir,'/HDP_Age.pdf'))
        print(ggplot(HDP.df) + geom_boxplot(aes(x=HDP,y=Age, fill=HDP)) + scale_fill_manual(values=cols) + theme_bw()
                             + theme(axis.text.x=element_text(angle=90)))
        dev.off()

        table(clusters.processed, useNA="always")

}
