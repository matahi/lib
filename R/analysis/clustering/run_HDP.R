run_HDP <- function(dat, 
                    out.dir="./results/clustering", 
                    cache=F, 
                    graphs=T, 
                    grouping=NULL, 
                    fast=T, 
                    col.colors=NULL,
                    legend.info=NULL)
{
        require(RColorBrewer)
        require(gplots)
        #### Libraries and data
        # require(hdp, lib.loc="src/lib/OSX/") ### This is the package that Nicola Roberts developed to extract mutation signatures - you can find it here https://github.com/nicolaroberts/hdp/tree/master
        require(hdp)
        # The following works with HDP version:


        # dat <- as.matrix(dat.binary)
        # out.dir <- "./results/clustering"
        # cache <- F
        # graphs <- F
        # fast <- T
        # grouping <- NULL
        # col.colors <- NULL
        # legend.info <- NULL


        system(paste0("mkdir -p ", out.dir,"/cache/"))

        print("Running HDP may take some time ...")

        ### Grouping format
        if (!is.null(grouping))
        {
                grouping$colors <- grouping$colors[grouping$features %in% colnames(dat)]
                grouping$features <- grouping$features[grouping$features %in% colnames(dat)]
                grouping$groups <- lapply(1:length(grouping$groups), function(n)
                                          {
                                                  return(grouping$groups[[n]][ grouping$groups[[n]] %in% colnames(dat)])
                                          })
        }

        #DPrun, cache=TRUE, cache.lazy=FALSE

        #DPsetup
        n <- ncol(dat)
        shape <- 1
        invscale <- 1
        hdp <- hdp_init(ppindex=0, #index of the parent DP for initial DP
                        cpindex=1, #index of alphaa and alphab for initial DP
                        hh=rep(1/n,n), #params for base distn (uniform Dirichlet)
                        alphaa=shape,
                        alphab=invscale)

        hdp <- hdp_adddp(hdp,
                         numdp=nrow(dat), # one DP for every sample in that cancer type
                         pp=1, # parent DP for group i is the i-th+1 overall DP because of the grandparent at position 1
                         cp=1) # index of alphaa and alphab for each DP

        # Assign the data from each patient to a child DP
        hdp <- hdp_setdata(hdp = hdp, dpindex=1:nrow(dat)+1, data=dat)

        # Activate the DPs with specified number of classes (signatures)
        hdp <- dp_activate(hdp, 1:(nrow(dat)+1), 7)

        #DP parameters
        # DPparam
        if (!fast)
        {
                burnin <- 5000
                postsamples <- 10000
                cpsamples <- 10
        } else 
        {
                burnin <- 2000
                postsamples <- 1500
                cpsamples <- 5 ## mm tries
        }
        spacebw <- 20

        if (!cache)
        {
                print("Calculating posterior distribution ...")
                output <- hdp_posterior(hdp, #activated hdp structure
                                        burnin=burnin,
                                        n=postsamples,
                                        space=spacebw,
                                        #doconparam=cpsamples,
                                        seed=42)

                #save(output, file=paste0("results/", out.dir,"/cache/HDP_output.RData"))
                save(output, file=paste0(out.dir,"/cache/HDP_output.RData"))
                print("Done!")
        } else 
        {
                print("Loading HDP output from cache...")
                output <- get(load(paste0(out.dir,"/cache/HDP_output.RData")))

        }

        # DP checks
        # plot_lik(output, bty="L") 
        # plot_numcluster(output, bty="L")
        # plot_data_assigned(output, bty="L")

        #posteriorMerged, cache=TRUE
        # posteriorMerged <- hdp_extract_signatures(output, prop.explained=0.99, cos.merge=0.99)


        if (!cache)
        {
                print("Merging posterior distribution ...")
                posteriorMerged <- hdp_extract_components(output) # cos.merge=0.9 # 2h
                save(posteriorMerged, file=paste0(out.dir,"/cache/HDP_posterior.RData"))
                #save(posteriorMerged, file=paste0("results/",out.dir,"/cache/HDP_posterior.RData"))
                print("Done!")

        } else
        {
                print("Loading posterior Merged output from cache...")
                posteriorMerged <- get(load(paste0(out.dir,"/cache/HDP_posterior.RData")))
        }

        # plot_comp_size(posteriorMerged, bty="L", lab=c(3, 5, 7))
        # n.components <- numcomp(posteriorMerged) + 1 
        #                                          + 1 (mut negative)

        # Find HDP classes
        ###################
        distn <- comp_dp_distn(posteriorMerged)$mean[-1,]


        HDP.class <- sapply(1:nrow(distn), function(n)
                            {
                                    if (all(is.na(distn[n,])))
                                    {
                                            return(1)
                                    } else
                                    {
                                            return(which.max(distn[n,])+1)
                                    }
                            })
        n.HDP <- length(unique(HDP.class))

        # table(HDP.class)[-1]/length(which(HDP.class!=1))
        # summary_distn

        HDP.post <- sapply(1:nrow(distn), function(n)
                           {
                                   if (all(is.na(distn[n,])))
                                   {
                                           return(1)
                                   } else
                                   {
                                           return(max(distn[n,]))
                                   }
                           })

        if (n.HDP>12)
        {
                source("./src/lib/R/plot/misc/gg_color_hue.R")
                cols <- gg_color_hue(n.HDP)
        } else
        {
                cols <- brewer.pal(max(n.HDP,3),"Set3")[1:n.HDP]
        }


        # Posterior distribution
        if (graphs)
        {
                #sig.title <- paste0("Signature HDP ", 2:(n.HDP))


                        # Signature associated with each class
                        #pdf(paste0('results/', out.dir,'/signatures.pdf'), width=16, height=9)
                        pdf(paste0(out.dir,'/signatures.pdf'), width=16, height=9)
                        # plot_comp_distn(posteriorMerged, cat_names=colnames(dat), plot_title=sig.title, col="skyblue3", col_nonsig="black")
                        plot_comp_distn(posteriorMerged, cat_names=colnames(dat), col="skyblue3", col_nonsig="black")
                        dev.off()

                #par(mfrow=c(ceiling(n.components/2 ),2), mar=c(3, 2, 2, 1))
                #plot_comp_distn(posteriorMerged, cat_names=colnames(dat), col="skyblue3", col_nonsig="grey70")



                ## Posterior probability distribution for each sample
                # plot_dp_comp_exposure(posteriorMerged, dpindices=2:(nrow(dat)+1), main_text="HDP classes", incl_numdata_plot=F,
                #                       col=RColorBrewer::brewer.pal(n.components, "Set3"))

                # names(cols) <- 1:n.HDP

                post.df <- data.frame(HDP.class=HDP.class, HDP.post=HDP.post)
                # Distribution of probability distribution per class assignment
                pdf(paste0(out.dir, '/distribution_assign.pdf'))
                #pdf(paste0('results/', out.dir, '/distribution_assign.pdf'))
                print(ggplot(post.df) + geom_boxplot(aes(x=factor(HDP.class),y=HDP.post, fill=factor(HDP.class))) + scale_fill_manual(name="HDP classes", values=cols)+ xlab("HDP classes") + ylab("Posterior probability") +  theme_bw())
                dev.off()


        }

        #### Reorder sample in classes 
        HDP.classes <- sort(unique(HDP.class))

        Index.bis <- lapply(1:length(HDP.classes), function(n)
                            {
                                    dat.bis <- data.frame(dat[which(HDP.class==HDP.classes[n]),,drop=F])

                                    order.gene <- order(colSums(dat.bis),decreasing=T)

                                    gene.order <- paste(paste0("-",colnames(dat.bis)[order.gene[1:10]]), collapse=",")

                                    Index.samples <- do.call(order,-dat.bis[,order.gene])
                            })

        final.order <- Reduce("c", lapply(1:length(HDP.classes), function(n)
                                          {
                                                  class.index <- which(HDP.class==HDP.classes[n])

                                                  index.final <- class.index[Index.bis[[n]]]
                                          }))


        ### Make a plot by classes
        for (k in 1:length(HDP.classes))
        {
                pdf(paste0(out.dir,'/probability_distrib_',k,'.pdf'), width=16, height=9)
                barplot(t(distn[HDP.class==HDP.classes[k],]), col= cols)
                dev.off()
        }

        # barplot(t(distn[order(HDP.class)[1:100],]), col= cols[2:length(cols)])

        # Plot final clustering
        if (graphs)
        {
                colmap <- colorRampPalette(c('grey','blue'))(100)
                clust.size <- as.vector(table(HDP.class))

                colors.HDP <- cols[match(HDP.class, sort(unique(HDP.class)))]

                if (is.null(col.colors))
                {
                        colors.final <- t(colors.HDP)
                } else
                {
                        colors.final <- rbind(colors.HDP, col.colors)
                }


                if (is.null(grouping))
                {
                        source("./lib/R/plot/heatmap/heatmap.3.R")
                        pdf(paste0(out.dir,'/HDP.pdf'))
                        heatmap.3(t(dat[final.order,]), scale='none', Rowv=NULL, Colv=NULL, dendrogram="none", 
                                  ColSideColors= t(colors.final[,final.order,drop=F]), trace="none", labCol=NA ,
                                  colsep=cumsum(clust.size)[-length(clust.size)], sepcolor="black",
                                  col=colmap)
                        if (!is.null(legend.info))
                        {
                                legend('bottomleft',fill= legend.info, legend = names(legend.info), border=NA) 
                        }
                        dev.off()

                } else
                {
                        if (!all(colnames(dat) %in% grouping$features))
                        {
                                print("Warning: Some features were missed in the grouping..")
                                missing.features <- colnames(dat)[!colnames(dat) %in% grouping$features]

                                grouping$features <- c(grouping$features, missing.features)
                                grouping$groups <- append(grouping$groups, list(Specific=missing.features))
                                grouping$colors <- c(grouping$colors, rep("red",length(missing.features)))
                        }



                        features.size <- sapply(grouping$groups,length)

                        source("./src/lib/R/plot/heatmap/heatmap.4.R")
                        pdf(paste0(out.dir,'/HDP_grouping.pdf'))
                        heatmap.4(t(dat[final.order,grouping$features]), scale='none', Rowv=NULL, Colv=NULL, dendrogram="none", 
                                  ColSideColors= t(colors.final[,final.order,drop=F]), trace="none", labCol=NA ,
                                  colsep=cumsum(clust.size)[-length(clust.size)], sepcolcolor="red",
                                  rowsep=cumsum(features.size)[-length(features.size)], seprowcolor="black",
                                  col.axis= grouping$colors,
                                  margins=c(5,10),
                                  col=colmap)
                        if (!is.null(legend.info))
                        {
                                legend('bottomleft',fill= legend.info, legend = names(legend.info), border=NA, bty="n") 
                                # legend('bottomleft',fill= legend.info, legend = names(legend.info), border=NA, bty="n", inset=-0.1) 
                        }
                        dev.off()


                }

        }

        return(list(class=HDP.class,order=final.order, prob=distn))

        # out.HDP <- list(class=HDP.class,order=final.order, prob=distn)
}
