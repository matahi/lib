####################################################
###  Hotspot analysis
####################################################

dat.genetics <- Reduce("rbind", dat.genetics.list[1:5])
gene.list <- unique(dat.genetics$Gene)

# gene.list.bis <- gene.list[-c(86,108,109,134,139)]       # 

# gene.list.bis <- gene.list[c(86,108,109,134,139)]       # 

all.consequences <- list("complex"=c("Complex"),
                         "frameshift"=c("frameshift_variant", "indel", "inframe_codon_gain", "inframe_codon_loss", "inframe_variant", "initiator_codon_change", "intron_variant", "transcript_variant"),
                         "missense"=c("non_synonymous_codon"),
                         "nonsense"=c("stop_gained", "stop_lost"),
                         "synonymous"=c("synonymous_codon"))

# gene.list <- gene.list.bis

for (n in 1:length(gene.list))
{
        print(n)
        gene <- gene.list[n]
        print(gene)
        #system(paste0("mkdir -p ./results/analysis/gene/",gene))

        dat.genetics.gene <- dat.genetics[which(dat.genetics$Gene==gene),]
        compare <- "disease"

        dat.clinical <- Reduce("rbind", dat.clinical.list[1:5])

        dat.genetics.gene$PROTEIN_POS <- sapply(1:nrow(dat.genetics.gene), function(n)
                                                {
                                                        PROTEIN_NUM <- strsplit(dat.genetics.gene$PROTEIN_CHANGE[n],"[^0-9]+")[[1]][2]
                                                        return(as.numeric(PROTEIN_NUM))
                                                })


        lollipop <- F

        if (lollipop)
        {
                #######
                # TODO:
                #######
                print("1. Working on lollipop plot")
                library(biomaRt)
                library(ggplot2)
                ## 1. get pfam
                # i. get transcriptID
                ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset= "hsapiens_gene_ensembl")
                cds.infos <- getBM(attributes=c("coding", "cds_length","ensembl_transcript_id"), filters = "hgnc_symbol", values = gene, mart= ensembl)
                # Removing non available sequences
                cds.infos <- cds.infos[!is.na(cds.infos$cds_length),]
                # Removing sequences that are not of length /3
                cds.infos <- cds.infos[(cds.infos$cds_length%%3==0),]
                if (nrow(cds.infos)!=0)
                {
                        # Keeping longest transcript
                        longest.transcript.list <- cds.infos$ensembl_transcript_id[(cds.infos$cds_length==max(cds.infos$cds_length))]
                        # If several keep the first one by default
                        transcript.taken <- longest.transcript.list[1] 
                        # Recover sequence
                        codingSeq <- cds.infos$coding[cds.infos$ensembl_transcript_id==transcript.taken]
                        # Transform sequence into AA
                        residueSeq <- GenVisR:::lolliplot_DNAconv(codingSeq, to = "residue")
                        # Filter out non AA codes
                        source("./src/lib/plot/lollipop/filter_residue.R")
                        residueSeq <- filter_residue(residueSeq)
                        # Get length
                        proteinLength <- length(residueSeq)
                        # Recover domains
                        protein_domain <- getBM(attributes=c("interpro_description", "interpro_start","interpro_end"), filters=c("ensembl_transcript_id"), values=transcript.taken, mart=ensembl)
                        if(all(is.na(protein_domain$interpro_description)))
                        {
                                protein_domain <- data.frame(interpro_description=gene, interpro_start=1, interpro_end=proteinLength)
                        }
                } else
                {
                        print("Careful: protein is invalid")
                        next
                }


                # # library(GenVisR)
                # source("./src/lib/plot/lollipop/constructGene.R")
                # gene_data <- constructGene(gene, protein_domain,proteinLength)

                ##################
                # TODO: Method 1
                ##################
                # Trackviewer
                library(trackViewer)

                # i) Mutation infos
                #####################
                #SNP <- c(20,200,500,500)
                dat.genetics.disease <- lapply(1:length(unique(dat.genetics.gene$disease)), function(n)
                                               {
                                                       return(dat.genetics.gene[dat.genetics.gene$disease==unique(dat.genetics.gene$disease)[n],,drop=F])
                                               })
                names(dat.genetics.disease) <- unique(dat.genetics.gene$disease)

                ###
                # check if only intronic mutations
                ###
                intronic.check <- all( sapply(1:length(dat.genetics.disease), function(n)
                                              {
                                                      PROTEIN_CHANGE <- dat.genetics.disease[[n]]$PROTEIN_CHANGE

                                                      PROTEIN_POS <- sapply(1:length(PROTEIN_CHANGE), function(n)
                                                                            {
                                                                                    PROTEIN_NUM <- strsplit(as.character(PROTEIN_CHANGE[n]),"[^0-9]+")[[1]][2]
                                                                                    return(as.numeric(PROTEIN_NUM))
                                                                            })

                                                      all(is.na(PROTEIN_POS))
                                              }))


                # Remove intronic mutations
                dat.genetics.disease <- lapply(1:length(unique(dat.genetics.disease)), function(n)
                                               {
                                                       PROTEIN_CHANGE <- dat.genetics.disease[[n]]$PROTEIN_CHANGE

                                                       PROTEIN_POS <- sapply(1:length(PROTEIN_CHANGE), function(n)
                                                                             {
                                                                                     PROTEIN_NUM <- strsplit(as.character(PROTEIN_CHANGE[n]),"[^0-9]+")[[1]][2]
                                                                                     return(as.numeric(PROTEIN_NUM))
                                                                             })

                                                       tmp <- dat.genetics.disease[[n]]
                                                       if (any(is.na(PROTEIN_POS)))
                                                       {
                                                               tmp <- tmp[-which(is.na(PROTEIN_POS)),]
                                                       }

                                                       return(tmp)
                                               })
                names(dat.genetics.disease) <- unique(dat.genetics.gene$disease)

                # Remove empty fields
                if (any(sapply(dat.genetics.disease,nrow)==0))
                        dat.genetics.disease <- dat.genetics.disease[-which(sapply(dat.genetics.disease,length)==0)]

                ###
                color.type <- "CONSEQUENCE"
                color.panel <- NULL
                if (color.type == "CONSEQUENCE")
                {
                        color.panel <- colors.consequence_type
                } else if (color.type == "type")
                {
                        color.panel <- colors.variant_type
                }

                cutoff <- 20
                gene.hotspots <- names(table(dat.genetics.gene$PROTEIN_CHANGE))[table(dat.genetics.gene$PROTEIN_CHANGE)>=cutoff]

                if (length(grep("\\?", gene.hotspots))>0)
                        gene.hotspots <- gene.hotspots[-grep("\\?", gene.hotspots)]

                if (length(dat.genetics.disease)>0)
                {
                        sample.gr.disease <- lapply(1:length(dat.genetics.disease), function(n)
                                                    {
                                                            dat.genetics.gene <- dat.genetics.disease[[n]]
                                                            mut.infos <- data.frame(table(PROTEIN_CHANGE=dat.genetics.gene$PROTEIN_CHANGE))

                                                            mut.infos$PROTEIN_POS <- sapply(1:nrow(mut.infos), function(n)
                                                                                            {
                                                                                                    PROTEIN_NUM <- strsplit(as.character(mut.infos$PROTEIN_CHANGE[n]),"[^0-9]+")[[1]][2]
                                                                                                    return(as.numeric(PROTEIN_NUM))
                                                                                            })

                                                            if (any(is.na(mut.infos$PROTEIN_POS )))
                                                                    mut.infos <- mut.infos[!is.na(mut.infos$PROTEIN_POS),]

                                                            sample.gr <- GRanges(gene, IRanges(mut.infos$PROTEIN_POS, width=1, names=mut.infos$PROTEIN_CHANGE))
                                                            # sample.gr <- GRanges(gene, IRanges(mut.infos$PROTEIN_POS, width=1))
                                                            sample.gr$disease <- names(dat.genetics.disease)[n]

                                                            if (is.null(color.type))
                                                            {
                                                                    sample.gr$color <- 1
                                                            } else if (color.type == "CONSEQUENCE")
                                                            {
                                                                    color.panel <- colors.consequence_type
                                                                    consequences <- dat.genetics.gene$CONSEQUENCE[match(names(sample.gr), dat.genetics.gene$PROTEIN_CHANGE)]
                                                                    sample.gr$type[consequences %in% all.consequences[["missense"]]] <- "missense"
                                                                    sample.gr$type[consequences %in% all.consequences[["nonsense"]]] <- "nonsense"
                                                                    sample.gr$type[consequences %in% all.consequences[["frameshift"]]] <- "frameshift"
                                                                    sample.gr$type[consequences %in% all.consequences[["complex"]]] <- "complex"
                                                                    sample.gr$type[consequences %in% all.consequences[["synonymous"]]] <- "synonymous"
                                                                    sample.gr$color <- color.panel[sample.gr$type]



                                                            } else if (color.type =="type")
                                                            {
                                                                    color.panel <- colors.variant_type
                                                                    #sample.gr$type <- 
                                                            }

                                                            sample.gr$score <- mut.infos$Freq
                                                            if (max(sample.gr$score) <=10) # circumvent the stacking
                                                                    sample.gr$score <- sample.gr$score + 0.0001
                                                            sample.gr$dashline.col <- sample.gr$color

                                                            # Remove those that are not hotspots

                                                            sample.gr$dashline.col[!(names(sample.gr) %in% gene.hotspots)] <- NA

                                                            # Remove names
                                                            #names(sample.gr)[!(names(sample.gr) %in% gene.hotspots)] <- NA
                                                            names(sample.gr)[!(names(sample.gr) %in% gene.hotspots)] <- " "

                                                            return(sample.gr)
                                                    })
                        names(sample.gr.disease) <- names(dat.genetics.disease)
                } else
                {
                        sample.gr.disease <- GRanges()
                }

                legends <- list(labels=names(color.panel), fill=color.panel)

                # ii) Protein infos
                ##################
                features.list <- lapply(1:length(dat.genetics.disease), function(n)
                                        {
                                                if (length(dat.genetics.disease)!=1 & (n!=1))
                                                {
                                                        IRange.info <- IRanges(protein_domain$interpro_start,
                                                                               width=protein_domain$interpro_end - protein_domain$interpro_start + 1)

                                                } else
                                                {
                                                        IRange.info <- IRanges(protein_domain$interpro_start,
                                                                               width=protein_domain$interpro_end - protein_domain$interpro_start + 1,
                                                                               names=protein_domain$interpro_description)
                                                }

                                                features <- GRanges(gene, IRange.info)

                                                if (length(unique(protein_domain$interpro_description)) <= 12)
                                                {
                                                        colors.domains <- brewer.pal(length(unique(protein_domain$interpro_description)),"Set3")
                                                } else
                                                {
                                                        source("./src/lib/plot/gg_color_hue.R")
                                                        colors.domains <- gg_color_hue(length(unique(protein_domain$interpro_description)))
                                                }
                                                names(colors.domains) <- unique(protein_domain$interpro_description)
                                                features$fill <- colors.domains[protein_domain$interpro_description]


                                                freq.domain <- data.frame(table(protein_domain$interpro_description))
                                                features$height <- sapply(1:nrow(protein_domain), function(n){1/(10*freq.domain$Freq[freq.domain$Var1==protein_domain$interpro_description[n]])})

                                                return(features)
                                        })

                if (length(dat.genetics.disease)==0)
                        features.list <- features.list[1]

                # iii) Plot infos
                #################
                if (proteinLength>=600)
                {
                        xaxis.infos <- c(1, seq(200, proteinLength, by=200)) # xticks
                } else if (proteinLength>=200)
                {
                        xaxis.infos <- c(1, seq(100, proteinLength, by=100)) # xticks
                } else if (proteinLength>=50)
                {
                        xaxis.infos <- c(1, seq(50, proteinLength, by=50)) # xticks
                } else
                {
                        xaxis.infos <- c(1, proteinLength) # xticks
                }

                pdf(paste0("./results/analysis/lollipop/", gene,".pdf"), width=7, height=14.5)
                # pdf(paste0("./results/analysis/gene/",gene,"/lollipop.pdf"), width=7, height=14.5)
                lolliplot(sample.gr.disease, 
                          features.list, 
                          ranges=GRanges(gene, IRanges(0,proteinLength+1)),
                          xaxis=xaxis.infos,
                          type= "circle",
                          legend=legends, # Add legend info
                          jitter="label")
                #dashline.col="red")
                # title
                grid.text(gene, x=.5, y=.98, just="top", 
                          gp=gpar(cex=1.5, fontface="bold"))
                dev.off()
        }



        ############################################################
        # Hotspot distribution difference
        ############################################################
        all.consequences <- list("complex"=c("Complex"),
                                 "frameshift"=c("frameshift_variant", "indel", "inframe_codon_gain", "inframe_codon_loss", "inframe_variant", "initiator_codon_change", "intron_variant", "transcript_variant"),
                                 "missense"=c("non_synonymous_codon"),
                                 "nonsense"=c("stop_gained", "stop_lost"),
                                 "synonymous"=c("synonymous_codon"))

        summary.consequences <- c("Complex"="complex",
                                  "frameshift_variant"= "frameshift",
                                  "indel"= "frameshift",
                                  "inframe_codon_gain"= "frameshift",
                                  "inframe_codon_loss"= "frameshift",
                                  "inframe_variant"= "frameshift",
                                  "initiator_codon_change"= "frameshift",
                                  "intron_variant"= "frameshift",
                                  "transcript_variant"= "frameshift",
                                  "non_synonymous_codon"= "missense",
                                  "stop_gained"= "nonsense",
                                  "stop_lost"= "nonsense",
                                  "synonymous_codon"="synonymous")


        system(paste0("mkdir -p ", paste0("./results/analysis/gene/",gene)))

        # Compare hotspot frequency
        hotspots <- T
        cutoff <- 20

        dat.genetics.gene$PROTEIN_POS <- sapply(1:nrow(dat.genetics.gene), function(n)
                                                {
                                                        PROTEIN_NUM <- strsplit(dat.genetics.gene$PROTEIN_CHANGE[n],"[^0-9]+")[[1]][2]
                                                        return(as.numeric(PROTEIN_NUM))
                                                })


        dat.genetics.gene$type_summary <- sapply(1:nrow(dat.genetics.gene), function(n)
                                                 {
                                                         return(summary.consequences[which(names(summary.consequences)==dat.genetics.gene$CONSEQUENCE[n])])
                                                 })

        ##
        # TODO: Gather hotspots

        gene.hotspots <- names(table(dat.genetics.gene$PROTEIN_VARIANT))[table(dat.genetics.gene$PROTEIN_VARIANT)>=cutoff]

        if (length(grep("\\?", gene.hotspots))>0)
                gene.hotspots <- gene.hotspots[-grep("\\?", gene.hotspots)]

        if (length(grep("\\-", gene.hotspots))>0)
                gene.hotspots <- gene.hotspots[-grep("\\-", gene.hotspots)]

        # Remove NA
        if (length(grep("\\_NA", gene.hotspots))>0)
                gene.hotspots <- gene.hotspots[-grep("\\_NA", gene.hotspots)]

        if (length(gene.hotspots) > 0)
        {
                hotspots.info <- lapply(1:length(gene.hotspots), function(n)
                                        {
                                                return(list(begin=dat.genetics.gene$PROTEIN_POS[match(gene.hotspots[n], dat.genetics.gene$PROTEIN_VARIANT)], 
                                                            end=dat.genetics.gene$PROTEIN_POS[match(gene.hotspots[n], dat.genetics.gene$PROTEIN_VARIANT)],
                                                            type=NULL))
                                        })
                names(hotspots.info) <- gene.hotspots
        } else
        {
                hotspots.info <- list() 
        }

        ## Manual annotation
        if (gene=="CBL")
        {
                #gene.hotspots <- "CBL_380_410"
                hotspots.info <- list("CBL_380_410"=list(begin=380, end=410, type=NULL))
        } else if (gene=="EZH2")
        {
                hotspots.info <- list("EZH2_610_740_missense"=list(begin=610, end=740, type=NULL),
                                      "EZH2_nonsense"=list(begin=1, end=100000, type="nonsense"))
        } else if (gene=="GATA2")
        {
                hotspots.info <- list("GATA2_300_410"=list(begin=300, end=410, type=NULL))
        } else if (gene=="KRAS")
        {
                hotspots.info <- list("KRAS_12_13"=list(begin=12, end=13, type=NULL))
        } else if (gene=="NRAS")
        {
                gene.hotspots <- c("NRAS_12_13", "NRAS_61")
                hotspots.info <- list("NRAS_12_13"=list(begin=12, end=13, type=NULL),
                                      "NRAS_61"=list(begin=61, end=61, type=NULL))

        } else if (gene=="PTPN11")
        {
                gene.hotspots <- c("PTPN11_60_76", "PTPN11_480_510")
                hotspots.info <- list("PTPN11_60_76"=list(begin=60, end=76, type=NULL),
                                      "PTPN11_480_510"=list(begin=480, end=510, type=NULL))
        } else if (gene=="SF3B1")
        {
                gene.hotspots <- c("SF3B1_622_625", "SF3B1_662_666", "SF3B1_700")
                hotspots.info <- list("SF3B1_622_625"=list(begin=622, end=625, type=NULL),
                                      "SF3B1_662_666"=list(begin=662, end=666, type=NULL),
                                      "SF3B1_700"=list(begin=700, end=700, type=NULL))
        } else if (gene=="WT1")
        {
                gene.hotspots <- c("WT1_380_382", "WT1_462")
                hotspots.info <- list("WT1_380_382"=list(begin=380, end=382, type=NULL),
                                      "WT1_462"=list(begin=462, end=462, type=NULL))
        } else if (gene=="RUNX1")
        {
                gene.hotspots <- c("RUNX1_107", "RUNX1_162_166", "RUNX1_198_204")
                hotspots.info <- list("RUNX1_107"=list(begin=107, end=107, type=NULL),
                                      "RUNX1_162_166"=list(begin=162, end=166, type=NULL),
                                      "RUNX1_198_204"=list(begin=198, end=204, type=NULL))
        } else if (gene=="NPM1")
        {
                gene.hotspots <- c("NPM1_288")
                hotspots.info <- list("NPM1_288"=list(begin=287, end=288, type=NULL))
        }





        # If no hotspots then remove
        if (length(hotspots.info)==0)
                hotspots <- F

        if (!hotspots)
        {
                print("No hotspots")
                dat.genetics.gene$hotspot <- "non-hotspots"
        } else if (hotspots)
        {
                print("Hotspot detected")
                print("2. Working on clinical distribution")

                dat.clinical$gene <- sapply(1:nrow(dat.clinical), function(n)
                                            {
                                                    sample.info <- grep(dat.clinical$ID[n],dat.genetics.gene$ID)
                                                    if (length(sample.info)==0)
                                                    {
                                                            return("WT")
                                                    } else
                                                    {

                                                            hotspots.validity <- sapply(1:length(hotspots.info), function(j)
                                                                                        {
                                                                                                if (is.null(hotspots.info[[j]]$type))
                                                                                                {
                                                                                                        validity <- any(sapply(1:length(sample.info), function(k) { return(dat.genetics.gene$PROTEIN_POS[sample.info[k]]>= hotspots.info[[j]]$begin &dat.genetics.gene$PROTEIN_POS[sample.info[k]]<= hotspots.info[[j]]$end)}), na.rm=T)
                                                                                                } else
                                                                                                {
                                                                                                        validity <- any(sapply(1:length(sample.info), function(k) { return(dat.genetics.gene$PROTEIN_POS[sample.info[k]]>= hotspots.info[[j]]$begin &dat.genetics.gene$PROTEIN_POS[sample.info[k]]<= hotspots.info[[j]]$end  & dat.genetics.gene$CONSEQUENCE[sample.info[k]] %in% all.consequences[[hotspots.info[[j]]$type]] ) }), na.rm=T)
                                                                                                }
                                                                                        })

                                                            if (all(hotspots.validity==F))
                                                            {
                                                                    return("non-hotspots")
                                                            } else 
                                                            {
                                                                    return(names(hotspots.info)[(hotspots.validity)][1])
                                                            }
                                                    }
                                            })

                dat.genetics.gene$hotspot <- sapply(1:nrow(dat.genetics.gene), function(n)
                                                    {
                                                            hotspots.validity <- sapply(1:length(hotspots.info), function(k)
                                                                                        {
                                                                                                if (is.null(hotspots.info[[k]]$type))
                                                                                                {
                                                                                                        validity <- any(dat.genetics.gene$PROTEIN_POS[n]>= hotspots.info[[k]]$begin &dat.genetics.gene$PROTEIN_POS[n]<= hotspots.info[[k]]$end, na.rm=T)
                                                                                                } else
                                                                                                {
                                                                                                        validity <- any(dat.genetics.gene$PROTEIN_POS[n]>= hotspots.info[[k]]$begin &dat.genetics.gene$PROTEIN_POS[n]<= hotspots.info[[k]]$end  & dat.genetics.gene$CONSEQUENCE[n] %in% all.consequences[[hotspots.info[[k]]$type]] , na.rm=T)
                                                                                                }
                                                                                        })

                                                            if (all(hotspots.validity==F))
                                                            {
                                                                    return("non-hotspots")
                                                            } else 
                                                            {
                                                                    return(names(hotspots.info)[(hotspots.validity)][1])
                                                            }
                                                    })


                library(RColorBrewer)
                color.panel <- brewer.pal(length(unique(dat.clinical$gene)), "Set3")

                # BM_Blasts 
                dat.clinical$BM_Blasts <- as.numeric(as.character(dat.clinical$BM_Blasts))
                if (any(dat.clinical$BM_Blasts>100, na.rm=T))
                        dat.clinical$BM_Blasts[which(dat.clinical$BM_Blasts>100)] <- NA

                library(gridExtra)

                ##
                p.Age <- ggplot(dat.clinical) + geom_boxplot(aes(x=gene, y=Age, fill=disease)) + scale_fill_manual(values=colors.myeloid) + theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1) ) 
                p.WBC <- ggplot(dat.clinical) + geom_boxplot(aes(x=gene, y=WBC, fill=disease)) + scale_fill_manual(values=colors.myeloid) + theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.text.x= element_text(angle=45, hjust=1) ) 
                p.BM_Blasts <- ggplot(dat.clinical) + geom_boxplot(aes(x=gene, y=BM_Blasts, fill=disease)) + scale_fill_manual(values=colors.myeloid) + theme_bw() + theme(plot.title = element_text(hjust = 0.5) , axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1)) 

                library(ggpubr)
                pp.clinical <- ggarrange(p.Age, p.WBC, p.BM_Blasts, ncol=3, nrow=1, common.legend=T, legend='right')
                pdf(paste0("./results/analysis/gene/",gene,"/clinical_hotspots.pdf"), width=12, height=3)
                plot(pp.clinical)
                dev.off()

                ##################
                # Survival
                ##################
                print("3. Working on survival analysis")
                dat.clinical$OS[which(dat.clinical$OS <0)] <- - dat.clinical$OS[which(dat.clinical$OS <0)]
                dat.clinical$gene <- factor(dat.clinical$gene)
                colors.hotspot <- brewer.pal(length(levels(dat.clinical$gene)),"Set1")
                names(colors.hotspot) <- levels(dat.clinical$gene)

                source("./src/lib/analysis/run_univariate_survival.R")

                ############
                # Option 1 
                ############
                library(GGally)
                surv.disease.list <- lapply(1:length(unique(dat.clinical$disease)), function(n)
                                            {
                                                    disease <- as.character(unique(dat.clinical$disease)[n])
                                                    dat.clinical.disease <- dat.clinical[dat.clinical$disease==disease,]
                                                    pp.disease <- try(univariate_survival(dat.clinical.disease, feature="gene", cols=color.panel, out.dir=paste0("./results/analysis/gene/",gene,"/", disease, "/")),silent=T)
                                                    if(class(pp.disease)!="try-error")
                                                    {
                                                            # how many infos

                                                            #surv.disease <- try(ggsurv(pp.disease[[2]], surv.col=colors.hotspot[pp.disease[[1]]$xlevels[[1]]], order.legend=T) + ggtitle(disease) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + guides(linetype = FALSE) + scale_colour_discrete(name="Gene") + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()))
                                                            #surv.disease <- try(ggsurv(pp.disease[[2]], order.legend=F) + ggtitle(disease) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + guides(linetype = FALSE) + scale_colour_discrete(name="Gene", values=colors.hotspot) + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()))

                                                            surv.disease <- ggsurv(pp.disease[[2]], order.legend=F) + ggtitle(disease) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + guides(linetype = FALSE) + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_manual(name   = 'Gene',values=colors.hotspot)


                                                            if (any(class(surv.disease)=="try-error"))
                                                                    surv.disease <- NA
                                                    } else
                                                    {
                                                            surv.disease <- NA
                                                    }
                                                    return(surv.disease)
                                            })
                names(surv.disease.list) <- as.character(unique(dat.clinical$disease))

                if (!all(sapply(1:length(surv.disease.list), function(n){is.na(surv.disease.list[n])})))
                {
                        pp.surv <- ggarrange(surv.disease.list[[1]], surv.disease.list[[2]], surv.disease.list[[3]], ncol=length(surv.disease.list), nrow=1, common.legend=T, legend='right')

                        pdf(paste0("./results/analysis/gene/",gene,"/survival_hotspots.pdf"), width=6*length(surv.disease.list), height=3)
                        plot(pp.surv)
                        dev.off()
                }

                ######################################################
                # Hotspot interaction differences
                ######################################################
                print("4. Working on comutation pattern")

                source("./src/lib/plot/plot_compare_comutation.R")
                dat.hotspots.list <- lapply(1:length(unique(dat.clinical$disease)), function(n)
                                            {
                                                    disease <- as.character(unique(dat.clinical$disease)[n])
                                                    dat.clinical.disease <- dat.clinical[dat.clinical$disease==disease,]

                                                    dat.binary.hotspots <- data.frame(matrix(0, nrow=nrow(dat.clinical.disease), ncol=length(hotspots.info)))
                                                    colnames(dat.binary.hotspots) <- names(hotspots.info)
                                                    rownames(dat.binary.hotspots) <- dat.clinical.disease$ID

                                                    for (i in 1:length(hotspots.info))
                                                    {
                                                            samples.hotspots <- dat.clinical$ID[which(dat.clinical$gene==names(hotspots.info)[i])]
                                                            samples.hotspots <- samples.hotspots[samples.hotspots %in% rownames(dat.binary.hotspots)]
                                                            if (length(samples.hotspots) >0)
                                                                    dat.binary.hotspots[samples.hotspots,names(hotspots.info)[i]] <- 1
                                                    }

                                                    #if (length(hotspots.info)==1) # unique hotspot add other
                                                    #{
                                                    samples.non_hotspots <- dat.clinical$ID[which(dat.clinical$gene %in% "non-hotspots")]
                                                    samples.non_hotspots <- samples.non_hotspots[samples.non_hotspots %in% rownames(dat.binary.hotspots)]

                                                    dat.binary.hotspots <- cbind(dat.binary.hotspots, non_hotspots=0)
                                                    colnames(dat.binary.hotspots)[ncol(dat.binary.hotspots)] <- paste0(gene,"_others")
                                                    dat.binary.hotspots[samples.non_hotspots,ncol(dat.binary.hotspots)] <- 1
                                                    #}

                                                    return(dat.binary.hotspots)
                                            })
                names(dat.hotspots.list) <- unique(dat.clinical$disease)

                source("./src/lib/plot/plot_compare_groups.R")
                pp <- plot_compare_groups(dat.hotspots.list)

                pdf(paste0("./results/analysis/gene/",gene,"/proportion.pdf"), width=12, height=6)
                print(pp$p + scale_fill_manual(values=colors.myeloid))
                dev.off()

                dat.hotspots.mutated.only <- lapply(1:length(dat.hotspots.list), function(n)
                                                    {
                                                            return(dat.hotspots.list[[n]][ dat.genetics.gene$ID[dat.genetics.gene$ID %in% rownames(dat.hotspots.list[[n]])], ] )
                                                    })
                names(dat.hotspots.mutated.only) <- names(dat.hotspots.list)

                source("./src/lib/plot/plot_compare_groups.R")
                pp <- plot_compare_groups(dat.hotspots.mutated.only)

                pdf(paste0("./results/analysis/gene/",gene,"/proportion_mutated.pdf"), width=12, height=6)
                print(pp$p + scale_fill_manual(values=colors.myeloid))
                dev.off()

                source('./src/fun/merge_datasets.R')
                source("./src/lib/plot/plot_compare_comutation.R")

                ##########################################
                ##########################################
                disease.infos <- list(1:2,3:4,5)
                names(disease.infos) <- c("AML", "MDS", "MPN")

                # AML

                comut.list <- lapply(1:length(disease.infos), function(n)
                                     {
                                             disease <- names(disease.infos)[n]
                                             out.merge <- merge_datasets(dat.binary.list[disease.infos[[n]]], dat.cytogenetics.list[disease.infos[[n]]], dat.clinical.list[disease.infos[[n]]], features=NULL, impute=T, cutoff.occurrence=0.01)
                                             dat.integrate <- out.merge$dat
                                             pp <- try(plot_compare_comutation(dat.hotspots.list[[n]], dat.integrate, out.dir=paste0("./results/analysis/gene/",gene), suffix=paste0("_",disease,"_hotspots"), cache=F, frequent=F, plot=F), silent=T)

                                             if (class(pp)!="try-error")
                                             {
                                                     pp <- data.frame(pp)
                                                     pp.m <- melt(pp)
                                                     pp.m <- cbind(pp.m, rownames(pp))
                                                     colnames(pp.m) <- c("genes", "OR", disease)
                                                     #colnames(pp.m) <- c("genes", "OR", disease)

                                                     # if (length(unique(pp.m$OR))==1)
                                                     # {
                                                     #         pp.m$OR[1] <- pp.m$OR[1]+0.01
                                                     # }

                                                     label.OR <- c("1e-4","1e-3","1e-2","0.1","1","10","100","1000","1e4")
                                                     comut <- ggplot(pp.m) + geom_tile(aes_string(x="genes", y=disease, fill="OR") ,colour="grey80") + theme_minimal() + theme(axis.title.x=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), axis.ticks=element_blank(), plot.title=element_text(size=14,face="bold") ,  axis.text.x=element_text(angle=90,hjust=1)) + scale_fill_gradientn(colours=brewer.pal(9,"RdBu"), name="Odd Ratio", labels=label.OR, breaks=1:9, limits=c(1,9))
                                             } else
                                             {
                                                     comut <- NA
                                             }
                                             return(comut)
                                     })
                names(comut.list) <- names(disease.infos)

                pp.comut <- ggarrange(comut.list[[1]], comut.list[[2]], comut.list[[3]], ncol=1, nrow=length(comut.list), common.legend=T, legend='right')

                pdf(paste0("./results/analysis/gene/",gene,"/comutation_hotspots.pdf"), width=12, height=3*length(comut.list))
                plot(pp.comut)
                dev.off()


        }


        ########################################
        ## Look at type if table(missense/nonsense) >25%
        dat.genetics.gene$type_summary <- factor(dat.genetics.gene$type_summary, levels=c("frameshift","nonsense","complex","synonymous","missense"))
        type.summary <- table(dat.genetics.gene$type_summary)/nrow(dat.genetics.gene)



        if (sum(type.summary[c("frameshift","complex","nonsense")]>0.2)&(type.summary[c("missense")]>0.2))
        {
                type.comparison <- T
        } else
        {
                type.comparison <- F
        }


        if (nrow(dat.genetics.gene) <50)
                type.comparison <- F

        if(type.comparison)
        {
                print("comparing types..")

                dat.clinical$gene.type <- sapply(1:nrow(dat.clinical), function(n)
                                                 {
                                                         sample.info <- grep(dat.clinical$ID[n],dat.genetics.gene$ID)
                                                         if (length(sample.info)==0)
                                                         {
                                                                 return("WT")
                                                         } else
                                                         {

                                                                 if (any(dat.genetics.gene$type_summary[sample.info] %in% c("frameshift","nonsense","complex")))
                                                                 {
                                                                         if (any(dat.genetics.gene$type_summary[sample.info] %in% c("missense")))
                                                                         {
                                                                                 return("both")
                                                                         } else
                                                                         {
                                                                                 return("nonsense")
                                                                         }
                                                                 } else 
                                                                 {
                                                                         return("missense")
                                                                 }
                                                         }
                                                 })

                library(RColorBrewer)
                color.panel <- brewer.pal(length(unique(dat.clinical$gene.type)), "Set3")

                # BM_Blasts 
                dat.clinical$BM_Blasts <- as.numeric(as.character(dat.clinical$BM_Blasts))
                if (any(dat.clinical$BM_Blasts>100,na.rm=T))
                        dat.clinical$BM_Blasts[which(dat.clinical$BM_Blasts>100)] <- NA

                library(gridExtra)

                ##
                p.Age <- ggplot(dat.clinical) + geom_boxplot(aes(x=gene.type, y=Age, fill=disease)) + scale_fill_manual(values=colors.myeloid) + theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1) ) 
                p.WBC <- ggplot(dat.clinical) + geom_boxplot(aes(x=gene.type, y=WBC, fill=disease)) + scale_fill_manual(values=colors.myeloid) + theme_bw() + theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_blank(), axis.text.x= element_text(angle=45, hjust=1) ) 
                p.BM_Blasts <- ggplot(dat.clinical) + geom_boxplot(aes(x=gene.type, y=BM_Blasts, fill=disease)) + scale_fill_manual(values=colors.myeloid) + theme_bw() + theme(plot.title = element_text(hjust = 0.5) , axis.title.x=element_blank(), axis.text.x=element_text(angle=45, hjust=1)) 

                library(ggpubr)
                pp.clinical <- ggarrange(p.Age, p.WBC, p.BM_Blasts, ncol=3, nrow=1, common.legend=T, legend='right')
                pdf(paste0("./results/analysis/gene/",gene,"/clinical_hotspots_type.pdf"), width=12, height=3)
                plot(pp.clinical)
                dev.off()

                ##################
                # Survival
                ##################
                print("3. Working on survival analysis")
                dat.clinical$OS[which(dat.clinical$OS <0)] <- - dat.clinical$OS[which(dat.clinical$OS <0)]
                dat.clinical$gene.type <- factor(dat.clinical$gene.type)
                colors.hotspot <- brewer.pal(length(levels(dat.clinical$gene.type)),"Set1")
                names(colors.hotspot) <- levels(dat.clinical$gene.type)

                source("./src/lib/analysis/run_univariate_survival.R")

                ############
                # Option 1 
                ############
                library(GGally)
                surv.disease.list <- lapply(1:length(unique(dat.clinical$disease)), function(n)
                                            {
                                                    disease <- as.character(unique(dat.clinical$disease)[n])
                                                    dat.clinical.disease <- dat.clinical[dat.clinical$disease==disease,]
                                                    pp.disease <- try(univariate_survival(dat.clinical.disease, feature="gene.type", cols=color.panel, out.dir=paste0("./results/analysis/gene/",gene,"/", disease, "/")),silent=T)
                                                    if(class(pp.disease)!="try-error")
                                                    {
                                                            surv.disease <- ggsurv(pp.disease[[2]], order.legend=F) + ggtitle(disease) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + guides(linetype = FALSE) + theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank()) + ggplot2::guides(linetype = FALSE) + ggplot2::scale_colour_manual(name   = 'Gene',values=colors.hotspot)
                                                            if (any(class(surv.disease)=="try-error"))
                                                                    surv.disease <- NA
                                                    } else
                                                    {
                                                            surv.disease <- NA
                                                    }
                                                    return(surv.disease)
                                            })
                names(surv.disease.list) <- as.character(unique(dat.clinical$disease))

                if (!all(sapply(1:length(surv.disease.list), function(n){is.na(surv.disease.list[n])})))
                {
                        pp.surv <- ggarrange(surv.disease.list[[1]], surv.disease.list[[2]], surv.disease.list[[3]], ncol=length(surv.disease.list), nrow=1, common.legend=T, legend='right')

                        pdf(paste0("./results/analysis/gene/",gene,"/survival_hotspots_type.pdf"), width=6*length(surv.disease.list), height=3)
                        plot(pp.surv)
                        dev.off()
                }

                ######################################################
                # Hotspot interaction differences
                ######################################################
                print("4. Working on comutation pattern")

                source("./src/lib/plot/plot_compare_comutation.R")
                dat.hotspots.list <- lapply(1:length(unique(dat.clinical$disease)), function(n)
                                            {
                                                    disease <- as.character(unique(dat.clinical$disease)[n])
                                                    dat.clinical.disease <- dat.clinical[dat.clinical$disease==disease,]

                                                    dat.binary.hotspots <- data.frame(matrix(0, nrow=nrow(dat.clinical.disease), ncol=2))
                                                    colnames(dat.binary.hotspots) <- paste0(gene,c("_missense","_nonsense"))
                                                    rownames(dat.binary.hotspots) <- dat.clinical.disease$ID

                                                    samples.hotspots <- dat.clinical$ID[which(dat.clinical$gene.type=="missense")]
                                                    samples.hotspots <- samples.hotspots[samples.hotspots %in% rownames(dat.binary.hotspots)]
                                                    if (length(samples.hotspots) >0)
                                                            dat.binary.hotspots[samples.hotspots,1] <- 1

                                                    samples.hotspots <- dat.clinical$ID[which(dat.clinical$gene.type %in% c("complex","frameshift", "nonsense"))]
                                                    samples.hotspots <- samples.hotspots[samples.hotspots %in% rownames(dat.binary.hotspots)]
                                                    if (length(samples.hotspots) >0)
                                                            dat.binary.hotspots[samples.hotspots,2] <- 1


                                                    return(dat.binary.hotspots)
                                            })
                names(dat.hotspots.list) <- unique(dat.clinical$disease)

                source("./src/lib/plot/plot_compare_groups.R")
                pp <- plot_compare_groups(dat.hotspots.list)

                pdf(paste0("./results/analysis/gene/",gene,"/proportion_type.pdf"), width=12, height=6)
                print(pp$p + scale_fill_manual(values=colors.myeloid))
                dev.off()

                dat.hotspots.mutated.only <- lapply(1:length(dat.hotspots.list), function(n)
                                                    {
                                                            return(dat.hotspots.list[[n]][ dat.genetics.gene$ID[dat.genetics.gene$ID %in% rownames(dat.hotspots.list[[n]])], ] )
                                                    })
                names(dat.hotspots.mutated.only) <- names(dat.hotspots.list)

                source("./src/lib/plot/plot_compare_groups.R")
                pp <- plot_compare_groups(dat.hotspots.mutated.only)

                pdf(paste0("./results/analysis/gene/",gene,"/proportion_mutated_type.pdf"), width=12, height=6)
                print(pp$p + scale_fill_manual(values=colors.myeloid))
                dev.off()

                source('./src/fun/merge_datasets.R')
                source("./src/lib/plot/plot_compare_comutation.R")

                ##########################################
                ##########################################
                disease.infos <- list(1:2,3:4,5)
                names(disease.infos) <- c("AML", "MDS", "MPN")

                # AML

                comut.list <- lapply(1:length(disease.infos), function(n)
                                     {
                                             disease <- names(disease.infos)[n]
                                             out.merge <- merge_datasets(dat.binary.list[disease.infos[[n]]], dat.cytogenetics.list[disease.infos[[n]]], dat.clinical.list[disease.infos[[n]]], features=NULL, impute=T, cutoff.occurrence=0.01)
                                             dat.integrate <- out.merge$dat
                                             pp <- try(plot_compare_comutation(dat.hotspots.list[[n]], dat.integrate, out.dir=paste0("./results/analysis/gene/",gene), suffix=paste0("_",disease,"_hotspots"), cache=F, frequent=F, plot=F), silent=T)

                                             if (class(pp)!="try-error")
                                             {
                                                     pp <- data.frame(pp)
                                                     pp.m <- melt(pp)
                                                     pp.m <- cbind(pp.m, rownames(pp))
                                                     colnames(pp.m) <- c("genes", "OR", disease)
                                                     #colnames(pp.m) <- c("genes", "OR", disease)

                                                     label.OR <- c("1e-4","1e-3","1e-2","0.1","1","10","100","1000","1e4")
                                                     comut <- ggplot(pp.m) + geom_tile(aes_string(x="genes", y=disease, fill="OR") ,colour="grey80") + theme_minimal() + theme(axis.title.x=element_blank(), panel.border=element_blank(), panel.grid=element_blank(), axis.ticks=element_blank(), plot.title=element_text(size=14,face="bold") ,  axis.text.x=element_text(angle=90,hjust=1)) + scale_fill_gradientn(colours=brewer.pal(9,"RdBu"), name="Odd Ratio", labels=label.OR, breaks=1:9, limits=c(1,9))
                                             } else
                                             {
                                                     comut <- NA 
                                             }
                                             return(comut)
                                     })
                names(comut.list) <- names(disease.infos)

                pp.comut <- ggarrange(comut.list[[1]], comut.list[[2]], comut.list[[3]], ncol=1, nrow=length(comut.list), common.legend=T, legend='right')

                pdf(paste0("./results/analysis/gene/",gene,"/comutation_hotspots_type.pdf"), width=12, height=3*length(comut.list))
                plot(pp.comut)
                dev.off()
        }

        write.table(dat.genetics.gene, file=paste0("./results/analysis/gene/",gene,"_genetics.txt"), sep="\t", quote=F, row.names=F)

}

dat.genetics.final <- Reduce("rbind", lapply(1:length(gene.list), function(n)
                                             { return(read.table(paste0("./results/analysis/gene/",gene.list[n],"_genetics.txt"), sep="\t", header=T, stringsAsFactors=F))}))

write.table(dat.genetics.final, file="~/Desktop/genetics_hotspots.txt", quote=F, sep="\t")


