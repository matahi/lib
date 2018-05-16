# wrapper
# Main function is plot_lollipop
filter_residue <- function (residueSeq) {
        if (any(residueSeq %in% c("OPAL", "OCHRE", "AMBER"))) {
                stopRes <- c("OPAL", "OCHRE", "AMBER")
                residueSeq <- residueSeq[-which(residueSeq %in% stopRes)]
        }

        return(residueSeq)
        
}

get_protein_info <- function (gene) {
        require(GenVisR)

        # TODO Recover good version
        ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset= "hsapiens_gene_ensembl", ensemblRedirect = F)
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
                residueSeq <- filter_residue(residueSeq)
                # Get length
                protein_length <- length(residueSeq)
                # Recover domains
                protein_domain <- getBM(attributes=c("interpro_description", "interpro_start","interpro_end"), filters=c("ensembl_transcript_id"), values=transcript.taken, mart=ensembl)
                if(all(is.na(protein_domain$interpro_description)))
                        protein_domain <- data.frame(interpro_description=gene, interpro_start=1, interpro_end=proteinLength)
        } else
        {
                return(NULL)
        }

        return(list(protein_length=protein_length,protein_domain=protein_domain))
        
}

my_lollipop_facet <- function (dat.genetics.unique, 
                               current.gene, 
                               facet_var=NULL,
                               protein_domain, 
                               protein_length, 
                               cutoff.hotspot=20, 
                               plot.hotspot = F,
                               colors.mutation_type = NULL,
                               plot.domains = "automatic",
                               size.text = 20,
                               plot.title=T
                               ) {

        require(RColorBrewer)
        require(ggthemr)
        ggthemr("fresh")
        
        
        ## 1. Protein
        protein_domain <- protein_domain %>% mutate(width = interpro_end  - interpro_start + 1,
                                                    protein_length = protein_length) %>%
                                             group_by(interpro_description) %>% 
                                             mutate(height = 1 / (10*n())) %>% 
                                             ungroup()

        # TODO: FIX COLOR
        if (plot.domains == "automatic")
        {
                plot.domains <- T
                n.domains <- length(unique(protein_domain$interpro_description))
                if (n.domains > 12)
                        plot.domains <- F
        }

        y.limits <- max(c(protein_domain$height *2, 0.1))

        if (plot.domains)
        {
                color.domains <- brewer.pal(max(3,n.domains), "Set3")[1:n.domains]


                p.protein <- ggplot(protein_domain) + 
                        geom_segment(aes(y=0,yend=0,x=1,xend=protein_length)) +  # 1. first plot the full protein
                                geom_rect(aes(xmin=interpro_start,
                                              xmax=interpro_end,
                                              ymin=-height,
                                              ymax=height,
                                              fill=interpro_description),
                                              alpha= 0.8) + # 2. plot domains
                        scale_fill_manual(name = "", values=color.domains) +
                        theme_void() +
                        theme(legend.position = "none") +
                        # theme(legend.position = "bottom") +
                        # guides(fill = guide_legend(ncol = 2)) +
                        xlim(0, protein_length+10) +
                        ylim(-y.limits, y.limits)
        } else
        {
                p.protein <- ggplot(protein_domain) + 
                        geom_segment(aes(y=0,yend=0,x=1,xend=protein_length)) +  # 1. first plot the full protein
                                geom_rect(aes(xmin=interpro_start,
                                              xmax=interpro_end,
                                              ymin=-height,
                                              ymax=height),
                                              alpha= 0.8) + # 2. plot domains
                        theme_void() +
                        xlim(0, protein_length+10) +
                        ylim(-y.limits, y.limits)
        }

        p.protein <- p.protein + theme(text = element_text(size=size.text))

        ## 2. lollipop
        dat.genetics.gene <- dat.genetics.unique %>% filter(Gene == current.gene, !is.na(PROTEIN_POS))

        dat.score <- dat.genetics.gene %>% 
                group_by_(facet_var, "PROTEIN_CHANGE")%>% 
                summarise(Gene = unique(Gene),
                          score = n(),
                          PROTEIN_POS = unique(PROTEIN_POS),
                          mutation_type = mutation_type[1]) %>% # just take 1
                          # mutation_type = unique(mutation_type)) %>%
                ungroup()

        colnames(dat.score) <- c("facet.var", colnames(dat.score)[2:ncol(dat.score)])

        all.mutation_types <- unique(dat.score$mutation_type)

        facets <- unique(dat.score$facet.var)

        # We add a fake mutation with NA to fix color scale
        dat.score <- rbind(dat.score,
                           data.frame(facet.var = facets[1],
                                      PROTEIN_CHANGE = NA,
                                      Gene = current.gene,
                                      score = NA,
                                      PROTEIN_POS = NA,
                                      mutation_type = all.mutation_types))
        # 
        p.lollipop <- lapply(1:length(facets), function(k)
                             {
                                     # plotting legend
                                     legend.info <- ifelse(k==1,"top","none")

                                     # plot
                                     tmp <- ggplot(dat.score %>% filter(facet.var == facets[k])) + 
                                         geom_segment(aes(x=PROTEIN_POS,
                                                          xend=PROTEIN_POS,
                                                          y=0,
                                                          yend=score),
                                                      alpha=0.5, na.rm=T) +
                                         geom_point(aes(x=PROTEIN_POS,
                                                        y=score,
                                                        color=mutation_type),
                                                        alpha=0.7, na.rm=T) + 
                                         ylab(facets[k]) + xlab("") + xlim(0, protein_length + 10) +
                                         theme(legend.position=legend.info)

                                 
                                 if (!is.null(colors.mutation_type))
                                         tmp <- tmp + scale_colour_manual(name = "", values=colors.mutation_type)

                                 ### modify size
                                 tmp <- tmp + theme(text = element_text(size=size.text))

                                 if (plot.hotspot)
                                         tmp <- tmp +  geom_text(data=dat.score %>% filter(facet.var == facets[k],
                                                                                           score >= cutoff.hotspot),
                                                                 aes(label = PROTEIN_CHANGE,
                                                                     x=PROTEIN_POS + 100,
                                                                     y= score ),
                                                                 position=position_jitter(width=5,height=2)) # TODO: this depends on the protein_length size and the max score



                                 return(tmp)
                             })


        p.final <- ggarrange(plotlist = p.lollipop,
                             p.protein,
                             ncol=1,
                             nrow= length(p.lollipop) +1,
                             common.legend = F,
                             heights =c(rep(4, length(p.lollipop)),1.5),
                             align="v")


        if (plot.title)
                p.final <- annotate_figure(p.final,
                                   top = text_grob(current.gene, color="black", face="bold", size=size.text))

        return(p.final)
}


plot_lollipop <- function (gene, dat.genetics, ...) {

        # gene
        # dat.genetics: VCF format with columns: Gene, PROTEIN_POS, PROTEIN_CHANGE, mutation_type (sorry for the inconsistency)

        # a. Recover protein
        prot.info <- get_protein_info(gene)

        # If no protein info then we cannot plot the gene
        if (is.null(prot.info)) 
                stop("No protein information for that gene")

        # b. plot info 
        pp <- my_lollipop_facet(dat.genetics, 
                                current.gene= gene, 
                                protein_domain=prot.info$protein_domain, 
                                protein_length=prot.info$protein_length, ...)

        # c. plot
        return(pp)
        # plot(pp)
}



