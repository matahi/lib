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


        if (length(p.lollipop)==3)
        {
         p.final <- ggarrange(
                              p.lollipop[[1]], 
                              p.lollipop[[2]], 
                              p.lollipop[[3]],
                              p.protein,
                              ncol=1, 
                              nrow=4, 
                              common.legend=F, 
                              heights=c(4,4,4,1.5), 
                              align="v")

        } else if (length(p.lollipop)==2)
        {

         p.final <- ggarrange(
                              p.lollipop[[1]], 
                              p.lollipop[[2]], 
                              p.protein,
                              ncol=1, 
                              nrow=3, 
                              common.legend=F, 
                              heights=c(4,4,1.5), 
                              align="v")

        } else if (length(p.lollipop)==1)
        {

         p.final <- ggarrange(
                              p.lollipop[[1]], 
                              p.protein,
                              ncol=1, 
                              nrow=2, 
                              common.legend=F, 
                              heights=c(4,1.5), 
                              align="v")
        }

        if (plot.title)
                p.final <- annotate_figure(p.final,
                                   top = text_grob(current.gene, color="black", face="bold", size=size.text))

        return(p.final)
}


