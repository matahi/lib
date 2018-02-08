

my_lollipop <- function (dat.genetics.unique, current.gene, diseases=NULL,protein_domain, protein_length, cutoff.hotspot=20, colors.mutation_type = NULL) {

        require(ggthemr)
        ggthemr("fresh")
        
        
        ## 1. Protein
        protein_domain <- protein_domain %>% mutate(width = interpro_end  - interpro_start + 1,
                                                    protein_length = protein_length) %>%
                                             group_by(interpro_description) %>% 
                                             mutate(height = 1 / (10*n())) %>% 
                                             ungroup()

        # TODO: FIX COLOR
        p.protein <- ggplot(protein_domain) + 
                geom_segment(aes(y=0,yend=0,x=1,xend=protein_length)) +  # 1. first plot the full protein
                        geom_rect(aes(xmin=interpro_start,
                                      xmax=interpro_end,
                                      ymin=-height,
                                      ymax=height),
                                      #fill=interpro_description),
                                      alpha= 0.8) + # 2. plot domains
                theme_void() +
                theme(legend.position = "bottom") +
                xlim(0, protein_length+10)

        ## 2. lollipop
        dat.genetics.gene <- dat.genetics.unique %>% filter(Gene == current.gene, !is.na(PROTEIN_POS))

        dat.score <- dat.genetics.gene %>% 
                group_by(DISEASE, PROTEIN_CHANGE)%>% 
                summarise(Gene = unique(Gene),
                          score = n(),
                          PROTEIN_POS = unique(PROTEIN_POS),
                          mutation_type = mutation_type[1]) %>% # just take 1
                          # mutation_type = unique(mutation_type)) %>%
                ungroup()

        # ## this is for testing
        # testing <- dat.genetics.gene %>% 
        #         group_by(DISEASE, PROTEIN_CHANGE)%>% 
        #         summarise(Gene = unique(Gene),
        #                   score = n(),
        #                   PROTEIN_POS = unique(PROTEIN_POS),
        #                   mutation_type = length(unique(mutation_type))) %>%
        #         ungroup()

        # weird.change <- testing  %>% filter(mutation_type==2) %>% .$PROTEIN_CHANGE
        # dat.genetics.gene %>% filter(PROTEIN_CHANGE == weird.change ) %>% select(mutation_type, dataset)

        if (is.null(diseases))
                diseases <- as.character(na.omit(unique(dat.score$DISEASE)))

        all.mutation_types <- unique(dat.score$mutation_type)

        # we add a fake row for legend's sake
        dat.score <- rbind(dat.score,
                           data.frame(DISEASE = diseases[1],
                                      PROTEIN_CHANGE = NA,
                                      Gene = current.gene,
                                      score = NA,
                                      PROTEIN_POS = NA,
                                      mutation_type = all.mutation_types))

        # 
        p.lollipop <- lapply(1:length(diseases), function(k)
                             {
                                     # plotting legend
                                     legend.info <- ifelse(k==1,"top","none")

                                     # plot
                                     tmp <- ggplot(dat.score %>% filter(DISEASE == diseases[k])) + 
                                         geom_segment(aes(x=PROTEIN_POS,
                                                          xend=PROTEIN_POS,
                                                          y=0,
                                                          yend=score),
                                                      alpha=0.5, na.rm=T) +
                                         geom_point(aes(x=PROTEIN_POS,
                                                        y=score,
                                                        color=mutation_type),
                                                        alpha=0.7, na.rm=T) + 
                                         geom_text(data=dat.score %>% filter(DISEASE == diseases[k],
                                                                             score >= cutoff.hotspot),
                                                         aes(label = PROTEIN_CHANGE,
                                                             x=PROTEIN_POS + 100,
                                                             y= score ),
                                                         position=position_jitter(width=5,height=2)) + # TODO: this depends on the protein_length size and the max score
                                         ylab(diseases[k]) + xlab("") + xlim(0, protein_length + 10) +
                                         theme(legend.position=legend.info)


                                 if (!is.null(colors.mutation_type))
                                         tmp <- tmp + scale_colour_manual(values=colors.mutation_type)

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
                              heights=c(4,4,4,1), 
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
                              heights=c(4,4,1), 
                              align="v")

        } else if (length(p.lollipop)==1)
        {

         p.final <- ggarrange(
                              p.lollipop[[1]], 
                              p.protein,
                              ncol=1, 
                              nrow=2, 
                              common.legend=F, 
                              heights=c(4,1), 
                              align="v")
        }

        p.final <- annotate_figure(p.final,
                                   top = text_grob(current.gene, color="black", face="bold", size=14))

        return(p.final)
}

# p.final <- ggarrange(
#                      p.protein,
#                      p.lollipop[[1]], 
#                      p.lollipop[[2]], 
#                      p.lollipop[[3]],
#                      ncol=1, 
#                      nrow=4, 
#                      common.legend=F, 
#                      heights=c(1,4,4,4), 
#                      align="v")


