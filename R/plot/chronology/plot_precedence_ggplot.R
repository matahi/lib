plot_precedence <- function(chronology.df,
                            total.range = c(-1,0,5,10,20,50,100,200),
                            Reordering=NULL)
{

        # Loading
        require(RColorBrewer)

        ########################
        # TESTING
        ########################

        # Reordering <- NULL
        # Reordering is a df with 3 columns:
        # - alterations
        # - groups
        # - colours

        # Reordering <- Reordering.CH_SP

        ## Gene ordering
        if (!is.null(Reordering))
        {
                Gene.levels <- levels(Reordering$alterations)[levels(Reordering$alterations) %in% unique(chronology.df$Gene1)]
        } else
        {
                Gene.levels <- unique(chronology.df$Gene1)
        }

        #
        precedence.df <- chronology.df %>% 
                mutate(
                       Gene1 = factor(Gene1, levels=Gene.levels),
                       Gene2 = factor(Gene2, levels=Gene.levels),
                       Ratio.win = ifelse(Win + Loss > 0, Win / (Win + Loss), NA),
                       Ratio.inconclusive = ifelse( Total > 0 , 1 - (Win + Loss) / Total, NA))

        ########################
        # TODO: TESTING default values
        ########################

        ########################
        # CALCULATING OR/Pvals/Pairs
        ########################

        # ###
        # logPInt 
        # odds
        # pairs
        # ###

        # logOdds == Ratio.win
        # logPInt == Ratio.inconclusive
        # pairs == Total 

        ###########################
        ## FORMATTING
        ###########################

        Ratio.inconclusive.df <- precedence.df %>% select(Gene1,Gene2,Ratio.inconclusive)
        Ratio.win.df <- precedence.df %>% select(Gene1,Gene2,Ratio.win)
        Total.df <- precedence.df %>% select(Gene1,Gene2,Total)

        ## Working on range
        # i. Breaks for Total 
        # TODO: improve total range and labelling
        max.val <- Total.df %>% .$Total %>% max(., na.rm=T)
        if ((max.val + 100)> max(total.range))
                total.range <- c(total.range, (max.val + 100))
        Total.labels <- formatC(total.range[-1])

        # ii. Breaks for ratio
        ratio.breaks <- c(-.Machine$double.eps, seq(0,1,0.25))
        ratio.labels <- scales::percent(seq(0,1,0.25)) 

        ## Formatting data.frames
        Total.upper <- Total.df %>% 
                filter(as.numeric(Gene1) > as.numeric(Gene2)) %>%
                #filter(as.numeric(Gene1) >= as.numeric(Gene2)) %>%
                mutate(
                       value=cut(Total, label=Total.labels, breaks=total.range))

        Ratio.win.lower <- Ratio.win.df %>% 
                filter(as.numeric(Gene1) < as.numeric(Gene2)) %>%
                mutate(
                       value=cut(Ratio.win, breaks = ratio.breaks, label= ratio.labels, include.lowest=TRUE))

        # Add dummy values
        Ratio.win.lower <- bind_rows( Ratio.win.lower, data.frame(Gene1=NA, Gene2=NA, value=levels(Ratio.win.lower$value)))
        Total.upper <- bind_rows( Total.upper, data.frame(Gene1=NA, Gene2=NA, value=levels(Total.upper$value)))

        # pairs colorscale
        colors.Total <- brewer.pal(length(Total.labels), "Greens")
        names(colors.Total) <- Total.labels

        # odds ratio colorscale
        colors.ratio <- c(brewer.pal(length(ratio.labels)-1, "RdBu"),"white")
        names(colors.ratio) <- c(setdiff(ratio.labels,"50%"), "50%") # TODO: rework with formatting
        colors.ratio <- colors.ratio[ratio.labels] # Reorder


        ## Plot
        pp <- ggplot() + 
                geom_tile(data=Total.upper, aes(x=Gene1, y=Gene2, fill=value), colour="white", show.legend=F) +
                geom_point(data=Total.upper, aes(x=NA, y=NA, color=value)) + # dummy plots for legend purposes
                geom_tile(data=Ratio.win.lower, aes(x=Gene1, y=Gene2, fill=value), colour="white", show.legend=F) +
                geom_point(data=Ratio.win.lower, aes(x=NA, y=NA, alpha=value)) + # dummy plots for legend purposes
                scale_fill_manual(values=c(colors.Total, colors.ratio), na.value="grey50") + 
                scale_x_discrete(limits=Gene.levels) +
                scale_y_discrete(limits=Gene.levels) +
                theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
                xlab("") + ylab("") +
                guides(color = guide_legend(title = "Co-mutation", order = 1, 
                                            override.aes = list(shape = 15, size = 5, color = colors.Total) ),
                       alpha = guide_legend(title = "Win Ratio", order = 2, 
                                            override.aes = list(shape = 15, size = 5, color = c(colors.ratio,"grey50"), alpha=1) )) +
                theme(legend.key = element_rect(fill = "white") )

        ## TODO: I AM HERE

      if (!is.null(Reordering))
      {
              cols.axis <- Reordering$colours[ match(Gene.levels, Reordering$alterations )] 
              features.sep <- Reordering %>% count(groups) %>% mutate(pos = cumsum(n)) %>% slice(1:(n()-1))
              pp <- pp +
                geom_hline(data= features.sep, aes(yintercept=pos+0.5), col="black") + # features separation
                geom_vline(data= features.sep, aes(xintercept=pos+0.5), col="black") + # features separation
                theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, colour= cols.axis),
                      axis.text.y = element_text(colour= cols.axis))
      }

        return(pp)

}


