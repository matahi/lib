plot_confidence_CI <- function(mutation.df, sample, out.dir="./")
{
        require(ggplot2)

        mutation.df.sample <- mutation.df[(mutation.df$ID==sample),]

        # pdf(paste0(out.dir, "/", sample,".pdf"))
        p <-  ggplot(mutation.df.sample) + geom_point(aes(x=Gene, y=VAF.corr),size=5) + ylim(0,1) + 
                geom_segment(aes(x=Gene, xend=Gene, y=LCI.corr, yend=UCI.corr)) + theme_bw() +
                        theme(legend.position="none")
                #scale_colour_manual(values=col.gene) + theme(legend.position="none")
                # dev.off()

              return(p)


}
