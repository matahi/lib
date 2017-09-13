plot_compare_burden <- function(dat.binary.list, dat.cytogenetics.list, features="all", out.dir="ALL")
{
        require(ggplot2)
        require(scales)
        require(dplyr)
        require(reshape2)

        # Number of mutations per samples
        dat.processed <- Reduce("rbind",lapply(1:length(dat.binary.list), function(n)
                                               {
                                                       data.frame(aberrations = c(rowSums(dat.binary.list[[n]],na.rm=T), rowSums(dat.binary.list[[n]],na.rm=T)) , 
                                                                  type= rep(c("Gene","CNAs"), each=nrow(dat.binary.list[[n]])),
                                                                  dataset=names(dat.binary.list)[n])
                                               }))

        p <- ggplot(dat.processed) + geom_histogram(aes(x=aberrations, fill=dataset), colour="grey", binwidth=1) + scale_fill_manual(values=colors.datasets) + facet_grid(dataset ~ type,scales="free_y") + theme_bw()

        pdf(paste0("results/comparison/", out.dir,"/genetic_burden.pdf"))
        print(p)
        dev.off()


        if (length(dat.binary.list)==2)
        {

        }





        return(p)
}
