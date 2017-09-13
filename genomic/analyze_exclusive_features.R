analyze_exclusive_features <- function(dat.list, Dataset.list)
{
        # dat.list <- karyotypes

        exclusive.features <- lapply(1:length(dat.list), function(n)
                                     {
                                             features <- setdiff(colnames(dat.list[[n]]), Reduce('c', lapply(dat.list[-n],colnames)))
                                     })



        dat.processed <- lapply(1:length(dat.common), function(n)
                                {
                                        n_event <- colSums(dat.list[[n]][,exclusive.features[[n]]])

                                        index <- order(n_event,decreasing=T)

                                        tmp.df <- data.frame(features=names(n_event)[index], n_event=n_event[index], frequency=n_event[index]/nrow(dat.list[[n]])*100)
                                        

                                        return(tmp.df)
                                })

        for (n in 1:length(dat.processed))
        {
                write.table(dat.processed[[n]], file=paste0("~/Documents/Work/Postdoc/projects/Myeloid/results/cytogenetics/events/",Dataset.list[[n]],"_exclusive_events.txt"), quote=F, sep="\t", row.names=F)
        }

        return(dat.processed)
}
