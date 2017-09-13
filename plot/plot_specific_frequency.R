plot_specific_frequency <- function(dat.list,features=NULL, add.significance=F, cutoff.occurrence= cutoff.occurrence, out.dir="ALL")
{
        # dat.list <- dat.binary.list

        require(ggplot2)
        require(scales)
        require(dplyr)
        require(reshape2)

        features.list <- lapply(dat.list,colnames)

        for (k in 1:length(features.list))
        {
                print(k)

                # Create the data.frame
                n_event <- colSums(dat.list[[k]],na.rm=T)
                dat.processed <- data.frame(features=names(n_event), n_event, type="error",  dataset=names(dat.list)[k], stringsAsFactors=F)


                # characterize features intersection
                common.features <-  Reduce("intersect",features.list)
                unique.features <-  setdiff(features.list[[k]], Reduce("union",features.list[-k]))
                type <- substr(names(dat.list)[k],1,3)
                type.index <- grep(type, names(dat.list))
                #type.features <-  setdiff(Reduce("intersect", features.list[type.index] ), Reduce("union",features.list[-type.index]))
                type.features <-  setdiff(Reduce("intersect", features.list[type.index]), common.features)
                intersect.features <-  setdiff( intersect(features.list[[k]], Reduce("union",features.list[-type.index])), union(common.features,type.features) )

                dat.processed$type[dat.processed$features %in% unique.features] <- "unique"
                dat.processed$type[dat.processed$features %in% common.features] <- "common"
                dat.processed$type[dat.processed$features %in% intersect.features] <- "intersect"
                dat.processed$type[dat.processed$features %in% type.features] <- "type"

                # Order
                order_var <- order(dat.processed$n_event,decreasing=T)
                dat.processed$features <- factor(dat.processed$features,levels=dat.processed$features[order_var])

                #p.specific <- ggplot() + geom_bar(dat=dat.processed,aes(x=features,y=n_event, fill=dataset), stat="identity", position="identity") + geom_hline(yintercept=cutoff.occurrence, linetype="longdash") +
                # scale_fill_manual(values=colors.datasets) +
                p.specific <- ggplot() + geom_bar(dat=dat.processed,aes(x=features,y=n_event, fill=type), stat="identity", position="identity") + # geom_hline(yintercept=cutoff.occurrence, linetype="longdash") +
                scale_fill_manual(values=colors.mut_type) +
                theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top") 

                pdf(paste0("results/comparison/",out.dir,"/specific_frequency_", names(dat.list)[k],".pdf"), width=16, height=9)
                print(p.specific + geom_hline(yintercept=10, color="black", linetype="longdash"))
                dev.off()

        }

        return(NA)
}
