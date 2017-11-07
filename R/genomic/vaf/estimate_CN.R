estimate_CN <- function (dat.cytogenetics, gene.list, gene.infos)
{
        # dat.cytogenetics <- dat.cytogenetics.list[[1]]

        ##
        # Default at 2 copies
        ##
        CN.df <- data.frame(matrix(2,nrow=nrow(dat.cytogenetics), ncol=length(gene.list)))
        rownames(CN.df) <- rownames(dat.cytogenetics)
        colnames(CN.df) <- gene.list

        for (k in 1:length(gene.list))
        {
                chr <- gene.infos[k,"Chr"]
                band <- substr(gene.infos[k,"Band"],1,1)

                # amplification
                feature.plus <- paste0("plus", chr)
                if (feature.plus %in% colnames(dat.cytogenetics))
                        CN.df[ which(dat.cytogenetics[,feature.plus]==1) ,k] <- 3

                # minus
                feature.minus <- paste0("minus", chr)
                if (feature.minus %in% colnames(dat.cytogenetics))
                        CN.df[ which(dat.cytogenetics[,feature.minus]==1) ,k] <- 1

                # add
                feature.add <- paste0("add_", chr, band)
                if (feature.add %in% colnames(dat.cytogenetics))
                        CN.df[ which(dat.cytogenetics[,feature.add]==1) ,k] <- 3

                # del
                feature.del <- paste0("del_", chr, band)
                if (feature.del %in% colnames(dat.cytogenetics))
                        CN.df[ which(dat.cytogenetics[,feature.del]==1) ,k] <- 1

        }

        return(CN.df)
       
}
