parse.exac <- function(exac.annot)
{
        # exac.annot <- exac.subs

        # prepare data.frames
        col.infos <- c("FIN","NFE","AFR","AMR","EAS","SAS","OTHER")

        prop.df <- data.frame(matrix(NA, nrow=nrow(exac.annot), ncol=length(col.infos)))
        n_pop.df <- data.frame(matrix(NA, nrow=nrow(exac.annot), ncol=2*length(col.infos)))
        colnames(prop.df) <- col.infos
        colnames(n_pop.df) <- unlist(lapply(col.infos, function(x){paste(x,c("Mut","Pop"))}))
        # colnames(n_pop.df) <- col.infos

        # Parse exac
        # parse.infos <- strsplit(exac.annot[,"SNP_MAF"],";")
        parse.infos <- strsplit(exac.annot$SNP_MAF,";")

        for (n in 1:length(parse.infos))
        {
                if (all(is.na(parse.infos[[n]]))) 
                {
                        prop.df[n,] <- NA
                        n_pop.df[n,] <- NA
                } else
                {
                        # proportions
                        prop.infos <- unlist(strsplit(parse.infos[[n]][1], "\\|"))
                        prop.df[n,] <- sapply(prop.infos, function(x){as.numeric(strsplit(x, ":")[[1]][2])})

                        # # n_pop
                        n_pop.infos <- unlist(strsplit(parse.infos[[n]][2], "\\|"))
                        n_pop.df[n,] <- as.numeric(unlist(strsplit(n_pop.infos, "/")))
                }
        }

        dat.df <- cbind(prop.df,n_pop.df)
        rownames(dat.df) <- rownames(exac.annot)

        return(dat.df)
}
