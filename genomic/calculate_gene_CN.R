calculate_gene_CN <- function(dat.cytogenetics, gene, chr, branch)
{

        # CN loss
        cols.neg <- c(paste0("-",chr), paste0("del(",chr,")([?]*",branch))

        Index.neg <- which(colnames(dat.cytogenetics) %in% cols.neg)

        if (length(Index.neg)==0)
        {
               CN.neg <- rep(0,nrow(dat.cytogenetics))
        } else if (length(Index.neg)==1)
        {
               CN.neg <- -dat.cytogenetics[,Index.neg]
        } else
        {
                cyto.infos <- dat.cytogenetics[, Index.neg, drop=F ]
                # CN.neg = or
                CN.neg <- -apply(cyto.infos,1,function(x){any(x>0,na.rm=T)})
        }

        # CN add
        cols.pos <- c(paste0("+",chr), paste0("add(",chr,")([?]*",branch))

        Index.pos <- which(colnames(dat.cytogenetics) %in% cols.pos)

        if (length(Index.pos)==0)
        {
                CN.pos <- rep(0,nrow(dat.cytogenetics))
        } else if (length(Index.pos)==1)
        {
                CN.pos <- dat.cytogenetics[,Index.pos]
        } else
        {
                cyto.infos <- dat.cytogenetics[, Index.pos, drop=F ]
                # CN.pos = or
                CN.pos <- apply(cyto.infos,1,function(x){any(x>0,na.rm=T)})
        }

        #
        CN.effect <- CN.neg + CN.pos

        return(CN.effect)

}
