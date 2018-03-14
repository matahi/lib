get_gene_infos <- function (gene.list) {

        library(biomaRt)
        ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset= "hsapiens_gene_ensembl")

        tmp <- listAttributes(mart=ensembl)

        gene.infos <- getBM(attributes=c("hgnc_symbol","chromosome_name","band"), filters = "hgnc_symbol", values = gene.list, mart= ensembl)

        gene.df <- data.frame(matrix(NA, nrow=length(gene.list), ncol=3))
        colnames(gene.df) <- c("Gene", "Chr", "Band")

        gene.df$Gene <- gene.list
        gene.df$Chr <- gene.infos$chromosome_name[match(gene.list,gene.infos$hgnc_symbol)]
        gene.df$Band <- gene.infos$band[match(gene.list,gene.infos$hgnc_symbol)]

        return(gene.df)
}

estimate_CN <- function (dat.cytogenetics, gene.list)
{

        # Get gene.infos from gene.list
        gene.infos <- get_gene_infos(gene.list)

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
