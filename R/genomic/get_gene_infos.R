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
