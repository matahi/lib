constructGene <- function (gene, domain_data, gene.length) {
        # domain_data <- protein_domain
        # gene.length <- proteinLength

        gene <- data.frame(Domain=gene, pos_from=1, pos_to=gene.length, nest=1) 
        if (nrow(na.omit(domain_data)) == 0) {
                gene$height_min <- 0.1/(as.numeric(gene$nest))
                gene$height_max <- -0.1/(as.numeric(gene$nest))
                gene$pos_from <- as.numeric(gene$pos_from)
                gene$pos_to <- as.numeric(gene$pos_to)
                return(gene)
        }

        ## Processing the domain data
        colnames(domain_data) <- c("description", "start", "end")
        domain_data$description <- as.character(domain_data$description)
        if (max(domain_data$end) > gene.length) {
                memo <- paste0("The end position of a domain: ", max(domain_data$end), 
                               " is exceeding the length of the protein:", gene.length)
                warning(memo)
        } else if (min(domain_data$start) < 1) {
                memo <- paste0("The start position of a domain:", min(domain_data$start), 
                               "is less than the start of the protein", 1)
                warning(memo)
        }
        if (any(domain_data$start >= domain_data$end)) {
                memo <- paste0("Found a start position greater than an end position", 
                               " in the protein features track. Check input to Z or", 
                               "results of the biomaRt query using dataOut==TRUE.")
                warning(memo)
        }

        domain_data$start <- as.numeric(domain_data$start)
        domain_data$end <- as.numeric(domain_data$end)
        domain_data <- domain_data[order(domain_data$start), ]
        nest <- vector("numeric")
        end <- vector("numeric")
        for (i in 1:nrow(domain_data)) {
                idx <- domain_data$start[i] < end
                end <- end[idx]
                nest <- c(nest, length(end))
                end <- c(end, domain_data$end[i])
        }
        domain_data$nest <- nest + 1
        colnames(domain_data) <- c("Domain", "pos_from", "pos_to", 
                                   "nest")
        gene <- rbind(gene, domain_data)
        gene$height_min <- 0.1/(as.numeric(gene$nest))
        gene$height_max <- -0.1/(as.numeric(gene$nest))
        gene$pos_from <- as.numeric(gene$pos_from)
        gene$pos_to <- as.numeric(gene$pos_to)
        return(gene)
}
