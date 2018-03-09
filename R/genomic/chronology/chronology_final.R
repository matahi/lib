chronology_final <- function (dat.genetics, dat.binary) {

        pairs.genes <- t(combn(colnames(dat.binary),2))

        source("./src/lib/R/genomic/chronology/analyze_chronology_bis.R")
        # tmp <- lapply(1:nrow(pairs.genes), function(k) {
        tmp <- lapply(1:nrow(pairs.genes), function(k) { analyze_chronology(dat.binary, dat.genetics, pairs.genes[k,1], pairs.genes[k,2])})


        chronology.df <- Reduce("rbind",lapply(1:length(tmp), function(k){wins1 <- length(grep(pairs.genes[k,1], tmp[[k]]));wins2 <- length(grep(pairs.genes[k,2], tmp[[k]])); draws <- sum(tmp[[k]]=="inconclusive");
                                               return(data.frame(Gene1=pairs.genes[k,1],
                                                                 Gene2=pairs.genes[k,2],
                                                                 Total = length(tmp[[k]]),
                                                                 Win=wins1,
                                                                 Loss=wins2,
                                                                 Inconclusive=draws,
                                                                 stringsAsFactors=F))})) %>% tbl_df()
        
        return(chronology.df)

}
