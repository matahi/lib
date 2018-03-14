chronology_final <- function (dat.genetics, dat.binary) {

        pairs.genes <- t(combn(setdiff(colnames(dat.binary),"ID"),2))

        source("./src/lib/R/genomic/chronology/analyze_chronology_bis.R")
        tmp <- lapply(1:nrow(pairs.genes), function(k) { analyze_chronology(dat.binary, dat.genetics, pairs.genes[k,1], pairs.genes[k,2])})

        chronology.df <- Reduce("rbind",lapply(1:length(tmp), function(k){wins1 <- length(grep(pairs.genes[k,1], tmp[[k]]));wins2 <- length(grep(pairs.genes[k,2], tmp[[k]])); draws <- sum(tmp[[k]]=="inconclusive");
                                               return(data.frame(Gene1=pairs.genes[k,1],
                                                                 Gene2=pairs.genes[k,2],
                                                                 Total = length(tmp[[k]]),
                                                                 Win=wins1,
                                                                 Loss=wins2,
                                                                 Inconclusive=draws,
                                                                 stringsAsFactors=F))})) %>% tbl_df()

        # Lets make the chronology.df symmetric (for future plotting purposes)
        chronology.reverse <- chronology.df %>%
                mutate(Gene1_bis = Gene2,
                       Gene2 = Gene1,
                       Gene1 = Gene1_bis,
                       Loss_bis = Win,
                       Win = Loss,
                       Loss = Loss_bis
                       ) %>%
                select(colnames(chronology.df))

        chronology.diag <- data.frame(Gene1= unique(union(chronology.df$Gene1, chronology.df$Gene2)),
                                      Gene2= unique(union(chronology.df$Gene1, chronology.df$Gene2)),
                                      Total = NA,
                                      Win = NA,
                                      Loss = NA,
                                      Inconclusive = NA, stringsAsFactors = F) %>% tbl_df()

        chronology.df.full <- bind_rows(chronology.df, 
                                        chronology.reverse,
                                        chronology.diag)

        
        return(chronology.df.full)
}
