filter_residue <- function (residueSeq) {
        if (any(residueSeq %in% c("OPAL", "OCHRE", "AMBER"))) {
                stopRes <- c("OPAL", "OCHRE", "AMBER")
                residueSeq <- residueSeq[-which(residueSeq %in% stopRes)]
        }

        return(residueSeq)
        
}
