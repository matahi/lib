normalize_counts <- function(Dat, ERCC, spikeins_only=T)
{
        require(statmod)
        require(gdata)
        require(genefilter)
        require(EBImage)
        require(rhdf5)
        require(DESeq)
        require(statmod)
        require(hom.Hs.inp.db)
        require(AnnotationDbi)
        require(org.Mm.eg.db)
        require(org.Hs.eg.db)


        load('data/Tcell/data_Tcells.Rdata')

        #2. calculate normalisation for counts
        # countsMmus <- Dat.processed[ which( geneTypes=="ENSM" ), ]
        # countsERCC <- Dat.processed[ which( geneTypes=="ERCC" ), ]

        countsMmus <- Dat
        countsERCC <- ERCC 

        sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
        if (spikeins_only)
        {
                sfMmus <- sfERCC #also use ERCC size factor for endogenous genes
        } else 
        {
                sfMmus <- estimateSizeFactorsForMatrix( countsMmus )
        }


        #normalise read counts
        nCountsERCC <- t( t(countsERCC) / sfERCC )
        nCountsMmus <- t( t(countsMmus) / sfMmus )

       return(nCountsMmus) 


}

