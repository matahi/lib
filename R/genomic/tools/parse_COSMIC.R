parse_COSMIC <- function(cosmic.info, DISEASE="HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", som.cutoff=1, all.cutoff=5, all.cutoff.hard=10, disease.cutoff=3)
{
        full.info <- unlist(strsplit(cosmic.info, "\\|"))

        # Verify if no COSMIC info: NA can be replaced by F but we put NA here to separate those that really dont pass the COSMIC vs those that are not in COSMIC
        if (all(is.na(full.info)))
                return(NA)

        # Verify if any genomic_EXACT + Somatic mutation
        # GENOMIC_EXACT
        genomic_exact.index <- grep("GENOMIC_EXACT", full.info)

        if (length(genomic_exact.index)>0)
        {
                disease.index <- grep(DISEASE, full.info[genomic_exact.index])
                if (length(disease.index)>0)
                        return(T)

                ALL.info.genomic_exact <- as.numeric(sub(".*;ALL\\=([0-9]+);.*",'\\1',full.info[genomic_exact.index]))  

                if (any(ALL.info.genomic_exact >= all.cutoff))
                        return(T)

                SOM.info.genomic_exact <- as.numeric(sub(".*;SOM\\=([0-9]+);.*",'\\1',full.info[genomic_exact.index]))  

                if (any(SOM.info.genomic_exact >= som.cutoff))
                        return(T)
        }

        # PROTEIN_EXACT 
        protein_exact.index <- grep("PROTEIN_EXACT", full.info)

        if (length(protein_exact.index)>0)
        {
                disease.index <- grep(DISEASE, full.info[protein_exact.index])
                if (length(disease.index)>0)
                        return(T)

                ALL.info.protein_exact <- as.numeric(sub(".*;ALL\\=([0-9]+);.*",'\\1',full.info[protein_exact.index]))  

                if (any(ALL.info.protein_exact >= all.cutoff))
                        return(T)

                SOM.info.protein_exact <- as.numeric(sub(".*;SOM\\=([0-9]+);.*",'\\1',full.info[protein_exact.index]))  

                if (any(SOM.info.protein_exact >= som.cutoff))
                        return(T)
        }

        # GENOMIC_POS
        genomic_pos.index <- grep("GENOMIC_POS", full.info)

        if (length(genomic_pos.index)>0)
        {
                disease.index <- grep(DISEASE, full.info[genomic_pos.index])
                if (length(disease.index)>0)
                        return(T)
        }


        # Otherwise grep the disease
        disease.index <- grep(DISEASE, full.info)

        if (length(disease.index)==0)
        {
                ALL.info <- as.numeric(sub(".*;ALL\\=([0-9]+);.*",'\\1',full.info))  

                if (any(ALL.info >= all.cutoff.hard))
                        return(T)

                SOM.info <- as.numeric(sub(".*;SOM\\=([0-9]+);.*",'\\1',full.info))  

                if (any(SOM.info >= all.cutoff))
                        return(T)

                return(F)
        } else
        {
                ALL.info <- as.numeric(sub(".*;ALL\\=([0-9]+);.*",'\\1',full.info[disease.index])) 
                SOM.info <- as.numeric(sub(".*;SOM\\=([0-9]+);.*",'\\1',full.info[disease.index])) 

                index <- any(((ALL.info>=all.cutoff)|(SOM.info>=som.cutoff)))
                if (index)
                        return(index)
        }

        # Recover how many of disease 
        # Verify if disease often annotated
        disease.index <- grep(DISEASE, full.info)
        if (length(disease.index)>=disease.cutoff)
        {
                return(T)
        } else
        {
                return(F)
        }



}
