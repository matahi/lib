venn_list <- function(gene.list, out.dir="")
{
        # Not optimal but at least its working
        index.list <- unlist(lapply(1:length(gene.list), function(k)
                                    {
                                     combn(names(gene.list), length(gene.list)-k+1, simplify = FALSE)
                                    }),recursive = FALSE)



        comparison.list <- list()

        for (i in 1:length(index.list))
        {
                comparison.list[[i]] <- setdiff(Reduce("intersect",gene.list[index.list[[i]]]), Reduce("union",comparison.list))
                names(comparison.list)[i] <- paste(index.list[[i]],collapse="_")
        }

        
        return(comparison.list)


}
