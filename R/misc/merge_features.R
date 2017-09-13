merge_features <- function(dat, feature1, feature2)
{
        merged.feature <- sapply(1:nrow(dat), function(n)
                                 {
                                         if (is.na(dat[n,feature1]))
                                         {
                                                 return(dat[n,feature2])
                                         } else if (is.na(dat[n,feature2]))
                                         {
                                                 return(dat[n,feature1])
                                         } else if ((dat[n,feature1]==1)|(dat[n,feature2]==1))
                                         {
                                                 return(1)
                                         } else
                                         {
                                                 return (0)
                                         }
                                 })

        # Place at the good position
        Index <- min(which(colnames(dat) %in% c(feature1,feature2)))
        Index.bis <- max(which(colnames(dat) %in% c(feature1,feature2)))
        dat[,Index] <- merged.feature
        dat <- dat[ ,-Index.bis]
        colnames(dat)[Index] <- paste0(feature1,"_",feature2)

        return(dat)

}
