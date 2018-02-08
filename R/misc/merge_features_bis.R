merge_features_bis <- function(dat, feature.list, feature.name)
{
        feature.present <- feature.list[feature.list %in% colnames(dat)]
        if (length(feature.present)>0)
        {
                merged.feature <- as.numeric(rowSums(dat[,feature.present, drop=F],na.rm=T)>=1)

                # Place at the good position
                dat[,feature.name] <- merged.feature
                dat[,feature.present] <- NULL
        }

        return(dat)

}
