predict_clinical_mutation <- function (dat.clinical, 
                                       dat.binary, 
                                       feature) {

        require(glmnet)
        require(lars)

        ##
        # Initialize
        # dat.clinical <- dat.clinical.analysis
        # dat.binary <- dat.binary.final
        # feature <- "HB"

        ###
        x <- as.matrix(dat.binary)
        y <- dat.clinical %>% pull(feature)

        ###
        if (any(is.na(y)))
        {
                x <- x[!is.na(y),]
                y <- y[!is.na(y)]
        }

        ###
        cv.fit <- cv.glmnet(x,y)

        coef.info <- as.matrix(coef(cv.fit))
        coef.df <- data.frame(features=rownames(coef.info),
                              coef = coef.info[,1],
                              stringsAsFactors = F) %>% tbl_df()

        # Formatting coefficient labels
        coef.breaks <- c(seq(-0.1, 0, by = 0.02), seq( .Machine$double.eps, 0.1, by = 0.02))
        coef.labels <- formatC(seq(-0.1,0.1, by = 0.02))

        coef.norm <- sum(abs(coef.df$coef))
        coef.df <- coef.df %>% 
                mutate(norm.coef = coef / coef.norm,
                       bounded.coef = case_when(is.na(norm.coef) ~ as.numeric(NA),
                                                norm.coef <= -0.1 ~ -0.1,
                                                norm.coef > 0.1 ~ 0.1,
                                                T ~ norm.coef
                                                ))

        coef.df <- coef.df %>%
                mutate(value = cut( bounded.coef, label = coef.labels, breaks = coef.breaks, right=F))

        # 
        colors.coef <- c(brewer.pal(length(coef.labels)-1, "RdBu"), "white")
        names(colors.coef) <- c(setdiff(coef.labels,"0"),"0")

        # add dummy values
        # coef.df <- bind_rows( coef.df, data.frame(value=levels(coef.df$value)))

        # recover
        coef.df <- coef.df %>%
                select(features, value)
        colnames(coef.df) <- c("features",feature)

        return(coef.df)

        # Reordering <- Reordering.features

        # features.limits <- setdiff(na.omit(coef.df$features), "(Intercept)")
        # cols.axis <- Reordering$colours[ match(features.limits, Reordering$alterations )] 

        # p.predict <- coef.df %>% ggplot() + 
        #         geom_tile(aes(x=features, y="HB", fill=value), colour="grey50") +
        #                 scale_fill_manual(name="coefficient", values=colors.coef) +
        #                 scale_x_discrete(limits = features.limits) +
        #                 theme(axis.text.x = element_text( angle = 90, hjust=1, vjust= 0.5, colour= cols.axis)) +
        #                 xlab("") + ylab("")

        # return(p.predict) 

}


