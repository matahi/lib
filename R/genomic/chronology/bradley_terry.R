bradley_terry <- function (chronology.df, group.genes = NULL) {

        # chronology.df

        library(BradleyTerry2)

        genes <- unique(chronology.df$Gene1)

        precedence.df <- chronology.df %>%
                mutate(Gene1 = factor(Gene1, levels=genes),
                       Gene2 = factor(Gene2, levels=genes)) %>%
                filter( as.numeric(Gene1) > as.numeric(Gene2))

        btModel <- BTm(cbind(Win, Loss), Gene1, Gene2, data = precedence.df)

        bt.score <- BTabilities(btModel)

        bt.df <- data.frame(Gene=rownames(bt.score), bt.score) %>% tbl_df() 

        bt.df <- bt.df %>% mutate(Gene = factor(Gene, levels=as.character(bt.df$Gene[order(bt.df$ability)])))

        bt.df.bis <- bt.df %>% filter(s.e. < 2) 
        bt.df.bis <- bt.df.bis %>% mutate(Gene = factor(Gene, levels=levels(bt.df$Gene)[levels(bt.df$Gene) %in% bt.df.bis$Gene]))
        cols.axis <- group.genes$colours[ match(levels(bt.df.bis$Gene), group.genes$alterations )] 

        p.bt <- bt.df.bis %>% ggplot() + 
                geom_point(aes(x=Gene, y=-ability)) +
                geom_segment(aes(x=Gene, xend=Gene, y = -ability - s.e., yend = -ability + s.e.)) +
                theme(axis.text.y = element_text(color=cols.axis)) +
                coord_flip()

        # TODO: bt with ties 
        bt.out <- list(model = bt.df, p=p.bt)
        return(bt.out)
}




# makeDesign <- function(I) {
#         w <- which(lower.tri(I), arr.ind=TRUE)
#         x <- matrix(0, nrow(w), nrow(I))
#         for(i in 1:nrow(w)){
#                 x[i,w[i,1]] <- 1
#                 x[i,w[i,2]] <- -1
#         }
#         return(x)
# }
# 
# btModel <- function(I){
#         y <- cbind(I[lower.tri(I)], t(I)[lower.tri(I)])
#         x <- makeDesign(I = I)
#         glm.fit(x=x[,-1],y=y, family=binomial())
# }
# 
# nCasesGene <- table(factor(unlist(sapply(plist, function(x) unique(unlist(x)))), levels=colnames(precedence)))
# w <- which(nCasesGene > 5)
# fit <- btModel(precedence[w,w]+.01)
# 
# ## Warning: non-integer counts in a binomial glm!
# 
# c <- c(0,coef(fit))
# names(c) <- colnames(precedence)[w]
# o <- rank(c)
# v <- pmin(2,sqrt(c(0,diag(chol2inv(fit$qr$qr))))) # this is the estimate of the Beta_i CI
# 
# l <- names(c)
# m <- paste("n=",nCasesGene[w], sep="")
# plot(-c, o, xlab="Relative time", yaxt="n", pch=19, col="grey", ylab="", xlim=range(-c+3*c(-v,v)))
# segments(-c-v, o,-c+v,o, col="grey")
# text(-c-v ,o,l, font=3, pos=2)
# text(-c+v ,o,m, font=1, pos=4)

