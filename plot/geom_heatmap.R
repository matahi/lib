geom_heatmap <- function(dat)
{
        #dat <- data.frame(dat)

        order.samples <- order(rowSums(dat), decreasing=T)
        order.genes <- order(colSums(dat))


        library(reshape2)
        dat.melt <- melt(dat)

        dat.melt$Var1 <- factor(dat.melt$Var1, levels=rownames(dat)[order.samples])
        dat.melt$Var2 <- factor(dat.melt$Var2, levels=colnames(dat)[order.genes])

        dat.melt$value[dat.melt$value==0] <- NA

        source("./src/lib/plot/pie_chart.R")
        p <- ggplot(dat.melt) + geom_tile(aes(x=Var1,y=Var2,fill=value), colour="black") + blank_theme +
              xlab("Samples") + ylab("Genes") +  theme(axis.text.x=element_text(angle=90,hjust=1), legend.position="none")

        return(p)

}
