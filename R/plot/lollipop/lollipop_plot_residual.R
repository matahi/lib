## # GenVisR:::lolliplot_buildMain
## gene_plot <- geom_rect(data = gene_data[1, ], mapping = aes_string(xmin = "pos_from", 
##                                                               xmax = "pos_to", ymin = "height_min", ymax = "height_max"), 
##                        fill = "#999999", colour = "#000000")
## if (nrow(gene_data) == 1) {
##         domain_plot <- geom_blank()
## } else {
##         domain_plot <- geom_rect(data = gene_data[-1, ], 
##                                  mapping = aes_string(xmin = "pos_from", xmax = "pos_to", 
##                                                       ymin = "height_min", ymax = "height_max", fill = "Domain"), 
##                                  alpha = 1, colour = "black")
## }
## #
## 
## title <- ggtitle(gene)
## x_label <- xlab("Amino Acid Position")
## theme_protein <- theme(legend.position = "bottom", legend.direction = "vertical", 
##                legend.box = "horizontal", axis.text.y = element_blank(), 
##                axis.ticks.y = element_blank(), axis.title.y = element_blank())
## guide <- guides(colour = guide_legend(ncol = 2), fill = guide_legend(ncol = 2))
## #
## 
## source("./src/lib/plot/lollipop/theme_protein.R")
## p.protein <- ggplot() + gene_plot + domain_plot + xlab("Amino Acid Position") + theme_protein() + theme(legend.position="none")
## 
## ## 2. Recover gene mutations
## dat.genetics.gene$PROTEIN_POS <- as.numeric(dat.genetics.gene$PROTEIN_POS)

###  ##################
###  # TODO: Method 1
###  ##################
###  # Trackviewer
###  library(trackViewer)
###  
###  # i) Mutation infos
###  #####################
###  #SNP <- c(20,200,500,500)
###  PROTEIN_CHANGE <- dat.genetics.gene$PROTEIN_CHANGE
###  
###  mut.infos <- data.frame(table(PROTEIN_CHANGE=dat.genetics.gene$PROTEIN_CHANGE, disease=dat.genetics.gene$disease))
###  mut.infos$PROTEIN_POS <- sapply(1:nrow(mut.infos), function(n)
###                                              {
###                                                      PROTEIN_NUM <- strsplit(as.character(mut.infos$PROTEIN_CHANGE[n]),"[^0-9]+")[[1]][2]
###                                                      return(as.numeric(PROTEIN_NUM))
###                                              })
###  
###  sample.gr <- GRanges(gene, IRanges(mut.infos$PROTEIN_POS, width=1, names=mut.infos$PROTEIN_CHANGE), disease=mut.infos$disease)
###  
###  color.type <- "CONSEQUENCE"
###  if (is.null(color.type))
###  {
###          color.panel <- NULL
###          sample.gr$color <- 1
###  } else if (color.type == "CONSEQUENCE")
###  {
###          color.panel <- colors.consequence_type
###          consequences <- dat.genetics.gene$CONSEQUENCE[match(names(sample.gr), dat.genetics.gene$PROTEIN_CHANGE)]
###  
###          #
###          sample.gr$type[consequences %in% c("non_synonymous_codon")] <- "missense"
###          sample.gr$type[!consequences %in% c("non_synonymous_codon")] <- "nonsense"
###          sample.gr$color <- color.panel[sample.gr$type]
###  } else if (color.type =="type")
###  {
###          # TODO
###          color.panel <- colors.variant_type
###  }
###  #sample.gr$color <- sample.int(6, nrow(mut.infos), replace=TRUE) # Colors depending of type (i.e missense vs nonsense)
###  
###  sample.gr$score <- mut.infos$Freq
###  sample.gr$dashline.col <- sample.gr$color
###  
###  # Possibility 1: double plot
###  sample.gr$SNPsideID[sample.gr$disease=="AML"] <- "top"
###  sample.gr$SNPsideID[sample.gr$disease=="MDS"] <- "bottom"
###  
###  # Possibility 2: multiple plot
###  # split sample.gr into sample.gr1 and sample.gr2
###  
###  # Create legend
###  if (is.null(color.panel))
###  {
###          legends <- NULL
###  } else
###  {
###          legends <- list(labels=names(color.panel), fill=color.panel)
###  }
###  
###  
###  # ii) Protein infos
###  ##################
###  features <- GRanges(gene, IRanges(c(1,100,190),
###                                    width=c(30,14,22),
###                                    names=paste0(gene,1:3)))
###  features$fill <- c("#FF8833", "#51C6E6", "#DFA32D")
###  features$height <- c(0.02, 0.05, 0.08) # Change the size of domains
###  
###  
###  # iii) Plot infos
###  #################
###  #xaxis.infos <- c(1, seq(200, proteinLength, by=200), proteinLength) # xticks
###  xaxis.infos <- c(1, seq(200, proteinLength, by=200)) # xticks
###  lolliplot(sample.gr, 
###            features, 
###            ranges=GRanges(gene, IRanges(0,proteinLength+1)),
###            xaxis=xaxis.infos,
###            legend=legends, # Add legend info
###            jitter="label")
###  
###  # title
###  grid.text(gene, x=.5, y=.98, just="top", 
###            gp=gpar(cex=1.5, fontface="bold"))

####################################################################
##################
# TODO: Method 2
##################
# Substitutions 
dat.subs <- dat.genetics.gene[dat.genetics.gene$type=="Sub",]
tmp.subs <- data.frame(table(dat.subs$PROTEIN_POS, dat.subs$disease))
tmp.subs$Var1 <- as.numeric(as.character(tmp.subs$Var1))
tmp.subs$Freq[tmp.subs$Var2=="MDS"] <- -tmp.subs$Freq[tmp.subs$Var2=="MDS"]

### if tmp.subs
if (nrow(tmp.subs)==0)
{
        tmp.subs <- data.frame(Freq=0, Var1=0)

        p.subs <- ggplot(tmp.subs) + 
                 geom_histogram(aes(x=Var1,y=Freq), colour="black", stat="identity") + 
                 geom_hline(aes(yintercept=0), colour="black") + 
                 theme_bw()  + ylim(-1,1)

} else
{
        p.subs <- ggplot(tmp.subs) + 
                 geom_histogram(aes(x=Var1,y=Freq), colour="black", stat="identity") + 
                 geom_hline(aes(yintercept=0), colour="black") + 
                 theme_bw()

}

# Insertions
dat.Ins <- dat.genetics.gene[dat.genetics.gene$type == "Ins",]

tmp.Ins <- data.frame(table(dat.Ins$PROTEIN_POS, dat.Ins$disease))
tmp.Ins$Freq[tmp.Ins$Var2=="MDS"] <- -tmp.Ins$Freq[tmp.Ins$Var2=="MDS"]
tmp.Ins$Var1 <- as.numeric(as.character(tmp.Ins$Var1))

if (nrow(tmp.Ins)==0)
{
        tmp.Ins <- data.frame(Freq=0, Var1=0)

        p.Ins <- ggplot(tmp.Ins) + 
                 geom_histogram(aes(x=Var1,y=Freq), colour="black", stat="identity") + 
                 geom_hline(aes(yintercept=0), colour="black") + 
                 theme_bw()  + ylim(-1,1)

} else
{
        p.Ins <- ggplot(tmp.Ins) + 
                 geom_histogram(aes(x=Var1,y=Freq), colour="black", stat="identity") + 
                 geom_hline(aes(yintercept=0), colour="black") + 
                 theme_bw()

}

# Deletions
dat.Del <- dat.genetics.gene[dat.genetics.gene$type == "Del",]

tmp.Del <- data.frame(table(dat.Del$PROTEIN_POS, dat.Del$disease))
tmp.Del$Freq[tmp.Del$Var2=="MDS"] <- -tmp.Del$Freq[tmp.Del$Var2=="MDS"]
tmp.Del$Var1 <- as.numeric(as.character(tmp.Del$Var1))

if (nrow(tmp.Del)==0)
{
        tmp.Del <- data.frame(Freq=0, Var1=0)

        p.Del <- ggplot(tmp.Del) + 
                 geom_histogram(aes(x=Var1,y=Freq), colour="black", stat="identity") + 
                 geom_hline(aes(yintercept=0), colour="black") + 
                 theme_bw()  + ylim(-1,1)

} else
{
        p.Del <- ggplot(tmp.Del) + 
                 geom_histogram(aes(x=Var1,y=Freq), colour="black", stat="identity") + 
                 geom_hline(aes(yintercept=0), colour="black") + 
                 theme_bw()

}

## 2. Recover gene infos

# tracks(substitutions=p.subs, p.protein, p.indels, heights=c(3,1,3), main=gene)
library(ggbio)
pdf(paste0("./results/analysis/gene/",gene,"/lollipop.pdf"))
print(tracks(substitutions=p.subs, p.protein, Insertions=p.Ins, Deletions=p.Del, heights=c(3,1,3,3), main=gene, xlim=c(-1,proteinLength + 5)))
dev.off()


