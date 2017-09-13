technoise <- function(Dat.processed, Path, ERCC.quant=0.8, normalize=T, spikeins_only=T, soft.filter=NULL)
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

        geneTypes <- factor( c( ENSM="ENSM", ERCC="ERCC" )[
                            substr( rownames(Dat.processed), 1, 4 ) ] )

        #2. calculate normalisation for counts
        countsMmus <- Dat.processed[ which( geneTypes=="ENSM" ), ]
        countsERCC <- Dat.processed[ which( geneTypes=="ERCC" ), ]
        lengthsMmus <- Dat.processed[ which( geneTypes=="ENSM" ), 1 ]
        lengthsERCC <- Dat.processed[ which( geneTypes=="ERCC" ), 1 ]

        if (normalize)
        {
                sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
                if (spikeins_only)
                {
                        sfMmus <- sfERCC #also use ERCC size factor for endogenous genes
                } else {
                        sfMmus <- estimateSizeFactorsForMatrix( countsMmus )
                }


                #normalise read counts
                nCountsERCC <- t( t(countsERCC) / sfERCC )
                nCountsMmus <- t( t(countsMmus) / sfMmus )
        } else {
                nCountsERCC <- t( t(countsERCC) )
                nCountsMmus <- t( t(countsMmus) )
        }

        #normalized counts (brennecke)
        meansMmus <- rowMeans( nCountsMmus )
        varsMmus <- rowVars( nCountsMmus )
        cv2Mmus <- varsMmus / meansMmus^2

        meansERCC <- rowMeans( nCountsERCC )
        varsERCC <- rowVars( nCountsERCC )
        cv2ERCC <- varsERCC / meansERCC^2

        # pdf(paste0("results/Brennecke/Var_Mean/",Path, ".pdf"))
        # plot(log(meansMmus), log(varsMmus))
        # points(log(meansERCC), log(varsERCC),col='red')
        # dev.off()

        #Do fitting of technical noise

        #normalised counts (with size factor)
        minMeanForFitA <- unname( quantile( meansERCC[ which( cv2ERCC > .3 ) ], ERCC.quant) )
        # minMeanForFitA <- unname( quantile( meansERCC[ which( cv2ERCC > .3 ) ], .4 ) )
        useForFitA <- meansERCC >= minMeanForFitA
        fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCC[useForFitA] ),
                           cv2ERCC[useForFitA] )

        #plot fit

        pdf(paste0("results/Brennecke/CV2/",Path, "_", ERCC.quant, "_ERCC.pdf"))
        # pdf(paste0("results/Brennecke/CV2/",Path,"_ERCC_bis.pdf"))
        plot( meansERCC, cv2ERCC, log="xy", col=1+useForFitA)
        xg <- 10^seq( -3, 5, length.out=100 )
        lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg )
        segments( meansERCC[useForFitA], cv2ERCC[useForFitA],
                 meansERCC[useForFitA], fitA$fitted.values, col="gray" )
        dev.off()

        # pdf(paste0("results/Brennecke/CV2/",Path, "_", ERCC.quant, "_ERCC_dem_1.pdf"))
        # # pdf(paste0("results/Brennecke/CV2/",Path,"_ERCC_bis.pdf"))
        # plot( meansERCC, cv2ERCC, log="xy")
        # dev.off()

        # pdf(paste0("results/Brennecke/CV2/",Path, "_", ERCC.quant, "_ERCC_dem_2.pdf"))
        # # pdf(paste0("results/Brennecke/CV2/",Path,"_ERCC_bis.pdf"))
        # plot( meansERCC, cv2ERCC, log="xy", col=1+useForFitA)
        # dev.off()

        # pdf(paste0("results/Brennecke/CV2/",Path, "_", ERCC.quant, "_ERCC_dem_3.pdf"))
        # # pdf(paste0("results/Brennecke/CV2/",Path,"_ERCC_bis.pdf"))
        # plot( meansERCC, cv2ERCC, log="xy", col=1+useForFitA)
        # xg <- 10^seq( -3, 5, length.out=100 )
        # lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg )
        # segments( meansERCC[useForFitA], cv2ERCC[useForFitA],
        #          meansERCC[useForFitA], fitA$fitted.values, col="gray" )
        # dev.off()






        #perfrom statistical test
        minBiolDisp <- .5^2
        xi <- mean( 1 / sfERCC )
        m <- ncol(countsMmus)
        psia1thetaA <- mean( 1 / sfERCC ) +
        ( coefficients(fitA)["a1tilde"] - xi ) * mean( sfERCC / sfMmus )
        cv2thA <- coefficients(fitA)["a0"] + minBiolDisp + coefficients(fitA)["a0"] * minBiolDisp
        testDenomA <- ( meansMmus * psia1thetaA + meansMmus^2 * cv2thA ) / ( 1 + cv2thA/m )
        pA <- 1 - pchisq( varsMmus * (m-1) / testDenomA, m-1 )
        padjA <- p.adjust( pA, "BH" )
        table( padjA < .1 )

        save(padjA, file=paste0('results/Brennecke/variable_genes/',Path, '_', ERCC.quant, '_padjA.RData'))
        write(table( padjA < .1 ),file=paste0('results/Brennecke/variable_genes/', Path, "_", ERCC.quant, '.txt'))

        #get cell cycle genes from GO 
        xxGO <- as.list(org.Mm.egGO2EG)
        cell_cycleEG <-unlist(xxGO['GO:0007049'])
        #get ENSEMBLE ids
        x <- org.Mm.egENSEMBL
        mapped_genes <- mappedkeys(x)
        xxE <- as.list(x[mapped_genes])
        ens_ids_cc<-unlist(xxE[cell_cycleEG])
        #cc_gene_indices <- na.omit(match(ens_ids_cc, rownames(dataMouse)))
        cc_gene_indices <- na.omit(match(ens_ids_cc, rownames(Dat.processed)))


        #plot mean/cv2 relationship and 
        pdf(paste0("results/Brennecke/CV2/",Path,"_", ERCC.quant, ".pdf"))
        # pdf(paste0("results/Brennecke/CV2/",Path,"_bis.pdf"))
        plot( meansMmus, cv2Mmus, log="xy", col=1+(padjA<0.1),ylim=c(0.1,95), xlab='Mean Counts', ylab='CV2 Counts')
        xg <- 10^seq( -3, 5, length.out=100 )
        lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg,lwd=2,col='blue' )
        points(meansERCC, cv2ERCC,col='blue',pch=15,cex=1.1)
        #points(meansMmus[cc_gene_indices], cv2Mmus[cc_gene_indices],col=rgb(0,255,0,100,maxColorValue=255),pch=2,cex=0.75)
        points(meansMmus[cc_gene_indices], cv2Mmus[cc_gene_indices],col="green",pch=2,cex=0.75)
        #points(meansMmus[ccCBall_gene_indices], cv2Mmus[ccCBall_gene_indices],col=rgb(0,255,0,20,maxColorValue=255),pch=2,cex=0.8)
        legend('bottomleft',c('Non-Variable Genes (padj >= 0.1)','Variable Genes (padj<0.1)','ERCC','Cell Cycle genes'),pch=c(1,1,15),col=c('black','red','blue','green'),cex=0.7)
        dev.off()

        # #  Supplementary plot with Wilson genes
        Wilson.ensembl <- get(load("data/Info/Bertie_ensembl_bis.RData")) # Wilson signature
        #wilson_gene_indices <- na.omit(match(Wilson.ensembl, rownames(dataMouse)))
        wilson_gene_indices <- match(Wilson.ensembl, rownames(Dat.processed))

        #plot mean/cv2 relationship and 
        pdf(paste0("results/Brennecke/CV2/",Path,"_", ERCC.quant, "_wilson.pdf"))
        plot( meansMmus, cv2Mmus, log="xy", col=1+(padjA<0.1),ylim=c(0.1,95), xlab='Mean Counts', ylab='CV2 Counts')
        #plot(meansMmus[wilson_gene_indices], cv2Mmus[wilson_gene_indices],ylim=c(0.1,95), xlab='Mean Counts', log="xy",ylab='CV2 Counts',col="green",pch=2,cex=0.75)
        points(meansMmus[wilson_gene_indices], cv2Mmus[wilson_gene_indices],col="green",pch=2,cex=0.75)
        xg <- 10^seq( -3, 5, length.out=100 )
        lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg,lwd=2,col='blue' )
        points(meansERCC, cv2ERCC,col='blue',pch=15,cex=1.1)
        points(meansMmus[wilson_gene_indices], cv2Mmus[wilson_gene_indices],col="green",pch=2,cex=0.75)
        legend('bottomleft',c('Non-Variable Genes (padj >= 0.1)','Variable Genes (padj<0.1)','ERCC','Wilson 48 gene signature'),pch=c(1,1,15),col=c('black','red','blue','green'),cex=0.7)
        dev.off()

        ##

        #4. Transform to log-space and propagate error
        eps=1
        LogNcountsMmus=log10(nCountsMmus+eps)
        dLogNcountsMmus=1/(meansMmus+eps)
        var_techMmus=(coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/meansMmus)*meansMmus^2
        LogVar_techMmus=(dLogNcountsMmus*sqrt(var_techMmus))^2 #error propagation 


        #gene names in the T-cell data.
        gene_names=rownames(nCountsMmus)
        gene_names_het=gene_names[which(padjA<0.1)]
        
        if (is.null(soft.filter))
        {
        } else
        {
                gene_names_het=c(gene_names_het, setdiff(soft.filter,gene_names_het))
        }


        #all Cycle base genes homologs (top 600 genes)
        hu2musAll=inpIDMapper(dataCB[1:600,3],'HOMSA','MUSMU',srcIDType='ENSEMBL',destIDType='ENSEMBL')
        ccCBall_gene_indices=match(unlist(hu2musAll),rownames(nCountsMmus))
        #lenCB=unlist(((lapply(hu2musAll,function(x){length(x)}))))

        #get cell cycle genes from GO 
        xxGO <- as.list(org.Mm.egGO2EG)
        cell_cycleEG <-unlist(xxGO['GO:0007049'])
        #get ENSEMBLE ids
        x <- org.Mm.egENSEMBL
        mapped_genes <- mappedkeys(x)
        xxE <- as.list(x[mapped_genes])
        ens_ids_cc<-unlist(xxE[cell_cycleEG])
        cc_gene_indices <- na.omit(match(ens_ids_cc, rownames(Dat.processed)))

        #ensemble IDs to gene symbols
        x <- org.Mm.egSYMBOL
        # Get the gene symbol that are mapped to an entrez gene identifiers
        mapped_genes <- mappedkeys(x)
        # Convert to a list
        xx <- as.list(x[mapped_genes])
        xxenseg <- as.list(org.Mm.egENSEMBL2EG)
        gene_syms=unlist(xx[unlist(xxenseg[gene_names])])
        gene_names_list<-(lapply(xxenseg[gene_names],function(x){if(is.null(x)){x=NA}else{x=x[1]}}))
        sym_names=unlist(lapply(xx[unlist(gene_names_list)],function(x){if(is.null(x)){x=NA}else{x=x[1]}}))
        sym_names[is.na(sym_names)]=gene_names[is.na(sym_names)]

        if (is.null(soft.filter))
        {
                sym_names_het=sym_names[which(padjA<0.1)] #gene symbols of variable genes
        } else
        {
                index <- unique(c(which(padjA<0.1),match(soft.filter,names(padjA))))
                sym_names_het=sym_names[index] #gene symbols of variable genes
        }

        #rename a few variables...
        cellcyclegenes <- ens_ids_cc
        cellcyclegenes_filter <- cc_gene_indices
        cell_names <- colnames(nCountsMmus)
        Y <- nCountsMmus
        genes_heterogen <- (padjA<0.1)*1
        if (is.null(soft.filter))
        {
        } else
        {
                genes_heterogen[soft.filter] <- 1
        }
        countsERCC_mat=as.matrix(countsERCC * 1)
        countsMmus_mat = as.matrix(countsMmus * 1)

        # h5save(ccCBall_gene_indices,gene_names,sym_names,sym_names_het,cellcyclegenes_filter,cellcyclegenes,cell_names,nCountsMmus,genes_heterogen,LogVar_techMmus,LogNcountsMmus,countsMmus_mat,sfERCC,countsERCC_mat,file=paste0('data/RNASeq/',Path,'.h5f'))
        h5save(ccCBall_gene_indices,gene_names,sym_names,sym_names_het,cellcyclegenes_filter,cellcyclegenes,cell_names,nCountsMmus,genes_heterogen,LogVar_techMmus,LogNcountsMmus,countsMmus_mat,sfERCC,countsERCC_mat,file=paste0('data/Brennecke/',Path,'_', ERCC.quant, '.h5f'))

        ##### Compare Variable Genes?
        variable_genes <- names(which(genes_heterogen==1))
        # save(variable_genes, file=paste0('results/Brennecke/variable_genes/',Path, '.RData'))
        save(variable_genes, file=paste0('results/Brennecke/variable_genes/',Path, '_', ERCC.quant, '.RData'))
        return(genes_heterogen)
}

