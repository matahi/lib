library(genefilter)
library(statmod)
require(ggplot2)
library(gplots)
require(DESeq2)
# library(scLVM)
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



setwd('~/Documents/Work/Postdoc/projects/scd/')

#We also need to set the limix path; if you installed, you can delete the last 2 lines, otherwise you need to adapt the path.
# limix_path = '/Users/flo/software/limix-master/build/release.darwin/interfaces/python'
# configLimix(limix_path)

# data(data_Tcells)                      # 
load('data/Tcell/data_Tcells.Rdata')
help(data_Tcells)

#########################################

dataMouse[ 1:5, 1:4 ]

geneTypes <- factor( c( ENSM="ENSM", ERCC="ERCC" )[
  substr( rownames(dataMouse), 1, 4 ) ] )

#2. calculate normalisation for counts
countsMmus <- dataMouse[ which( geneTypes=="ENSM" ), ]
countsERCC <- dataMouse[ which( geneTypes=="ERCC" ), ]
lengthsMmus <- dataMouse[ which( geneTypes=="ENSM" ), 1 ]
lengthsERCC <- dataMouse[ which( geneTypes=="ERCC" ), 1 ]

sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
sfMmus <- sfERCC #also use ERCC size factor for endogenous genes


#normalise read counts
nCountsERCC <- t( t(countsERCC) / sfERCC )
nCountsMmus <- t( t(countsMmus) / sfERCC )

#get technical noise
techNoise = fitTechnicalNoise(nCountsMmus,nCountsERCC=nCountsERCC, fit_type = 'counts')  

techNoiseLogFit = fitTechnicalNoise(nCountsMmus, fit_type = 'log', use_ERCC = FALSE, plot=FALSE) 

techNoiseLogVarFit = fitTechnicalNoise(nCountsMmus, fit_type = 'logvar', use_ERCC = FALSE, plot=FALSE) 

#call variable genes
is_het = getVariableGenes(nCountsMmus, techNoise$fit, method = "fdr", 
                          threshold = 0.1, fit_type="counts",sfEndo=sfMmus, sfERCC=sfERCC)
table(is_het)

#we an also do this for the other fits

is_hetLog = getVariableGenes(nCountsMmus, techNoiseLogFit$fit, plot=TRUE)
table(is_hetLog)

is_hetLogVar = getVariableGenes(nCountsMmus, techNoiseLogVarFit$fit, plot=TRUE)
table(is_hetLogVar)

#get cell cycle genes from GO 
ens_ids_cc <- getEnsembl('GO:0007049')


#rename a few variables
Y = t(log10(nCountsMmus+1)) #normalised trandformed read counts
genes_het_bool = as.vector(is_het) #variable genes
geneID = rownames(nCountsMmus) #gene IDs
tech_noise = as.vector(techNoise$techNoiseLog) #technical noise


