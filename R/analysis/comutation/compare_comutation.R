compare_comutation <- function (dat1, dat2) {

                # Odds ratio
                odds <- sapply(1:ncol(dat2), 
                               function(i) sapply(1:ncol(dat1), 
                                                  function(j) {f<- try(fisher.test(table(dat2[,i], 
                                                                                         dat1[,j])), 
                                                                       silent=TRUE); 
                                                  if(class(f)=="try-error") 
                                                          f=NA 
                                                  else f$estimate} ))


                # Pvals
                logPInt <- sapply(1:ncol(dat2), 
                                  function(i) sapply(1:ncol(dat1), 
                                                     function(j) {f<- try(fisher.test(dat2[,i], 
                                                                                      dat1[,j]), 
                                                                          silent=TRUE); 
                                                     if(class(f)=="try-error") 
                                                             0 
                                                     else ifelse(f$estimate>1, 
                                                                 -log10(f$p.val),
                                                                 log10(f$p.val))} ))

                # 
                rownames(logPInt) <- rownames(odds) <- colnames(dat1)
                colnames(logPInt) <- colnames(odds) <- colnames(dat2)

                # Thresholding the odd ratios to stand between  -4 and 4
                odds[is.infinite(odds)] <- NA
                # odds[odds<1e-3] = 1e-4
                # odds[odds>1e3] = 1e4

                # odds[odds<1e-2] = 1e-3
                # odds[odds>1e2] = 1e3

                odds[odds<1e-1] = 1e-2
                odds[odds>1e1] = 1e2

                odds[10^-abs(logPInt) > 0.05] = 1

                ## 
                logOdds=log10(odds)
                M <-  matrix( NA, ncol=ncol(odds), nrow=nrow(odds))
                # TODO: More breaks
                # M <- cut(logOdds, breaks = c(-4:0-.Machine$double.eps,0:4), include.lowest=TRUE)  
                # M <- cut(logOdds, breaks = c(-3:0-.Machine$double.eps,0:3), include.lowest=TRUE)  
                M <- cut(logOdds, breaks = c(-2:0-.Machine$double.eps,0:2), include.lowest=TRUE)  
                M <- matrix(as.numeric(M), ncol=ncol(odds), nrow=nrow(odds))

                rownames(M) <- rownames(odds)
                colnames(M) <- colnames(odds)

                # Output
                return(M)
}
