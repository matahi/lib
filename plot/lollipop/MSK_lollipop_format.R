MSK_lollipop_format <- function (dat.genetics) {

        old.infos <- c("Gene", "PROTEIN_CHANGE", "type") 
        new.infos <- c("Hugo_Symbol", "Protein_Change", "Mutation_Type") 

        MSK.format <- data.frame(matrix(0, nrow=nrow(dat.genetics), ncol=length(new.infos))) 
        colnames(MSK.format) <- new.infos


        ##########
        # 



        ##
        # At least Hugo_Symbol and Protein_Change have to be defined
        MSK.format[,"Hugo_Symbol"] <- dat.genetics[,"Gene"]

        tmp <- gsub("p\\.","",dat.genetics[,"PROTEIN_CHANGE"])
        MSK.format[,"Protein_Change"] <- gsub("p\\.","",dat.genetics[,"PROTEIN_CHANGE"])

        ##
        #
        
}
