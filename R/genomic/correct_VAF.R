correct_VAF <- function(dat.df, CN.effect=data.frame())
{
        source("./src/lib/genomic/calculate_VAF_95CI.R")
        # dat.df <- mutation.df
        # CN.effect <- CN.effect

        out <- Reduce("rbind", lapply(1:nrow(dat.df), function(n)
                                      {
                                              # print(n)
                                              Sample <- as.character(dat.df[n,"ID"])
                                              VAF <- dat.df[n,"VAF"]
                                              Gene <- as.character(dat.df[n,"Gene"])
                                              Chr <- as.character(dat.df[n,"Chr"])
                                              Ntot <- dat.df[n,"TotRead"]
                                              Nmut <- round(VAF*Ntot)

                                              # Gender <- dat.df[n,"gender"] # 1 is male 2 is female?
                                              Gender <- 2 # 1 is male 2 is female?

                                              # if (is.na(Gender))
                                              # {
                                              #         Gender <- 2 
                                              # }


                                              # 1- calculate_VAF_95CI
                                              if (is.na(Ntot))
                                              {
                                                      VAF_95_CI_low <- NA
                                                      VAF_95_CI_upp <- NA
                                              } else
                                              {
                                                      VAF_95_CI_low <- calculate_VAF_95CI(Nmut,Ntot)[1]
                                                      VAF_95_CI_upp <- calculate_VAF_95CI(Nmut,Ntot)[2]
                                              }

                                              if (Chr=="X"&Gender==1) # if chr is X and male sample same as deletion
                                              {
                                                      VAF.corrected <- VAF/2
                                                      correction <- "chr_X"

                                                      # Similar case as deletion
                                                      new_LCI <- calculate_VAF_95CI(Nmut,2*Ntot)[1]
                                                      new_UCI <- calculate_VAF_95CI(Nmut,2*Ntot)[2]
                                              } else if (Chr=="Y") # if chr is Y only one copy
                                              {
                                                      VAF.corrected <- VAF/2
                                                      correction <- "chr_Y"

                                                      # Similar case as deletion
                                                      new_LCI <- calculate_VAF_95CI(Nmut,2*Ntot)[1]
                                                      new_UCI <- calculate_VAF_95CI(Nmut,2*Ntot)[2]
                                              } else if (!(Gene %in% colnames(CN.effect))) # if Gene is not in CN.effect just verify LOH
                                              {
                                                      if (is.na(VAF_95_CI_low))
                                                      {
                                                              VAF.corrected <- VAF/2
                                                              correction <- "None"
                                                              new_LCI <- NA
                                                              new_UCI <- NA
                                                      } else if (VAF_95_CI_upp<0.65)# else if (VAF_95_CI_low<0.65)
                                                      {
                                                              VAF.corrected <- VAF
                                                              correction <- "None"
                                                              new_LCI <- VAF_95_CI_low
                                                              new_UCI <- VAF_95_CI_upp

                                                      } else
                                                      {
                                                              VAF.corrected <- VAF/2
                                                              correction <- "LOH"

                                                              # In this case: Normal copy was lost and mut copy was amplified. So we were supposed to have Ntot but Nmut/2
                                                              new_LCI <- calculate_VAF_95CI(round(Nmut/2),Ntot)[1]
                                                              new_UCI <- calculate_VAF_95CI(round(Nmut/2),Ntot)[2]
                                                      }
                                              } else if (Gene %in% colnames(CN.effect))
                                              {
                                                      if (is.na(CN.effect[Sample,Gene]))# Same as if no CN        
                                                      {
                                                              if (is.na(VAF_95_CI_low))
                                                              {
                                                                      VAF.corrected <- VAF/2
                                                                      correction <- "None"
                                                                      new_LCI <- NA
                                                                      new_UCI <- NA
                                                              } else if (VAF_95_CI_upp<0.65)# else if (VAF_95_CI_low<0.65)
                                                              {
                                                                      VAF.corrected <- VAF
                                                                      correction <- "None"
                                                                      new_LCI <- VAF_95_CI_low
                                                                      new_UCI <- VAF_95_CI_upp

                                                              } else
                                                              {
                                                                      VAF.corrected <- VAF/2
                                                                      correction <- "LOH"

                                                                      # In this case: Normal copy was lost and mut copy was amplified. So we were supposed to have Ntot but Nmut/2
                                                                      new_LCI <- calculate_VAF_95CI(round(Nmut/2),Ntot)[1]
                                                                      new_UCI <- calculate_VAF_95CI(round(Nmut/2),Ntot)[2]
                                                              }
                                                      } else if (CN.effect[Sample,Gene]==0) # no CN 
                                                      {
                                                              if (is.na(VAF_95_CI_low))
                                                              {
                                                                      VAF.corrected <- VAF/2
                                                                      correction <- "None"
                                                                      new_LCI <- NA
                                                                      new_UCI <- NA
                                                              } else if (VAF_95_CI_upp<0.65)# else if (VAF_95_CI_low<0.65)
                                                              {
                                                                      VAF.corrected <- VAF
                                                                      correction <- "None"
                                                                      new_LCI <- VAF_95_CI_low
                                                                      new_UCI <- VAF_95_CI_upp

                                                              } else
                                                              {
                                                                      VAF.corrected <- VAF/2
                                                                      correction <- "LOH"

                                                                      # In this case: Normal copy was lost and mut copy was amplified. So we were supposed to have Ntot but Nmut/2
                                                                      new_LCI <- calculate_VAF_95CI(round(Nmut/2),Ntot)[1]
                                                                      new_UCI <- calculate_VAF_95CI(round(Nmut/2),Ntot)[2]
                                                              }
                                                      } else if (CN.effect[Sample,Gene]==-1) # deletion
                                                      {
                                                              VAF.corrected <- VAF/2
                                                              correction <- "Deletion"

                                                              # In this case we were supposed to have depth=2*Ntot (because half was lost) and we detected Nmut
                                                              new_LCI <- calculate_VAF_95CI(Nmut,2*Ntot)[1]
                                                              new_UCI <- calculate_VAF_95CI(Nmut,2*Ntot)[2]

                                                      } else if (CN.effect[Sample,Gene]==1) # amplification
                                                      {
                                                              # First verify if it is CHR.X
                                                              if (is.na(VAF_95_CI_low))
                                                              {
                                                                      VAF.corrected <- VAF
                                                                      correction <- "None"
                                                                      new_LCI <- VAF_95_CI_low
                                                                      new_UCI <- VAF_95_CI_upp
                                                              } else if (VAF_95_CI_upp <0.35) # TODO: See if correct
                                                              {
                                                                      VAF.corrected <- VAF *3/2
                                                                      correction <- "Norm Amplification"

                                                                      # In this case we were supposed to have depth= Ntot*2/3 (one less normal copy)
                                                                      new_LCI <- calculate_VAF_95CI(Nmut,round(Ntot*2/3))[1]
                                                                      new_UCI <- calculate_VAF_95CI(Nmut,round(Ntot*2/3))[2]
                                                              } else
                                                              {

                                                                      VAF.corrected <- VAF/2
                                                                      correction <- "Mut Amplification"

                                                                      # TODO: See if correct
                                                                      new_LCI <- calculate_VAF_95CI(Nmut,2*Ntot)[1]
                                                                      new_UCI <- calculate_VAF_95CI(Nmut,2*Ntot)[2]

                                                              }
                                                      }
                                              }

                                              return(data.frame(VAF.corr=VAF.corrected,LCI.corr=new_LCI, UCI.corr=new_UCI,correction))
                                      }))
        return(out)

}

