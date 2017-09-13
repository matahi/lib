# ### Add age
# Dat.surv.Age <- cbind(Dat.surv, Age=dat[,"Age_at_diagnosis"])
# f <- paste0("Surv(followUp,status) ~", "Age")
# s <- survfit(as.formula(f), data=Dat.surv.Age)
# c <- coxph(as.formula(f), data= Dat.surv.Age)
# summary(c)
# 
# library(survMisc)
# opt.age <- cutp(c)$Age$Age[1]
# f <- paste0("Surv(followUp,status) ~", "Age>",opt.age)
# s <- survfit(as.formula(f), data=Dat.surv.Age)
# c <- coxph(as.formula(f), data= Dat.surv.Age)
# summary(c)


