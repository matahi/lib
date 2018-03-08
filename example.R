library(survival) # survival package
library(survminer) # survival plots with ggplot
library(ggpubr) # arrange plots for publication
library(ggthemr)

# Colors
ggthemr("fresh")
color.code <- swatch() # ggthemr color code

fit <- survfit(Surv(time, status) ~ x, data=aml)

p.surv <- ggsurvplot(fit, 
                     data = aml,
                     palette = color.code[2:3],
                     pval =T, # add pval
                     pval.coord = c(100, 0.4), # position of pval
                     risk.table=T) # risk.table

#
p.surv

# publication-ready multisurvival plot
ggarrange(p.surv$plot, 
          p.surv$plot,
          labels=c("A.", "B."))

# publication-ready multipanel plot
p.boxplot <- ggplot(iris) + geom_boxplot(aes(x=Species, y=Sepal.Length, fill=Species)) + theme(legend.position="top")

ggarrange(p.surv$plot, 
          p.boxplot,
          labels=c("A.", "B."))




