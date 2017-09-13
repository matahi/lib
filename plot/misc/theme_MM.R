theme_MM <- function()
{
        theme(panel.grid=element_blank(),
              text = element_text(size=20),
              legend.position="top",
              panel.background=element_rect(fill="white"),
              axis.text=element_text(colour="black", size=rel(0.8)),
              axis.ticks=element_line(colour="black"),
              panel.border=element_rect(fill=NA, colour="black", size=0.7),
              axis.title.y=element_text(vjust=0.35),
              strip.background=element_rect(colour="black", fill="white",size=0.7),
              axis.line = element_line(colour="black",size=0.7))
}
