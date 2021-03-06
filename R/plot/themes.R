# blank_theme
blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )

# heatmap_theme
theme_heatmap <- theme(axis.text.x = element_blank(),
                      axis.ticks = element_blank(),
                      panel.border = element_blank(),
                      panel.grid = element_blank(),
                      axis.line = element_blank())

