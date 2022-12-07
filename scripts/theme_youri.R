#!/usr/bin/env R

theme_youri <- theme_bw() +
  theme(
    strip.background = element_rect(colour="white",fill="white"),
    axis.title = element_text(face = "bold",size = rel(1)),
    
    #axis.text.x = element_blank(),
    #axis.ticks.x = element_blank(),
    
    legend.position = 'bottom',
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1.25),
    strip.text = element_text(size = 8) # facet size
  )


theme_youri_barplot <- theme_youri +
  theme(
    # show text at x-axis centered and 90 deg rotated
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    
    # hide xlab, redundant
    axis.title.x = element_blank()
    )


