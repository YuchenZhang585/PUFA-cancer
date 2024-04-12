library("readxl")
library(forestplot)

o6_forest <- read_excel("/Users/yuchen/Desktop/Cancer/figures/Figures_0824/con_forest.xlsx")
o6_forest$events <- format(o6_forest$events, big.mark=",") 
o6_forest$events[o6_forest$events == "    NA"] = ""
o6_forest$`P value` = formatC(o6_forest$`P value`, digits = 3, format = "f")
o6_forest$`P value`[o6_forest$`P value` == "  NA"] = " "


o6_forest$`Cancer type`[o6_forest$`Cancer type` == "Additionally adjusted"] = "Additionally adjusted*"


#o6_forest$estimate = round(o6_forest$estimate, 2)
#o6_forest$low = round(o6_forest$low, 2)
#o6_forest$high = round(o6_forest$high, 2)



# fill in the hazard ratios
#o6_forest$`Hazard ratio \r\n(95% CI)` = paste0(round(o6_forest$estimate,2)," (",round(o6_forest$low,2),"-",round(o6_forest$high,2),")")
#o6_forest$`Hazard ratio \r\n(95% CI)` [o6_forest$`Hazard ratio \r\n(95% CI)` == "NA (NA-NA)" ] = ""

## Labels defining subgroups are a little indented!
subgps <- c(2,3,5,6,8,9,10,12,13,15,16,17,19,20,21,23,24,26,27,28,30,31,32,34,35,36,38,39,41,42,43,45,46,47,49,50,51,53,54,55,57,58,60,61,63,64,66,67,69,70)
o6_forest$`Cancer type`[subgps] <- paste("   ",o6_forest$`Cancer type`[subgps]) 


## The rest of the columns in the table. 
tabletext <- cbind(c("Cancer type","\n",o6_forest$`Cancer type`), 
                   c("Events","\n",o6_forest$events), 
                   c("HR (95% CI)","\n",o6_forest$HR), 
                   c("P value","\n",o6_forest$`P value`))


# draw the plot
# "darkgoldenrod1"

styles <- fpShapesGp(
  lines = list(  gpar(col = "white"),  gpar(col = "white"), 
                 gpar(lwd=10, col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"),
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"),
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), gpar(col = "darkorange3"), 
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"),
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), gpar(col = "darkorange3"), 
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), gpar(col = "darkorange3"), 
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"),
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), gpar(col = "darkorange3"), 
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), gpar(col = "darkorange3"), 
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), gpar(col = "darkorange3"), 
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"),
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), gpar(col = "darkorange3"), 
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), gpar(col = "darkorange3"), 
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), gpar(col = "darkorange3"), 
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), gpar(col = "darkorange3"), 
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"),
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"),
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"),
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"),
                 gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D")
  ),
  box = list( gpar(fill = "white"), gpar(fill = "white"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), gpar(fill = "darkorange"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), gpar(fill = "darkorange"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), gpar(fill = "darkorange"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), gpar(fill = "darkorange"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), gpar(fill = "darkorange"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), gpar(fill = "darkorange"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), gpar(fill = "darkorange"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), gpar(fill = "darkorange"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), gpar(fill = "darkorange"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), gpar(fill = "darkorange"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180")
  ) 
)
tiff(file.path("/Users/yuchen/Desktop/or.tiff"), units="in", width=10, height=18, res=400)
forestplot(labeltext=tabletext, graph.pos=3, 
           mean=c(NA,NA,o6_forest$estimate), 
           lower=c(NA,NA,o6_forest$low), upper=c(NA,NA,o6_forest$high),
           txt_gp=fpTxtGp(label=gpar(fontfamily="Times", cex=1),
                          ticks=gpar(cex=1),
                          xlab=gpar(cex = 1)),
           shapes_gp = styles,
           is.summary = c(T,T,T,F,F,T,F,F,T,F,F,F,T,F,F,T,F,F,F,T,F,F,F,T,F,F,T,F,F,F,T,F,F,F,T,F,F,F,T,F,F,T,F,F,F,T,F,F,F,T,F,F,F,T,F,F,F,
                          T,F,F,T,F,F,T,F,F,T,F,F,T,F,F),
           #col=fpColors(box="black", lines="black", zero = "gray50"),
           clip = c(0.9, 1.05),
           xticks = c(0.90, 0.95, 1, 1.05),
           zero=1, cex=1.2, lineheight = "auto", boxsize=0.55, 
           lwd.ci=2, ci.vertices=TRUE, ci.vertices.height = 0.2
)

dev.off()
