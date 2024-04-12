### forest plot ###

library("readxl")
library(forestplot)

o6_forest <- read_excel("/Users/yuchen/Desktop/Cancer/figures/Figures_0824/o6_forest_data.xlsx")
o6_forest$events <- format(o6_forest$events, big.mark=",") 
o6_forest$events[o6_forest$events == "    NA"] = ""
o6_forest$`P value` = round(o6_forest$`P value`, 3)
#o6_forest$estimate = round(o6_forest$estimate, 2)
#o6_forest$low = round(o6_forest$low, 2)
#o6_forest$high = round(o6_forest$high, 2)
o6_forest$`P value`[o6_forest$`P value` < 0.001] = "<0.001"


# fill in the hazard ratios
#o6_forest$`Hazard ratio \r\n(95% CI)` = paste0(round(o6_forest$estimate,2)," (",round(o6_forest$low,2),"-",round(o6_forest$high,2),")")
#o6_forest$`Hazard ratio \r\n(95% CI)` [o6_forest$`Hazard ratio \r\n(95% CI)` == "NA (NA-NA)" ] = ""

## Labels defining subgroups are a little indented!
subgps <- c(2,3,5,6,8,9,11,12,14,15,17,18,20,21,23,24,26,27,29,30,32,33,35,36,38,39,41,42,44,45,47,48,50,51,53,54,56,57,59,60)
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
    gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), 
    gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"),
    gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), 
    gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), 
    gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"),
    gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), 
    gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), 
    gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), 
    gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"),
    gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), 
    gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), 
    gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), 
    gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"), 
    gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"),
    gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"),
    gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"),
    gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D"),
    gpar(col = "white"), gpar(col = "#836F88"), gpar(col = "#8D8D5D")
  ),
  box = list( gpar(fill = "white"), gpar(fill = "white"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), 
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), 
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), 
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), 
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), 
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), 
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), 
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), 
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), 
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"), 
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180"),
              gpar(fill = "white"), gpar(fill = "#E0C7E4"), gpar(fill = "#C6D180")
  ) 
)
tiff(file.path("/Users/yuchen/Desktop/o6_cat.tiff"), units="in", width=8, height=14, res=300)
forestplot(labeltext=tabletext, graph.pos=3, 
           mean=c(NA,NA,o6_forest$estimate), 
           lower=c(NA,NA,o6_forest$low), upper=c(NA,NA,o6_forest$high),
           hrzl_lines=list("3" = gpar(lwd=2, col="#99999922"), 
                        "8" = gpar(lwd=54, lineend="butt", columns=c(1:5), col="#99999922"),
                       "14" = gpar(lwd=54, lineend="butt", columns=c(1:5), col="#99999922"),
                        "20" = gpar(lwd=54, lineend="butt", columns=c(1:5), col="#99999922"),
                       "26" = gpar(lwd=54, lineend="butt", columns=c(1:5), col="#99999922"),
                      "32" = gpar(lwd=54, lineend="butt", columns=c(1:5), col="#99999922"),
                     "38" = gpar(lwd=54, lineend="butt", columns=c(1:5), col="#99999922"),
           "44" = gpar(lwd=54, lineend="butt", columns=c(1:5), col="#99999922"),
           "50" = gpar(lwd=54, lineend="butt", columns=c(1:5), col="#99999922"),
           "56" = gpar(lwd=54, lineend="butt", columns=c(1:5), col="#99999922"),
           "62" = gpar(lwd=54, lineend="butt", columns=c(1:5), col="#99999922")),
           
           txt_gp=fpTxtGp(label=gpar(fontfamily="Times", cex=1),
                        ticks=gpar(cex=1),
                        xlab=gpar(cex = 1)),
           shapes_gp = styles,
           is.summary = c(T,T,T,F,F,T,F,F,T,F,F,T,F,F,T,F,F,T,F,F,T,F,F,T,F,F,T,F,F,T,F,F,T,F,F,T,F,F,T,F,F,T,F,F,T,F,F,
                          T,F,F,T,F,F,T,F,F,T,F,F,T,F,F),
           #col=fpColors(box="black", lines="black", zero = "gray50"),
           clip = c(0.84, 1.05),
           xticks = c(0.85, 0.90, 0.95, 1, 1.05),
           zero=1, cex=1.2, lineheight = "auto", boxsize=0.55, 
           lwd.ci=2, ci.vertices=TRUE, ci.vertices.height = 0.2
)

dev.off()


##### Additional adjusted models
o6_forest <- read_excel("/Users/yuchen/Desktop/Cancer/figures/Figures_0824/o6_add.xlsx")
o6_forest$events <- format(o6_forest$events, big.mark=",") 
o6_forest$`P value` = round(o6_forest$`P value`, 3)
o6_forest$estimate = round(o6_forest$estimate, 2)
#o6_forest$low = round(o6_forest$low, 2)
#o6_forest$high = round(o6_forest$high, 2)
#o6_forest$`P value for \r\ntrend`[o6_forest$`P value for \r\ntrend` < 0.001] = "<0.001"
o6_forest$`P value`[o6_forest$`P value` < 0.001] = "<0.001"


## The rest of the columns in the table. 
tabletext <- cbind(c("Cancer type",o6_forest$`Cancer type`), 
                   c("Events",o6_forest$events), 
                   c("HR",o6_forest$`HR`), 
                   c("P value",o6_forest$`P value`))

tiff(file.path("/Users/yuchen/Desktop/o6_add.tiff"), units="in", width=7.5, height=6, res=300)
forestplot(labeltext=tabletext, graph.pos=3, 
           mean=c(NA,o6_forest$estimate), 
           lower=c(NA,o6_forest$low), upper=c(NA,o6_forest$high),
           txt_gp=fpTxtGp(label=gpar(fontfamily="Times", cex=1.1),
                          ticks=gpar(cex=1),
                          xlab=gpar(cex = 1)),
           #hrzl_lines=list("2" = gpar(lwd=2, col="#99999922")),
          # shapes_gp = styles,
           is.summary = c(T,F,F,F,F,F,F,F,F,F,F,F),
           hrzl_lines=list("7" = gpar(lwd=440, lineend="butt", columns=1, col="#99999922")),
           col=fpColors(box="#E0C7E4", lines="#836F88", zero = "gray50"),
           clip = c(0.89, 1.05),
           xticks = c(0.9, 0.95, 1.0, 1.05),
           zero=1, cex=1.2, lineheight = "auto", boxsize=0.33, 
           lwd.ci=1.6, ci.vertices=TRUE, ci.vertices.height = 0.1
)

dev.off()


##### Additional adjusted models for o3
o6_forest <- read_excel("/Users/yuchen/Desktop/Cancer/figures/Figures_0824/o3_add.xlsx")
o6_forest$events <- format(o6_forest$events, big.mark=",") 
o6_forest$`P value` = round(o6_forest$`P value`, 3)
o6_forest$estimate = round(o6_forest$estimate, 2)
#o6_forest$low = round(o6_forest$low, 2)
#o6_forest$high = round(o6_forest$high, 2)
#o6_forest$`P value for \r\ntrend`[o6_forest$`P value for \r\ntrend` < 0.001] = "<0.001"
o6_forest$`P value`[o6_forest$`P value` < 0.001] = "<0.001"


## The rest of the columns in the table. 
tabletext <- cbind(c("Cancer type",o6_forest$`Cancer type`), 
                   c("Events",o6_forest$events), 
                   c("HR",o6_forest$`HR`), 
                   c("P value",o6_forest$`P value`))

tiff(file.path("/Users/yuchen/Desktop/o3_add.tiff"), units="in", width=7.5, height=6, res=300)
forestplot(labeltext=tabletext, graph.pos=3, 
           mean=c(NA,o6_forest$estimate), 
           lower=c(NA,o6_forest$low), upper=c(NA,o6_forest$high),
           txt_gp=fpTxtGp(label=gpar(fontfamily="Times", cex=1.1),
                          ticks=gpar(cex=1),
                          xlab=gpar(cex = 1)),
           #hrzl_lines=list("2" = gpar(lwd=2, col="#99999922")),
           # shapes_gp = styles,
           is.summary = c(T,F,F,F,F,F,F,F,F,F,F,F),
           hrzl_lines=list("7" = gpar(lwd=440, lineend="butt", columns=1, col="#99999922")),
           col=fpColors(box="#C6D180", lines="#8D8D5D", zero = "gray50"),
           clip = c(0.89, 1.05),
           xticks = c(0.9, 0.95, 1.0, 1.05),
           zero=1, cex=1.2, lineheight = "auto", boxsize=0.33, 
           lwd.ci=1.6, ci.vertices=TRUE, ci.vertices.height = 0.1
)

dev.off()





