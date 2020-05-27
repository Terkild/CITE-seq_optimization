require("ggplot2")

text.size <- 7
text.axis.size <- 6
panel.label_size <- 10
panel.label_vjust <- 0.98
panel.label_hjust <- 0
figure.resolution <- 600
figure.antialias <- "cleartype"
figure.width.full <- 7
figure.unit <- "in"


theme_set(theme_bw(base_size=text.size) + 
            theme(
              text=element_text(size=text.size),
              axis.text.y=element_text(size=text.axis.size),
              axis.text.x=element_text(angle=45, hjust=1, size=text.axis.size), 
              panel.grid.minor=element_blank(), 
              strip.background=element_blank(), 
              strip.text=element_text(face="bold", size=text.size), 
              legend.position = "bottom",
              plot.margin = unit(c(1,1,1,1),"mm")))

update_geom_defaults("line", list(size=0.35))
update_geom_defaults("bar", list(size=0.25))
update_geom_defaults("tile", list(size=0.25))
update_geom_defaults("rect", list(size=0.25))
update_geom_defaults("density", list(size=0.25))
update_geom_defaults("vline", list(size=0.25))
update_geom_defaults("hline", list(size=0.25))
update_geom_defaults("point", list(size=1))

library("ggalluvial")
update_geom_defaults("stratum", list(size=0.25))
update_geom_defaults("flow", list(size=0.25))
