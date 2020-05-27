feature_rankplot_hist <- function(data, 
                                  group=NULL, 
                                  color=NULL, 
                                  linetype=NULL, 
                                  wrap=NULL, 
                                  barcodeGroup=NULL, 
                                  draw.histogram=TRUE, 
                                  trans="biexp", 
                                  add.n=0, 
                                  histogram.colors=NULL,
                                  title="",
                                  gates=NULL, 
                                  legend=TRUE, 
                                  yaxis.text=FALSE, ...){
  library("cowplot")
  
  ## Make a combined data.matrix
  data.combined <- data.frame(value=data)
  if(!is.null(group)) data.combined$group <- group
  if(!is.null(color)) data.combined$color <- color
  if(!is.null(linetype)) data.combined$linetype <- linetype
  if(!is.null(wrap)){
    data.combined$wrap <- wrap
  } else {
    data.combined$wrap <- 1
  }
  if(!is.null(barcodeGroup)) data.combined$barcodeGroup <- barcodeGroup
  
  ## Calculate (UMI) sum values
  data.combined.sum <- data.combined %>% 
    group_by(wrap=wrap, group=group) %>% 
    summarise(sum=sum(value)) %>% 
    arrange(wrap, sum)
  
  ## Make "nice" labels with group name and UMI sum for each wrap
  data.combined.sum.label <- data.combined.sum %>% 
    group_by(wrap) %>% 
    summarise(label=paste(paste0(group,": ",sprintf("%05s",as.character(sum))),collapse="\n"))

  ## Make histograms  
  if(draw.histogram == TRUE){
    p.hist <- ggplot(data.combined, aes(x=value)) + 
      scale_x_continuous(trans=trans,limits=c(-1,max(data.combined$value)), expand=c(0.01,0.01)) + 
      geom_density(aes(y=..density.. ,linetype=group, fill=group), alpha=0.5, bw=0.35) + 
      guides(fill=guide_legend(reverse = TRUE), linetype=guide_legend(reverse = TRUE)) + 
      scale_y_continuous(expand=c(0,0)) + scale_fill_manual(values=histogram.colors) + 
      theme(axis.title=element_blank(),
            axis.text.y=element_blank(), 
            axis.text.x=element_blank(), 
            axis.ticks=element_blank(), 
            panel.border=element_blank(), 
            panel.grid=element_blank(), 
            plot.margin=unit(c(0,0,0,0),"cm"), 
            legend.direction="vertical", 
            legend.title=element_blank(),
            legend.background=element_blank(),
            legend.box.margin=unit(c(0,0,0,0),"mm"), 
            legend.key.width=unit(0.15,"cm"),
            legend.key.height=unit(0.10,"cm"),
            legend.position=c(0.4,2), 
            legend.justification=c(0,1))
    
    if(!is.null(wrap)){
      p.hist <- p.hist + geom_text(data=data.combined.sum.label, x = Inf, y = Inf, hjust=1, vjust=1.5, aes(label=label), size=1.5)
      p.hist <- p.hist + facet_wrap( ~wrap)
    } else {
      scale_label <- with(data.combined.sum[order(factor(data.combined.sum$group, levels=levels(data.combined$group))),],paste0(group,": ",sprintf("%05s",as.character(sum))))
      p.hist <- p.hist + 
        scale_linetype_discrete(labels=scale_label) + 
        scale_fill_manual(values=histogram.colors, labels=scale_label) + 
        theme(legend.position=c(1,1.5), legend.justification=c(1,1), legend.text.align=1, plot.margin=unit(c(0.3,0,0,0),"cm"))
    }
  }
  
  ## Draw feature_rankplot
  p.feature_rankplot <- feature_rankplot(data=data.combined$value, 
                                         group=data.combined$group, 
                                         linetype=data.combined$group, 
                                         wrap=data.combined$wrap, 
                                         barcodeGroup=data.combined$barcodeGroup, 
                                         barcode.stepSize=0.4, 
                                         draw.points = F, 
                                         draw.barcode = T, 
                                         draw.line = T, 
                                         trans=trans, ...) + 
  theme(strip.text = element_blank(), 
        plot.margin = unit(c(0,0.3,0,0),"cm"), 
        legend.direction = "vertical",
        legend.position = c(0.90,0.02), 
        legend.key.size=unit(0.2,"cm"),
        legend.title=element_blank(),
        legend.justification=c(1,0),
        axis.title=element_blank(),
        axis.text.y=element_blank()) + 
  guides(linetype=F, col=guide_legend(override.aes = list(shape = 15)), group=F) + ylab("UMI count") + xlab("Cell ranking")

  if(legend == FALSE){
    p.feature_rankplot <- p.feature_rankplot + theme(legend.position="none")
  }
  
  if(yaxis.text == TRUE){
    p.feature_rankplot <- p.feature_rankplot + theme(axis.title.y=element_text(size=6))
  }
  
  if(!is.null(gates)){
    p.feature_rankplot <- p.feature_rankplot + geom_vline(data=gates,aes(xintercept=gate), col="red", alpha=0.5, linetype="dashed")
  }
  
  p <- plot_grid(p.hist, p.feature_rankplot, ncol=1, align="v", axis="lr", label_size=7, labels=c(title,""), hjust = 0, vjust=1.1, rel_heights=c(5,15,2))
  
  return(p)
}