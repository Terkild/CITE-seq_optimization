fractile.line <- function(x, Q=0.9, trans="log2",add.n=0.5){
  q <- quantile(x,probs=c(Q))
  transrev <- ggforce::trans_reverser(trans)$inverse
  data.frame(y=q, yend=q,label=signif(transrev(-q),2)-add.n)
}

#' Feature-Rankplot
#' 
#' Draws a normalized rank plot for a given feature. 
#' Can include barcode plot to annotate cells along the ranking (such as cell type).
#' 
#' @param data `Double` used for ranking the cells
#' @param group `Factor` used for grouping comparitors (each having independt ranking)
#' @param color `Factor` used for coloring the graphs
#' @param linetype `Factor` used for annotating line and smooth line plots
#' @param wrap `Factor` used for wrapping the plot
#' @param draw.points `Boolean` of whether to draw ranked points
#' @param draw.line `Boolean` of whether to connect ranked points with a line
#' @param draw.smooth `Boolean` of whether to draw a smoothed rank plot line
#' @param draw.fractile `Boolean` of whether to draw a fractile lines and values
#' @param draw.barcode `Boolean` of whether to draw a barcode plot by clusters
#' @param trans Transformation method for visualizing data
#' @param add.n `Double` constant to add for showing "0" in log transformation
#' @param colors Named vector of colors for manual coloring
#' @param fractile.upper `Double` fractile for upper line
#' @param fractile.lower `Double` fractile for lower line
#' @param barcodeGroup `Factor` used grouping the barcode plot
#' @param barcode.stepSize `Double` of the relative size of each barcode "column" (until a better way is determined)
#' @param barcode.downsample=500 `Integer` number of cells to be included in random downsampling within each barcodeGroup
#' @param barcode.alpha `Double` opacity of each "bar" in the barcode plot
#' @param barcode.refGroups `Character` from which groups (from "group" input) should barcodes be shown. (If NULL all groups will be included)
#'
#' @return A ggplot object

feature_rankplot <- function(data, 
                             group=NULL, 
                             color=NULL, 
                             linetype=NULL, 
                             wrap=NULL, 
                             draw.points=TRUE, 
                             draw.line=FALSE, 
                             draw.smooth=FALSE, 
                             draw.fractile=FALSE,
                             draw.barcode=FALSE, 
                             trans="biexp", 
                             add.n=0, 
                             colors=NULL,
                             fractile.upper=0.8, 
                             fractile.lower=0.2, 
                             barcodeGroup=NULL, 
                             barcode.stepSize=0.3, 
                             barcode.downsample=500,
                             barcode.alpha=0.3,
                             barcode.refGroups=NULL,
                             barcode.colors=NULL
                            ){
  require("ggrepel")
  require("ggplot2")
  
  ## Make a combined dataframe before ordering by rank
  data.combined <- data.frame(value=data,include=1)
  if(!is.null(group)) data.combined$group <- group
  if(!is.null(color)) data.combined$color <- color
  if(!is.null(linetype)) data.combined$linetype <- linetype
  if(!is.null(wrap)) data.combined$wrap <- wrap
  if(!is.null(barcodeGroup)) data.combined$barcodeGroup <- barcodeGroup
  
  ## Split data into groups for which ranks should be independent
  if(!is.null(group) | !is.null(wrap)){
    if(is.null(group)){
      data.grouped <- split(data.combined,list(wrap))
    } else if(is.null(wrap)){
      data.grouped <- split(data.combined,list(group))
    } else {
      data.grouped <- split(data.combined,list(group,wrap))
    }
  } else {
    data.grouped <- list(data.combined)
  }
  
  ## Calculate normalized ranks for each group
  rankedList <- lapply(data.grouped,FUN=function(x){
    x <- x[order(x$value),]
    x$rank <- seq_along(x$value)/length(x$value)
    return(x)
  })
  
  ## Define transformation functions
  transFun <- ggforce::trans_reverser(trans)$transform
  transFun.inverse <- ggforce::trans_reverser(trans)$inverse
  
  ## Merge the rankedList into a single dataframe for plotting
  plotData <- do.call("rbind",rankedList)
  plotData$value <- plotData$value+add.n

  p <- ggplot(plotData,aes(x=rank,y=value))
  
  if(draw.smooth == TRUE) p <- p + geom_line(stat="smooth",method="auto",se=FALSE, aes(color=color, group=group, linetype=linetype), alpha=0.5)
  if(draw.points == TRUE) p <- p + geom_point(aes(color=color),alpha=1,pch=19)
  if(draw.line == TRUE) p <- p + geom_line(aes(color=color, group=group, linetype=linetype),alpha=0.5)
  
  ## Include "barcode plot"
  if(draw.barcode == TRUE & !is.null(barcodeGroup)){
    value.max <- max(plotData$value)
    
    # Would be nice to have the barcode plot as a seperate plot
    # But at the same time, we would like to keep the ability to
    # make facets and keep the barcode information for each facet. 
    # To achieve this, we need to define a set "step" size for each
    # barcode line. As this varies depending on the range of values
    # we are using a bit of a hack (seems to work for our data):
    ## Step denotes the "column" width of each barcode column.
    step <- barcode.stepSize*(log(max(plotData$value))/log(300))
    
    ## Where should the first barcode column start (x-axis value)
    step.max <- transFun(value.max)+(step/2)
    barcodeGroups <- unique(plotData$barcodeGroup)
    
    ## Transform steps according to transformation function
    steps <- transFun.inverse((step.max-step+seq_along(barcodeGroups)*step))
    names(steps) <- barcodeGroups
    
    subset <- lapply(barcodeGroups,FUN=function(x){
      ## We allow basing barcodes on specific groups as this allows us to only show
      ## cells from DF1 (as these are most likely to show expression - if detectable)
      if(!is.null(barcode.refGroups)){
        subset <- which(plotData$barcodeGroup == x & plotData$group %in% barcode.refGroups)
      } else {
        subset <- which(plotData$barcodeGroup == x)
      }
      
      ## Barcodes easily saturate with, to avoid this, we do random downsampling within each barcodeGroup
      downsample <- ifelse(length(subset) > barcode.downsample, barcode.downsample, length(subset))
      subset <- subset[sample(x=length(subset), size=downsample, replace=FALSE)]
    })
    subset <- do.call("c",subset)
    plotData.barcode <- plotData[subset,]
    
    ## Allow another color scale
    p <- p + ggnewscale::new_scale_color()
    p <- p + geom_point(data=plotData.barcode, aes(y=steps[barcodeGroup], col=barcodeGroup, alpha=barcode.alpha),shape="-", size=2)
    
    ## Set manual color scheme for barcode groups
    if(!is.null(barcode.colors)){
      p <- p + scale_color_manual(values=barcode.colors)
    }
  }
  
  ## ADD fractile STATS
  if(draw.fractile == TRUE){
    if(fractile.upper > 0){
      ## Add line segments for upper fractile
      # a bit of a hack to get positions to align? is there a better solution?
      p <- p + stat_summary_bin(geom = "segment", binwidth=2, fun.data = fractile.line,  fun.args=list(Q=fractile.upper,trans=trans,add.n=add.n), aes(x=1, xend=fractile.upper, group=group), linetype="dashed")
      
      ## Add text labels for upper fractile
      p <- p + stat_summary_bin(geom = "text_repel", binwidth=2, fun.data = fractile.line,  fun.args=list(Q=fractile.upper,trans=trans,add.n=add.n), aes(x=fractile.upper, group=group), position=position_nudge(x=(1-fractile.upper)), col="black", direction="y",hjust=1,nudge_x=fractile.upper, fontface="bold",segment.alpha=0.25)
      
      ## a bit of a hack to get the lines in all facets - not sure why its needed?
      #,purpose=unique(plotData[,wrap.by])
      p <- p + geom_vline(data=data.frame(expand.grid(list(q=c(fractile.lower,fractile.upper)))),aes(xintercept=q),alpha=0.35,linetype="dotted")
    }
    
    if(fractile.lower > 0){
      ## Add line segments for lower fractile
      p <- p + stat_summary_bin(geom = "segment", binwidth=2, fun.data = fractile.line,  fun.args=list(Q=fractile.lower,trans=trans,add.n=add.n), aes(x=1, xend=fractile.lower, group=group),linetype="dashed",alpha=0.5)
      
      
      ## Add text labels for lower fractile
      p <- p + stat_summary_bin(geom = "text", binwidth=2, fun.data = fractile.line,  fun.args=list(Q=fractile.lower,trans=trans,add.n=add.n), aes(x=0.95, group=group), col="black", hjust=0, vjust=-0.5, fontface="italic")
    }
  }
  
  ## Scale
  if(trans == "biexp"){
    p <- p + scale_y_continuous(trans=trans, limits=c(-1,max(plotData$value)), expand=c(0.01,0.01))
  } else {
    p <- p + scale_y_continuous(trans=trans, expand=c(0.01,0.01))
  }
  
  p <- p + scale_x_continuous(expand=c(0.01,0.01))
  
  ## Facet
  if(!is.null(wrap)){
    p <- p + facet_grid(~wrap)
  }
  
  ## Layout
  p <- p + labs(col="Sample") + theme_bw() + ylab("Count") + xlab("Rank fraction") + guides(alpha=FALSE)
  p <- p + theme_get() + coord_flip()
  
  ## Manual colors
  if(!is.null(colors))p <- p + scale_color_manual(values=colors)
  
  return(p)
}