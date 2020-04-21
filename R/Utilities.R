## A bit of a hacked version of the foot plot for fast plotting of individual markers.
foot_plot_density_custom <- function(marker,conc=NULL, gates=NULL, wrap="tissue", group="dilution", color=color.dilution, trans="biexp", legend=TRUE){
  library("cowplot")
    curMarker <- FetchData(object, vars=c(marker,"supercluster",group,wrap), slot = "counts")
    colnames(curMarker)[1:3] <- c("count","supercluster","group")
    
    color.manual <- color
    if(is.null(wrap)){
      curWrap <- NULL
      curMarker$wrap <- 1
    } else {
      colnames(curMarker)[4] <- "wrap"
    }
    
    if(group == "dilution"){
      if(!is.null(conc)){
        curMarker$conc <- conc
        curMarker$conc[curMarker$group == "DF4"] <- conc/4
        curMarker$conc <- factor(curMarker$conc, levels=rev(sort(unique(curMarker$conc))))
        levels(curMarker$conc) <- sprintf("%2.2fug/mL",as.double(levels(curMarker$conc)))
      } else {
        curMarker$conc <- curMarker$group
      }
      
      curMarker$group <- curMarker$conc
      curWrap <- curMarker$wrap
      names(color.manual) <- levels(curMarker$group)
    }
    
    curMarker.sum <- curMarker %>% group_by(wrap=wrap, group=group) %>% summarise(sum=sum(count)) %>% arrange(wrap, sum)
    curMarker.sum.label <- curMarker.sum %>% group_by(wrap) %>% summarise(label=paste(paste0(group,": ",sprintf("%05s",as.character(sum))),collapse="\n"))
    
    p.hist <- ggplot(curMarker, aes(x=count)) + 
      scale_x_continuous(trans=trans,limits=c(-1,max(curMarker$count)), expand=c(0.01,0.01)) + 
      geom_density(aes(y=..density.. ,linetype=group, fill=group), alpha=0.5, bw=0.35) + 
      guides(fill=guide_legend(reverse = TRUE), linetype=guide_legend(reverse = TRUE)) + 
      scale_y_continuous(expand=c(0,0)) + scale_fill_manual(values=color.manual) + 
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
      p.hist <- p.hist + geom_text(data=curMarker.sum.label, x = Inf, y = Inf, hjust=1, vjust=1.5, aes(label=label), size=1.5)
      p.hist <- p.hist + facet_wrap( ~wrap)
    } else {
      scale_label <- with(curMarker.sum[order(factor(curMarker.sum$group, levels=levels(curMarker$group))),],paste0(group,": ",sprintf("%05s",as.character(sum))))
      p.hist <- p.hist + 
        scale_linetype_discrete(labels=scale_label) + 
        scale_fill_manual(values=color.manual, labels=scale_label) + 
        theme(legend.position=c(1,1.5), legend.justification=c(1,1), legend.text.align=1, plot.margin=unit(c(0.3,0,0,0),"cm"))
    }
    #geom_text_repel(data=curMarker.sum, 
    #                aes(label=paste0(group,": ",sprintf("%05s",as.character(sum)))), 
    #                y=Inf, x=Inf, hjust=0, vjust=0, direction="y", segment.alpha=0, seed=114, size=3)
    
    p.foot <- foot_plot_flip(curMarker$count, group=curMarker$group, linetype=curMarker$group, wrap=curWrap, barcodeGroup=curMarker$supercluster, barcode.stepSize=0.4, draw.points = F, draw.barcode = T, draw.line = T, barcode.refGroups=levels(curMarker$group)[1], trans=trans) + 
      theme(strip.text = element_blank(), 
            plot.margin = unit(c(0,0.3,0,0),"cm"), 
            legend.direction = "vertical",
            legend.position = c(0.90,0.02), 
            legend.key.size=unit(0.2,"cm"),
            legend.title=element_blank(),
            legend.justification=c(1,0),
            axis.title=element_blank(),
            axis.text.y=element_blank()) + 
      guides(linetype=F, col=guide_legend(override.aes = list(shape = 15)), group=F) + ylab("UMI count")
    #geom_text_repel(data=curMarker.sum, aes(label=paste0(group,": ",sprintf("%06s",as.character(sum)))), y=-ggforce::trans_reverser(trans)$transform(max(curMarker$count)*0.5), x=-Inf, xlim=c(0.02,1), hjust=1, vjust=-1, direction="y", segment.alpha=0, seed=114, size=3)
    
    
    if(legend == FALSE){
      p.foot <- p.foot + theme(legend.position="none")
    }
    
    if(!is.null(gates)){
      p.foot <- p.foot + geom_vline(data=gates,aes(xintercept=gate), col="red", alpha=0.5, linetype="dashed")
    }
    
    p.hist.legend <- get_legend(p.hist)
    #p.hist <- p.hist + theme(legend.position="none")
    p.foot.legend <- get_legend(p.foot)
    #p.foot <- p.foot + theme(legend.position="none")
    #p.legends <- plot_grid(p.hist.legend, p.foot.legend, nrow=1)
    p <- plot_grid(p.hist, p.foot, ncol=1, align="v", axis="lr", label_size=7, labels=c(curMarker.name,""), hjust = 0, vjust=1.1, rel_heights=c(5,15,2))
    
    return(p)
}


biexp_trans <- function(lim = 5, decade.size = lim){
    trans <- function(x){
        ifelse(x <= lim,
               x,
               lim + decade.size * (suppressWarnings(log(x, 10)) -
                                        log(lim, 10)))
    }
    inv <- function(x) {
        ifelse(x <= lim,
               x,
               10^(((x-lim)/decade.size) + log(lim,10)))
    }
    breaks <- function(x) {
        if (all(x <= lim)) {
            scales::pretty_breaks()(x)
        } else if (all(x > lim)) {
            scales::breaks_log(10)(x)
        } else {
            unique(c(scales::pretty_breaks()(c(x[1],lim)),
                     scales::breaks_log(10)(c(lim, x[2]))))
        }
    }
    scales::trans_new(paste0("biexp-",format(lim)), trans, inv, breaks)
}

biexp1000_trans <- function(lim = 800, decade.size = lim){
    trans <- function(x){
        ifelse(x <= lim,
               x,
               lim + decade.size * (suppressWarnings(log(x, 10)) -
                                        log(lim, 10)))
    }
    inv <- function(x) {
        ifelse(x <= lim,
               x,
               10^(((x-lim)/decade.size) + log(lim,10)))
    }
    breaks <- function(x) {
        if (all(x <= lim)) {
            scales::pretty_breaks()(x)
        } else if (all(x > lim)) {
            scales::breaks_log(10)(x)
        } else {
            unique(c(scales::pretty_breaks()(c(x[1],lim)),
                     scales::breaks_log(10)(c(lim, x[2]))))
        }
    }
    scales::trans_new(paste0("biexp-",format(lim)), trans, inv, breaks)
}

makeBiexp <- function(x){
    x <- x + scale_x_continuous(trans="biexp")
    x
}
themeTBB <- function(x){
    x <- x + guides(fill=F)
    x <- x + theme_bw() + theme(axis.ticks.x=element_line(),axis.text.x=element_text(angle=45,hjust=1),axis.title.y=element_blank())
    x
}
ridgeTBB <- function(p){
    o <- lapply(p,FUN=function(x)x+theme(axis.text.y=element_blank(),axis.title.y=element_blank()))
    o[[1]] <- p[[1]]
    CombinePlots(o)
}
UtilityLoaded <- TRUE

# Slightly modified from BUSpaRse, just to avoid installing a few dependencies not used here
read_count_output <- function(dir, name) {
    dir <- normalizePath(dir, mustWork = TRUE)
    m <- Matrix::readMM(paste0(dir, "/", name, ".mtx"))
    m <- Matrix::t(m)
    m <- as(m, "dgCMatrix")
    # The matrix read has cells in rows
    ge <- ".genes.txt"
    genes <- readLines(file(paste0(dir, "/", name, ge)))
    barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
    colnames(m) <- barcodes
    rownames(m) <- genes
    return(m)
}
#' Knee plot for filtering empty droplets
#' 
#' Visualizes the inflection point to filter empty droplets. This function plots 
#' different datasets with a different color. Facets can be added after calling
#' this function with `facet_*` functions. Will be added to the next release
#' version of BUSpaRse.
#' 
#' @param bc_rank A `DataFrame` output from `DropletUtil::barcodeRanks`.
#' @return A ggplot2 object.
knee_plot <- function(bc_rank) {
    library("ggplot2")
    knee_plt <- tibble(rank = bc_rank[["rank"]],
                       total = bc_rank[["total"]]) %>% 
        distinct() %>% 
        dplyr::filter(total > 0)
    annot <- tibble(inflection = S4Vectors::metadata(bc_rank)[["inflection"]],
                    rank_cutoff = max(bc_rank$rank[bc_rank$total > S4Vectors::metadata(bc_rank)[["inflection"]]]),
                    knee = S4Vectors::metadata(bc_rank)[["knee"]],
                    knee_cutoff = max(bc_rank$rank[bc_rank$total > S4Vectors::metadata(bc_rank)[["knee"]]]))
    p <- ggplot(knee_plt, aes(total, rank)) +
        geom_line() +
        geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2) +
        geom_vline(aes(xintercept = inflection), data = annot, linetype = 2) +
        geom_hline(aes(yintercept = knee_cutoff), data = annot, linetype = "dotted", col="red") +
        geom_vline(aes(xintercept = knee), data = annot, linetype = "dotted", col="red") +
        geom_label(aes(y=0,x=knee,label=knee), data = annot, hjust=0, vjust=0, col="red") + 
        geom_label(aes(y=0,x=inflection,label=inflection), data = annot, hjust=1, vjust=0) + 
        geom_label(aes(y=knee_cutoff,x=Inf,label=knee_cutoff), data = annot, hjust=1, vjust=1, col="red") + 
        geom_label(aes(y=rank_cutoff,x=Inf,label=rank_cutoff), data = annot, hjust=1, vjust=0) + 
        scale_x_log10() +
        scale_y_log10() +
        annotation_logticks() +
        labs(y = "Rank", x = "Total UMIs")
    return(p)
}

knee_plot_auc <- function(bc_rank) {
  library("ggplot2")
  knee_plt <- tibble(rank = bc_rank[["rank"]],
                     total = bc_rank[["total"]]) %>% 
    distinct() %>% 
    dplyr::filter(total > 0)
  annot <- tibble(inflection = S4Vectors::metadata(bc_rank)[["inflection"]],
                  rank_cutoff = max(bc_rank$rank[bc_rank$total > S4Vectors::metadata(bc_rank)[["inflection"]]]),
                  knee = S4Vectors::metadata(bc_rank)[["knee"]],
                  knee_cutoff = max(bc_rank$rank[bc_rank$total > S4Vectors::metadata(bc_rank)[["knee"]]]))
  p <- ggplot(knee_plt, aes(total, rank)) +
    geom_line() +
    geom_ribbon(aes(xmin = 0, xmax = total, fill = rank > annot$rank_cutoff), alpha=0.5) + 
    geom_hline(data=annot,aes(yintercept = rank_cutoff), linetype = 2) +
    geom_label(data=annot,aes(y=rank_cutoff,x=Inf,label=rank_cutoff), hjust=1, vjust=0) + 
    scale_fill_manual(values=c("black","grey"), labels=c("Cell","EmptyDrop")) + 
    scale_x_log10(expand=c(0,0,0.05,0)) +
    scale_y_log10(expand=c(0,0,0.05,0)) +
    annotation_logticks() +
    labs(y = "Rank", x = "Total UMIs") + 
    guides(fill=guide_legend(override.aes=list(alpha=1, color="black"))) + 
    theme(legend.position=c(1,.99), 
          legend.justification=c(1,1), 
          legend.title=element_blank(),
          legend.direction="vertical")
  return(p)
}

knee_plot_highlight <- function(bc_rank, highlight=c()) {
  library("ggplot2")
  knee_plt <- tibble(rank = bc_rank[["rank"]],
                     total = bc_rank[["total"]], 
                     barcode=rownames(bc_rank)) %>% 
    distinct() %>% 
    dplyr::filter(total > 0)
  
  annot <- tibble(inflection = S4Vectors::metadata(bc_rank)[["inflection"]],
                  rank_cutoff = max(bc_rank$rank[bc_rank$total > S4Vectors::metadata(bc_rank)[["inflection"]]]),
                  knee = S4Vectors::metadata(bc_rank)[["knee"]],
                  knee_cutoff = max(bc_rank$rank[bc_rank$total > S4Vectors::metadata(bc_rank)[["knee"]]]))
  
  cutoff <- 18000
  data.highlight <- knee_plt[knee_plt$barcode %in% highlight,]
  data.highlight <- rbind(data.highlight[data.highlight$rank <= cutoff,],data.highlight[sample(nrow(data.highlight[data.highlight$rank > cutoff,]),1000),])
  
  p <- ggplot(knee_plt, aes(total, rank)) +
    #geom_point(data=data.highlight,aes(x=1), shape="-", size=1, alpha=.2) + 
    geom_segment(data=data.highlight[data.highlight$rank > cutoff,],aes(x=0, xend=total, yend=rank), size=0.001, color="black", alpha=1, show.legend=FALSE) + 
    geom_rect(data=data.highlight[data.highlight$rank <= cutoff,],aes(xmin=0, xmax=total, ymin=rank-0.5, ymax=rank+0.5), fill="black", alpha=1) + 
    geom_line(color="grey") +    
    #geom_ribbon(data=data.highlight, aes(xmin = 0, xmax = total), fill="black", alpha=0.5) +
    #geom_segment(data=knee_plt[-c(knee_plt$barcode %in% highlight | knee_plt$rank > max(data.highlight$rank)),],aes(x=0, xend=total, yend=rank, size=(10/(10^rank))), color="white", alpha=0.05, show.legend=FALSE) + 
    # to show emptyDroplets within dense cell area
    #geom_segment(data=knee_plt[-c(knee_plt$rank > max(highlight$rank) | knee_plt$barcode %in% highlight),],aes(x=0, xend=total, yend=rank), color="white", alpha=0.1, size=0.05) +
    #geom_hline(data=annot,aes(yintercept = rank_cutoff), linetype = 2) +
    #geom_label(data=annot,aes(y=rank_cutoff,x=Inf,label=rank_cutoff), hjust=1, vjust=0) + 
    geom_hline(yintercept=length(highlight), linetype="dashed", color="red", size=0.25, alpha=0.5) + 
    geom_label(data=annot,aes(y=length(highlight),x=Inf,label=length(highlight)), hjust=1, vjust=1) + 
    scale_x_log10(expand=c(0,0,0.05,0)) +
    scale_y_log10(expand=c(0,0,0.05,0)) + 
    #scale_color_manual(values=c("black","white")) + 
    annotation_logticks() +
    labs(y = "Rank", x = "Total UMIs") + 
    #guides(color=guide_legend(override.aes=list(alpha=1, size=1))) + 
    theme(legend.position=c(1,.99), 
          legend.justification=c(1,1), 
          legend.title=element_blank(),
          legend.direction="vertical")
  return(p)
}

## nth function extracts the value at a set fractile or median if fractile "rank" is less than a set "nth" threshhold
nth <- function(value, nth=10, fractile=0.9){
  if(length(value)*(1-fractile) <= nth){
    newvalue <- median(value)
  } else {
    newvalue <- quantile(value, probs=c(fractile))
  }
  return(newvalue)
}