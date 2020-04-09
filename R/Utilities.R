## A bit of a hacked version of the foot plot for fast plotting of individual markers.
foot_plot_density_custom <- function(marker,conc=NULL, gates=NULL, trans="biexp", legend=TRUE){
  library("cowplot")
    curMarker <- FetchData(object, vars=c(marker,"tissue","dilution", "supercluster"), slot = "counts")
    colnames(curMarker)[1] <- "count"
    if(!is.null(conc)){
        curMarker$conc <- conc
        curMarker$conc[curMarker$dilution == "DF4"] <- conc/4
        curMarker$conc <- factor(curMarker$conc, levels=rev(sort(unique(curMarker$conc))))
        levels(curMarker$conc) <- sprintf("%2.2fµg/mL",as.double(levels(curMarker$conc)))
    } else {
        curMarker$conc <- curMarker$dilution
    }
    
    color.manual <- color.dilution
    names(color.manual) <- levels(curMarker$conc)
    
    curMarker.sum <- curMarker %>% group_by(tissue=tissue, conc=conc) %>% summarise(sum=sum(count)) %>% arrange(tissue, sum)
    curMarker.sum.label <- curMarker.sum %>% group_by(tissue) %>% summarise(label=paste(paste0(conc,": ",sprintf("%05s",as.character(sum))),collapse="\n"))
    
    p.hist <- ggplot(curMarker, aes(x=count)) + 
      scale_x_continuous(trans=trans,limits=c(-1,max(curMarker$count)), expand=c(0.01,0.01)) + 
      geom_density(aes(y=..density.. ,linetype=conc, fill=conc), alpha=0.5, bw=0.35) + 
      guides(fill=guide_legend(reverse = TRUE), linetype=guide_legend(reverse = TRUE)) + 
      scale_fill_manual(values=color.manual) + scale_y_continuous(expand=c(0,0)) + facet_wrap( ~tissue) + 
      theme(axis.title=element_blank(),
            axis.text.y=element_blank(), 
            axis.text.x=element_blank(), 
            axis.ticks=element_blank(), 
            panel.border=element_blank(), 
            panel.grid = element_blank(), 
            plot.margin = unit(c(0,0,0,0),"cm"), 
            legend.direction="vertical", 
            legend.title=element_blank(),
            legend.background = element_blank(),
            legend.box.margin = unit(c(0.3,0,0,0),"cm"), 
            legend.key.size=unit(0.3,"cm"),
            legend.position=c(0.5,1.5), 
            legend.justification=c(0,1)) + 
      geom_text(data=curMarker.sum.label, x = Inf, y = Inf, hjust=1, vjust=1.5, aes(label=label), size=3)
      #geom_text_repel(data=curMarker.sum, 
      #                aes(label=paste0(group,": ",sprintf("%05s",as.character(sum)))), 
      #                y=Inf, x=Inf, hjust=0, vjust=0, direction="y", segment.alpha=0, seed=114, size=3)

    p.foot <- foot_plot_flip(curMarker$count, group=curMarker$conc, linetype=curMarker$conc, wrap=curMarker$tissue, barcodeGroup=curMarker$supercluster, barcode.stepSize=0.4, draw.points = F, draw.barcode = T, draw.line = T, barcode.refGroups=levels(curMarker$conc)[1], trans=trans) + 
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
        p.foot <- p.foot + geom_vline(data=gates,aes(xintercept=gate), col="red", size=0.25, alpha=0.5, linetype="dashed")
    }

    p.hist.legend <- get_legend(p.hist)
    #p.hist <- p.hist + theme(legend.position="none")
    p.foot.legend <- get_legend(p.foot)
    #p.foot <- p.foot + theme(legend.position="none")
    #p.legends <- plot_grid(p.hist.legend, p.foot.legend, nrow=1)
    p <- plot_grid(p.hist, p.foot, ncol=1, align="v", axis="lr", labels=c(curMarker.name,""), hjust = 0, vjust=1, rel_heights=c(5,15,2))
    
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
