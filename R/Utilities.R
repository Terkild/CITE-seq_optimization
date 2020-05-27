#' Knee plot for filtering empty droplets
#' 
#' Visualizes the inflection point to filter empty droplets. This function plots 
#' different datasets with a different color. 
#' 
#' @param bc_rank A `DataFrame` output from `DropletUtil::barcodeRanks`.
#' @return A ggplot2 object.
#' 
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
                  rank_cutoff = length(which(bc_rank$total > S4Vectors::metadata(bc_rank)[["inflection"]])),
                  knee = S4Vectors::metadata(bc_rank)[["knee"]],
                  knee_cutoff = length(which(bc_rank$total > S4Vectors::metadata(bc_rank)[["knee"]])))
  p <- ggplot(knee_plt, aes(total, rank)) +
    geom_line() +
    geom_ribbon(aes(xmin = 0, xmax = total, fill = rank > annot$rank_cutoff), alpha=0.5) + 
    geom_hline(data=annot,aes(yintercept = rank_cutoff), linetype = 2) +
    geom_label(data=annot,aes(y=rank_cutoff,x=Inf,label=rank_cutoff), hjust=1, vjust=1) + 
    scale_fill_manual(values=c("black","grey"), labels=c("Cell","EmptyDrop")) + 
    scale_x_log10(expand=c(0,0,0.05,0)) +
    scale_y_log10(expand=c(0,0,0.05,0)) +
    annotation_logticks() +
    labs(y = "Rank", x = "Total UMIs") + 
    guides(fill=guide_legend(override.aes=list(alpha=1, color="black"))) + 
    theme(legend.position=c(1,1), 
          legend.justification=c(1,1), 
          legend.title=element_blank(),
          legend.direction="vertical",
          legend.key.size=unit(0.3,"cm"),
          legend.background=element_blank())
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
    geom_line(color="black") +    
    geom_hline(yintercept=length(highlight), linetype="dashed", color="red", size=0.25, alpha=0.5) + 
    geom_label(data=annot,aes(y=length(highlight),x=Inf,label=length(highlight)), hjust=1, vjust=1) + 
    scale_x_log10(expand=c(0,0,0.05,0)) +
    scale_y_log10(expand=c(0,0,0.05,0)) + 
    annotation_logticks() +
    labs(y = "Rank", x = "Total UMIs") + 
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

## Biexponential transformation (inspired by flowJo)
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