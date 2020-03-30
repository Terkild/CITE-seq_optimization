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
