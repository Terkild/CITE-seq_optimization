CITE-seq optimization - Antibody concentration titration
================
Terkild Brink Buus
30/3/2020

## Load utilities

Including libraries, plotting and color settings and custom utility
functions

``` r
set.seed(114)
require("Seurat", quietly=T)
require("tidyverse", quietly=T)
library("Matrix", quietly=T)
library("patchwork", quietly=T)

## Load ggplot theme and defaults
source("R/ggplot_settings.R")

## Load helper functions
source("R/Utilities.R")

## Load predefined color schemes
source("R/color.R")

## Load feature_rankplot functions
source("R/feature_rankplot.R")
source("R/feature_rankplot_hist.R")
source("R/feature_rankplot_hist_custom.R")

outdir <- "figures"
data.Seurat <- "data/5P-CITE-seq_Titration.rds"
data.abpanel <- "data/Supplementary_Table_1.xlsx"
data.markerStats <- "data/markerByClusterStats.tsv"

## Make a custom function for formatting the concentration scale
scaleFUNformat <- function(x) sprintf("%.2f", x)
```

## Load Seurat object

Subset to only focus on conditions with 50µl staining volume (thus
comparing DF1 to DF4) in both PBMC and lung samples.

``` r
object <- readRDS(file=data.Seurat)

## Show number of cells from each sample
table(object$group)
```

    ## 
    ## PBMC_50ul_1_1000k PBMC_50ul_4_1000k PBMC_25ul_4_1000k  PBMC_25ul_4_200k 
    ##              1777              1777              1777              1777 
    ##  Lung_50ul_1_500k  Lung_50ul_4_500k           Doublet          Negative 
    ##              1681              1681                 0                 0

``` r
object <- subset(object, subset=volume == "50µl")
object
```

    ## An object of class Seurat 
    ## 33572 features across 6916 samples within 3 assays 
    ## Active assay: RNA.kallisto (33514 features)
    ##  2 other assays present: HTO.kallisto, ADT.kallisto
    ##  3 dimensional reductions calculated: pca, tsne, umap

## Load Ab panel annotation and concentrations

Marker stats is reused in other comparisons and was calculated in the
end of the preprocessing vignette.

``` r
abpanel <- data.frame(readxl::read_excel(data.abpanel))
rownames(abpanel) <- abpanel$Marker

markerStats <- read.table(data.markerStats)
markerStats.PBMC <- markerStats[markerStats$tissue == "PBMC",]
rownames(markerStats) <- paste(markerStats$marker,markerStats$tissue,sep="_")

## Make a ordering vector ordering markers per concentration and total UMI count
marker.order <- markerStats.PBMC$marker[order(markerStats.PBMC$conc_µg_per_mL, markerStats.PBMC$UMItotal, decreasing=TRUE)]

head(abpanel)
```

    ##        Marker Category     Alias   Clone Isotype_Mouse Corresponding_gene
    ## CD103   CD103        B      <NA> BerACT8          IgG1              ITGAE
    ## CD107a CD107a        B     LAMP1    H4A3          IgG1              LAMP1
    ## CD117   CD117        E     C-kit   104D2          IgG1                KIT
    ## CD11b   CD11b        B      <NA>  ICRF44          IgG1              ITGAM
    ## CD123   CD123        E      <NA>     6H6          IgG1              IL3RA
    ## CD127   CD127        E IL7Ralpha  A019D5          IgG1               IL7R
    ##        TotalSeqC_Tag BioLegend_Cat Stock_conc_µg_per_mL conc_µg_per_mL
    ## CD103           0145        350233                  500          1.250
    ## CD107a          0155        328649                  500          2.500
    ## CD117           0061        313243                  500          2.500
    ## CD11b           0161        301359                  500          0.625
    ## CD123           0064        306045                  500          0.500
    ## CD127           0390        351356                  500          1.250
    ##        dilution_1x
    ## CD103          400
    ## CD107a         200
    ## CD117          200
    ## CD11b          800
    ## CD123         1000
    ## CD127          400

``` r
head(markerStats)
```

    ##             marker tissue fineCluster nCells UMIsum   nth median   f90 UMItotal
    ## CD103_PBMC   CD103   PBMC           1    638   1740   5.0      2   5.0     5082
    ## CD103_Lung   CD103   Lung          16    132   7084 187.0      5 187.0    60252
    ## CD107a_PBMC CD107a   PBMC           1    638   7757  26.3      8  26.3    13396
    ## CD107a_Lung CD107a   Lung          12    260  11674  99.2     15  99.2    23273
    ## CD117_PBMC   CD117   PBMC           1    638   1318   4.0      2   4.0     3316
    ## CD117_Lung   CD117   Lung          21     32   1695  41.0     41 130.4     5878
    ##             Category Alias   Clone Isotype_Mouse Corresponding_gene
    ## CD103_PBMC         B  <NA> BerACT8          IgG1              ITGAE
    ## CD103_Lung         B  <NA> BerACT8          IgG1              ITGAE
    ## CD107a_PBMC        B LAMP1    H4A3          IgG1              LAMP1
    ## CD107a_Lung        B LAMP1    H4A3          IgG1              LAMP1
    ## CD117_PBMC         E C-kit   104D2          IgG1                KIT
    ## CD117_Lung         E C-kit   104D2          IgG1                KIT
    ##             TotalSeqC_Tag BioLegend_Cat Stock_conc_µg_per_mL conc_µg_per_mL
    ## CD103_PBMC            145        350233                  500           1.25
    ## CD103_Lung            145        350233                  500           1.25
    ## CD107a_PBMC           155        328649                  500           2.50
    ## CD107a_Lung           155        328649                  500           2.50
    ## CD117_PBMC             61        313243                  500           2.50
    ## CD117_Lung             61        313243                  500           2.50
    ##             dilution_1x marker.y DSB.cutoff positive count   pct
    ## CD103_PBMC          400    CD103          7       14  1777  0.79
    ## CD103_Lung          400    CD103          7      501  1681 29.80
    ## CD107a_PBMC         200   CD107a          7      122  1777  6.87
    ## CD107a_Lung         200   CD107a          7      150  1681  8.92
    ## CD117_PBMC          200    CD117          7        3  1777  0.17
    ## CD117_Lung          200    CD117          7       32  1681  1.90

## Cell type and tissue overview

Make tSNE plots colored by cell type, cluster and tissue of origin.

``` r
p.tsne.tissue <- DimPlot(object, group.by="tissue", reduction="tsne", pt.size=0.1, combine=FALSE)[[1]] + theme_get() + facet_wrap(~"Tissue") + scale_color_manual(values=color.tissue)

p.tsne.cluster <- DimPlot(object, group.by="supercluster", reduction="tsne", pt.size=0.1, combine=FALSE)[[1]] + theme_get() + scale_color_manual(values=color.supercluster) + facet_wrap(~"Cell types")

p.tsne.finecluster <- DimPlot(object, label=TRUE, label.size=3, reduction="tsne", pt.size=0.1, group.by="fineCluster", combine=FALSE)[[1]] + theme_get() + facet_wrap( ~"Clusters") + guides(col=F)

p.tsne.cluster + p.tsne.finecluster + p.tsne.tissue
```

![](Antibody-titration_files/figure-gfm/tsnePlots-1.png)<!-- -->

## Overall ADT counts

Extract UMI data and calculate UMI sum per marker within each condition.

``` r
## Get the data
ADT.matrix <- data.frame(GetAssayData(object, assay="ADT.kallisto", slot="counts"))
ADT.matrix$marker <- rownames(ADT.matrix)
ADT.matrix$conc <- abpanel[ADT.matrix$marker,"conc_µg_per_mL"]
ADT.matrix <- ADT.matrix %>% pivot_longer(c(-marker,-conc))

## Get cell annotations
cell.annotation <- FetchData(object, vars=c("dilution", "tissue"))

## Calculate marker sum from each dilution within both tissues
ADT.matrix.agg <- ADT.matrix %>% group_by(dilution=cell.annotation[name,"dilution"], tissue=cell.annotation[name,"tissue"], marker, conc) %>% summarise(sum=sum(value))

## Order markers by concentration
ADT.matrix.agg$marker.byConc <- factor(ADT.matrix.agg$marker, levels=marker.order)

## Extract marker annotation
ann.markerConc <- abpanel[marker.order,]
ann.markerConc$Marker <- factor(marker.order, levels=marker.order)

ADT.matrix.agg.total <- ADT.matrix.agg
```

## Plot overall ADT counts by conditions

Samples stained with diluted Ab panel have reduced ADT counts.

``` r
p.UMIcountsPerCondition <- ggplot(ADT.matrix.agg.total[order(-ADT.matrix.agg$conc, -ADT.matrix.agg$sum),], aes(x=dilution, y=sum/10^6, fill=conc)) + 
  geom_bar(stat="identity", col=alpha(col="black",alpha=0.05)) + 
  scale_fill_viridis_c(trans="log2", labels=scaleFUNformat, breaks=c(0.0375,0.15,0.625,2.5,10)) + 
  scale_y_continuous(expand=c(0,0,0,0.05)) + 
  labs(fill="DF1\nµg/mL", y=bquote("ADT UMI counts ("~10^6~")")) + 
  guides(fill=guide_colourbar(reverse=T)) + 
  facet_wrap(~tissue) + 
  theme(panel.grid.major=element_blank(), axis.title.x=element_blank(), panel.border=element_blank(), axis.line = element_line(), legend.position="right")

p.UMIcountsPerCondition
```

![](Antibody-titration_files/figure-gfm/UMIcountsPerCondition-1.png)<!-- -->

## Compare total UMI counts per marker

Plot total UMI counts for each marker at the investigated dilution
factors (DF1 vs. DF4). To ease readability, we place dashed lines
between each concentration.

``` r
## Calculate "breaks" where concentration change.
lines <- length(marker.order)-cumsum(sapply(split(ann.markerConc$Marker,ann.markerConc$conc_µg_per_mL),length))+0.5
lines <- data.frame(breaks=lines[-length(lines)])

## Make a marker by concentration "heatmap"
p.markerByConc <- ggplot(ann.markerConc, aes(x=1, y=Marker, fill=conc_µg_per_mL)) + 
  geom_tile(col=alpha(col="black",alpha=0.2)) + 
  geom_hline(data=lines,aes(yintercept=breaks), linetype="dashed", alpha=0.5) + 
  scale_fill_viridis_c(trans="log2") + 
  labs(fill="µg/mL") + 
  theme_get() + 
  theme(axis.ticks.x=element_blank(), axis.title = element_blank(), axis.text.x=element_blank(), panel.grid=element_blank(), legend.position="right", plot.margin=unit(c(0.1,0.1,0.1,0.1),"mm")) + scale_x_continuous(expand=c(0,0))
  
## Make UMI counts per Marker plot
p.UMIcountsPerMarker <- ggplot(ADT.matrix.agg, aes(x=marker.byConc,y=log2(sum))) + 
  geom_line(aes(group=marker), size=1.2, color="#666666") + 
  geom_point(aes(group=dilution, fill=dilution), pch=21, size=0.7) + 
  geom_vline(data=lines,aes(xintercept=breaks), linetype="dashed", alpha=0.5) + 
  scale_fill_manual(values=color.dilution) + 
  facet_wrap(~tissue) + 
  scale_y_continuous(breaks=c(9:17)) + 
  ylab("log2(UMI sum)") + 
  guides(fill=guide_legend(override.aes=list(size=1.5), reverse=TRUE)) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position="bottom", legend.justification="left", legend.title.align=0, legend.key.width=unit(0.2,"cm"), legend.title=element_blank()) + 
  coord_flip()

## Combine plot with markerByConc annotation heatmap
plotUMIcountsPerMarker <- p.markerByConc + guides(fill=F) + p.UMIcountsPerMarker + guides(fill=F) + plot_spacer() + guide_area() + plot_layout(ncol=4, widths=c(1,30,0.1), guides='collect')

plotUMIcountsPerMarker
```

![](Antibody-titration_files/figure-gfm/plotUMIcountsPerMarker-1.png)<!-- -->

## Compare change in UMI/cell within expressing cluster

Using a specific percentile may be prone to outliers in small clusters
(i.e. the 90th percentile of a cluster of 30 will be the \#3 higest cell
making it prone to outliers). We thus set a threshold of the value to
only be the 90th percentile if cluster contains more than 100 cells. For
smaller clusters, the median is used. Expressing cluster is identified
in the “preprocessing” vignette.

``` r
## Get the data
ADT.matrix <- data.frame(GetAssayData(object, assay="ADT.kallisto", slot="counts"))
ADT.matrix$marker <- rownames(ADT.matrix)
ADT.matrix$conc <- abpanel[ADT.matrix$marker,"conc_µg_per_mL"]
ADT.matrix <- ADT.matrix %>% pivot_longer(c(-marker,-conc))

## Get cell annotations
cell.annotation <- FetchData(object, vars=c("dilution", "tissue", "fineCluster"))

## Calculate marker statistics from each dilution within each cluster
ADT.matrix.agg <- ADT.matrix %>% group_by(dilution=cell.annotation[name,"dilution"], tissue=cell.annotation[name,"tissue"], fineCluster=cell.annotation[name,"fineCluster"], marker, conc) %>% summarise(sum=sum(value), median=quantile(value, probs=c(0.9)), nth=nth(value))

## Use data for the previously determined expressing cluster.
Cluster.max <- markerStats[,c("marker","tissue","fineCluster")]
Cluster.max$fineCluster <- factor(Cluster.max$fineCluster)
Cluster.max$tissue <- factor(Cluster.max$tissue, levels=c("PBMC","Lung"))

ADT.matrix.aggByClusterMax <- Cluster.max %>% left_join(ADT.matrix.agg)
ADT.matrix.aggByClusterMax$marker.byConc <- factor(ADT.matrix.aggByClusterMax$marker, levels=marker.order)

p.UMIinExpressingCells <- ggplot(ADT.matrix.aggByClusterMax, aes(x=marker.byConc, y=log2(nth))) + 
  geom_line(aes(group=marker), size=1.2, color="#666666") + 
  geom_point(aes(group=dilution, fill=dilution), pch=21, size=0.7) + 
  geom_vline(data=lines,aes(xintercept=breaks), linetype="dashed", alpha=0.5) + 
  geom_text(aes(label=paste0(fineCluster," ")), y=Inf, adj=1, size=1.5) + 
  facet_wrap(~tissue) + 
  scale_fill_manual(values=color.dilution) + 
  scale_y_continuous(breaks=c(0:11), labels=2^c(0:11), expand=c(0.05,0.5)) + 
  ylab("90th percentile UMI of expressing cluster") + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), legend.position="right", legend.justification="left", legend.title.align=0, legend.key.width=unit(0.2,"cm")) + 
  coord_flip()

## Combine plot with markerByConc annotation heatmap
UMIinExpressingCells <- p.markerByConc + theme(legend.position="none") + p.UMIinExpressingCells + theme(legend.position="none") + plot_spacer() + plot_layout(ncol=4, widths=c(1,30,0.1), guides='collect')

UMIinExpressingCells
```

![](Antibody-titration_files/figure-gfm/UMIinExpressingCells-1.png)<!-- -->

## ADT expression per cluster expression overview

Get data and calculate expression level and frequencies within positive
cells from each cluster

``` r
## Get the data
data.ADT.DSB <- GetAssayData(object, assay="ADT.kallisto", slot="data")
data.meta <- FetchData(object, vars=c("fineCluster","supercluster","dilution","tissue"))

## Determine expression level within positive fraction of each cluster
data.ADT.DSB.above <- as.data.frame(data.ADT.DSB) %>% 
                          mutate(marker=rownames(.)) %>% 
                          pivot_longer(-marker) %>% 
                          filter(data.meta[name,"dilution"]=="DF1" & 
                            value >= markerStats[paste(marker,data.meta[name,"tissue"],sep="_"),"DSB.cutoff"]) %>%
                          group_by(marker, 
                                   fineCluster=data.meta[name,"fineCluster"]) %>% 
                          summarize(count=length(value), 
                                    mean=mean(value), 
                                    sd=sd(value), 
                                    median=median(value), 
                                    nth=nth(value)) %>% 
                          arrange(fineCluster, marker)

## Convert from long into wide format using "mean" ADT signal within expression cells
data.ADT.DSB.above.matrix <- data.ADT.DSB.above %>% 
                                select(marker, fineCluster, mean) %>% 
                                pivot_wider(names_from=marker, values_from=mean, id_cols=fineCluster) %>%
                                replace(is.na(.), 0) %>% 
                                tibble::column_to_rownames("fineCluster")

## Get number of cells in each cluster
data.clusterCount <- as.numeric(table(data.meta$fineCluster))
names(data.clusterCount) <- levels(data.meta$fineCluster)

## Calculate positive frequency within each cluster
data.ADT.DSB.above.freq.matrix <- data.ADT.DSB.above %>% 
                                select(marker, fineCluster, count) %>% 
                                mutate(count=count/data.clusterCount[as.character(fineCluster)]) %>%
                                pivot_wider(names_from=marker, values_from=count, id_cols=fineCluster) %>%
                                replace(is.na(.), 0) %>% 
                                tibble::column_to_rownames("fineCluster")
```

Calculate annotation data

``` r
### Calculate annotation data
ann.conc <- log2(abpanel[colnames(data.ADT.DSB.above.matrix),"conc_µg_per_mL"])

## Determin
ann.fineCluster <- data.meta %>% 
  group_by(fineCluster, supercluster) %>% 
  summarize(PBMC.count = sum(tissue=="PBMC"), 
            Lung.count = sum(tissue=="Lung"), 
            count = length(tissue)) %>%
  mutate(tissue=ifelse((PBMC.count-Lung.count)>0,"PBMC","Lung")) %>% 
  filter(count > 1 & fineCluster %in% rownames(data.ADT.DSB.above.matrix)) %>%
  tibble::column_to_rownames("fineCluster")

ann.fineCluster <- ann.fineCluster %>% mutate(count.scaled=scales::rescale(count,to=c(0, 0.7), from=c(0,max(count))))
ann.fineCluster <- data.frame(ann.fineCluster)

## Calclu
ann.UMIuse <- ADT.matrix.agg.total %>% filter(dilution=="DF1")
ann.UMIuse.PBMC <- ann.UMIuse %>% filter(tissue=="PBMC") %>% as.data.frame()
rownames(ann.UMIuse.PBMC) <- ann.UMIuse.PBMC$marker
ann.UMIuse.Lung <- ann.UMIuse %>% filter(tissue=="Lung") %>% as.data.frame()
rownames(ann.UMIuse.Lung) <- ann.UMIuse.Lung$marker
```

Heatmap showing the fraction of positive cells (circle size) and their
expression level (color) of each antibody within each cluster. Including
additional annotation such as UMI counts within each tissue and antibody
concentration.

``` r
library(ComplexHeatmap)

## Set color schemes
col_fun = circlize::colorRamp2(c(0, 7.5, 15, 22.5, 30), c(viridis::inferno(5)))
col_fun.count = circlize::colorRamp2(c(0, 0.25, 0.5, 1), c("blue", "white", "yellow", "red"))
col_fun.conc = circlize::colorRamp2(c(-4, -2.25, -0.5, 1.25, 3), c(viridis::viridis(5)))

## Make row annotation object
row_ha = rowAnnotation("Tissue"=ann.fineCluster[,"tissue"], 
                       "Cell Type"=ann.fineCluster[,"supercluster"],
                       "Cluster size"=anno_simple(ann.fineCluster$count.scaled, 
                                                col=col_fun.count, 
                                                pch = 21, 
                                                pt_gp = gpar(col="black",
                                                             fill=col_fun.count(ann.fineCluster$count.scaled),
                                                             lwd=0.5), 
                                                pt_size = unit(ann.fineCluster$count.scaled+0.2, "npc"), 
                                                gp = gpar(fill=NA, 
                                                          lwd=0.5, 
                                                          col="black")),
                       col=list("Cell Type"=color.supercluster, "Tissue"=color.tissue),
                       annotation_legend_param = list(legend_direction = "horizontal"), 
                       show_legend=FALSE,
                       gp=gpar(color="black", lwd=0.5),
                       show_annotation_name=TRUE,
                       annotation_name_gp=gpar(fontsize=5, fontface="bold"),
                       annotation_width=unit(c(3,3,3),"mm"))


## Make column annotation object
column_ha = HeatmapAnnotation("UMI count (PBMC)"=ann.UMIuse.PBMC[colnames(data.ADT.DSB.above.matrix),"sum"]/10^3, 
                              "UMI count (Lung)"=ann.UMIuse.Lung[colnames(data.ADT.DSB.above.matrix),"sum"]/10^3,
                              "DF1 conc."=ann.conc, 
                              col=list("UMI count (PBMC)"=col_fun,
                                       "UMI count (Lung)"=col_fun,
                                       "DF1 conc."=col_fun.conc), 
                              gp=gpar(color="black", 
                                      lwd=0.5), 
                              annotation_legend_param = list(legend_direction = "horizontal"), 
                              show_legend=FALSE, 
                              show_annotation_name=TRUE,
                              annotation_name_gp=gpar(fontsize=5, fontface="bold"),
                              simple_anno_size=unit(3,"mm"), 
                              annotation_name_side="left",
                              which="column")

## Circle function for ComplexHeatmap
drawCircle <- function(j, i, x, y, width, height, fill){
  selected <- which(Cluster.max[,"fineCluster"] == rownames(data.ADT.DSB.above.matrix)[i] & Cluster.max[,"marker"]==colnames(data.ADT.DSB.above.matrix)[j])
  if(length(selected) > 0){
    grid.rect(x=x, y=y, width=width, height=height, gp=gpar(col=color.tissue[Cluster.max[selected,"tissue"]], 
                                                            fill=NA 
                                                           ))
  }
  radius <- data.ADT.DSB.above.freq.matrix[i, j]*1.8
  
  ## If positive population is less than a few percent percent, do not show it
  # What is the actual cutoff? this threhold is based on arbitrary values with arbitrary cutoff
  radius <- ifelse(radius > 0.02*1.8,radius+0.2, 0)
  grid.circle(x = x, y = y, r = radius*min(unit.c(width, height)),gp = gpar(fill = col_fun(data.ADT.DSB.above.matrix[i, j]), col=alpha("black",0.5), size=0.1))
}

## Make the heatmap
ht_list <- Heatmap(data.ADT.DSB.above.matrix, 
                   cell_fun=drawCircle, 
                   rect_gp = gpar(type = "none", lwd=0.5, col="black"), 
                   col = col_fun, 
                   top_annotation=column_ha, 
                   left_annotation=row_ha, 
                   border = FALSE,
                   name = "DSB", 
                   heatmap_legend_param = list(legend_direction = "horizontal"), 
                   show_heatmap_legend=FALSE,
                   heatmap_width = unit(7, "in"),
                   heatmap_height = unit(3, "in"),
                   row_names_gp = gpar(fontsize = 6, fontface="bold"),
                   row_names_side = "left",
                   column_names_gp = gpar(fontsize = 6),
                   row_dend_width = unit(3, "mm"),
                   column_dend_height = unit(3, "mm"),
                   column_dend_gp = gpar(lwd=0.5),
                   row_dend_gp = gpar(lwd=0.5)
                   )

## Draw the heatmap
draw_heatmap_now <- function(){
  draw(ht_list) + decorate_heatmap_body("DSB", {grid.grill(h=(0:21+0.5)/22, v=(0:52+0.5)/52, gp=gpar(alpha=0.05, col="black"))})
  popViewport()
}

draw_heatmap_now()
```

![](Antibody-titration_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Final plot

``` r
### First row of plots: tSNE overview
A <- p.tsne.cluster + 
  theme(axis.text.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.title.x=element_blank(), 
        axis.title.y=element_text(), 
        legend.position=c(0.5,0.01), 
        legend.justification=c(0.5,1), 
        legend.margin=margin(0,0,0,0), 
        legend.title = element_blank(), 
        legend.background=element_blank(), 
        plot.margin = unit(c(1,1,5,1),"mm")) + 
  guides(color=guide_legend(ncol=2, 
                            override.aes=list(shape=19, size=1.5), 
                            keyheight=unit(0.25, "cm"), 
                            keywidth=unit(0.15, "cm")))

B <- p.tsne.finecluster + 
  theme(axis.text.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.title.y=element_blank(), 
        axis.title.x=element_text())

C <- p.tsne.tissue + 
  theme(axis.text.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.title.y=element_blank(), 
        axis.title.x=element_blank(), 
        legend.position=c(0.5,0.01), 
        legend.justification=c(0.5,2), 
        legend.margin=margin(0,0,0,0), 
        legend.title = element_blank(), 
        legend.direction="horizontal", 
        legend.background=element_blank()) + 
  guides(color=guide_legend(override.aes=list(shape=19, size=1.5), 
                            keyheight=unit(0.25, "cm"), 
                            keywidth=unit(0.15, "cm")))

D <- p.UMIcountsPerCondition + 
  theme(legend.key.width=unit(0.3,"cm"), 
        legend.key.height=unit(0.4,"cm"), 
        legend.text=element_text(size=unit(5,"pt")))

## Combine A-D panels.
AD <- cowplot::plot_grid(A, B, C, D, 
                         labels=c("A", "B", "C", "D"), 
                         label_size=panel.label_size, 
                         nrow=1, 
                         vjust=panel.label_vjust, 
                         hjust=panel.label_hjust, 
                         align="h", 
                         axis="tb", 
                         rel_widths = c(2.2,2,2,2.2))

### Second row: Heatmap

## Make ComplexHeatmap output into a grid object for cowplot compatibility
E <- grid::grid.grabExpr(draw_heatmap_now(), wrap=TRUE)

## Make custom legend for the E (heatmap) panel. As other figures use "ggplot" legend style, we will remake the heatmap color bar in this style
E.legend <- cowplot::get_legend(ggplot(data.frame(values=c(0, 7.5, 15, 22.5, 30)), 
                                               aes(x=values,y=values,col=values)) + 
                                          scale_color_gradientn(colours=col_fun(c(0, 7.5, 15, 22.5, 30)), 
                                                                values=c(0, 7.5, 15, 22.5, 30)/30, 
                                                                breaks=c(0,30),
                                                                labels=c("Low","High")) + 
                                          geom_point() + 
                                          theme(legend.title=element_blank(), 
                                                legend.background=element_blank(), 
                                                legend.box.background=element_blank(), 
                                                legend.key=element_blank(), 
                                                legend.key.height=unit(2,"mm"), 
                                                legend.key.width=unit(3,"mm")))


### Third row: Titration "line" plots.

F1 <- p.markerByConc + theme(legend.position="none")
F2 <- p.UMIcountsPerMarker + theme(legend.position="none")
G1 <- p.UMIinExpressingCells + theme(legend.position="none")

F.legend <- cowplot::get_legend(p.UMIcountsPerMarker + 
                                   theme(legend.position="bottom", 
                                         legend.direction="horizontal", 
                                         legend.background=element_blank(), 
                                         legend.box.background=element_blank(), 
                                         legend.key=element_blank()))

## Combine F and G panels
FG <- cowplot::plot_grid(F1, F2, F1, G1, 
                         nrow=1, 
                         align="h", 
                         axis="tb", 
                         labels=c("F", "", "G", ""), 
                         label_size=panel.label_size, 
                         vjust=panel.label_vjust, 
                         hjust=panel.label_hjust, 
                         rel_widths=c(2,10,2,10))

## Combine A-D, E and F-G panels.
p.figure1 <- cowplot::ggdraw(cowplot::plot_grid(AD, E, FG, 
                                   ncol=1, 
                                   rel_heights=c(2.2,3.4,4.8), 
                                   align="v", 
                                   axis="l", 
                                   labels=c("", "E", ""), 
                                   label_size=panel.label_size, 
                                   vjust=panel.label_vjust, 
                                   hjust=panel.label_hjust)) + 
  ## Add custom legend for F and G panel
  cowplot::draw_plot(F.legend,0.08,0.017,0.2,0.00001) + 
  ## Add custom legend for E (heatmap) panel
  cowplot::draw_plot(E.legend,1,0.47,0.12,0.001,hjust=1, vjust=0)


## Make PNG of final figure
png(file=file.path(outdir,"Figure 1.png"), 
    width=figure.width.full, 
    height=9.5, 
    units = figure.unit, 
    res=figure.resolution, 
    antialias=figure.antialias)

  p.figure1

dev.off()
```

    ## png 
    ##   2

``` r
p.figure1
```

![](Antibody-titration_files/figure-gfm/Figure1-1.png)<!-- -->

## Titration examples

Reduction of antibody concentration had different inpacts on different
markers. We divided the markers into five categories and show examples
of each category below:

``` r
## Make helper function for plotting titration plots
titrationPlot <- function(marker, gate.PBMC=NULL, gate.Lung=NULL, y.axis=FALSE, show.gate=TRUE){
  curMarker.name <- marker
  
  ## Get antibody concentration for legends
  curMarker.DF1conc <- abpanel[curMarker.name, "conc_µg_per_mL"]
  if(show.gate==TRUE){
    ## Load gating percentages from manually set DSB thresholds
    gate <- markerStats[markerStats$marker == curMarker.name,c("tissue","pct")]
    colnames(gate) <- c("wrap","gate")
    gate$gate <- 1-(gate$gate/100)
    rownames(gate) <- gate$wrap
    ## Allow manual gating
    if(!is.null(gate.PBMC)) gate["PBMC","gate"] <- gate.PBMC
    if(!is.null(gate.Lung)) gate["Lung","gate"] <- gate.Lung
    gate$wrap <- factor(gate$wrap, levels=c("PBMC","Lung"))
  } else {
    gate <- NULL
  }

  p <- feature_rankplot_hist_custom(data=object, 
                                    marker=paste0("adt_",curMarker.name),      
                                    group="dilution",
                                    wrap="tissue",
                                    barcodeGroup="supercluster",
                                    conc=curMarker.DF1conc, 
                                    legend=FALSE, 
                                    yaxis.text=y.axis, 
                                    gates=gate,
                                    histogram.colors=color.dilution, 
                                    title=curMarker.name)
  
  return(p)
}

## Make titration plots for selected markers from each "category"
titrationPlots <- list()
# Category A
titrationPlots[["CD30"]] <- titrationPlot("CD30", y.axis=TRUE)
titrationPlots[["CD183"]] <- titrationPlot("CD183")
titrationPlots[["TCRgd"]] <- titrationPlot("TCRgd")
# Category B
titrationPlots[["HLA-DR"]] <- titrationPlot("HLA-DR", y.axis=TRUE, gate.PBMC=0.47, gate.Lung=0.295)
titrationPlots[["CD4"]] <- titrationPlot("CD4", gate.PBMC=0.28, gate.Lung=0.50)
titrationPlots[["CD19"]] <- titrationPlot("CD19", gate.PBMC=0.92, gate.Lung=0.75)
# Category C
titrationPlots[["CD3"]] <- titrationPlot("CD3", y.axis=TRUE, gate.PBMC=0.56, gate.Lung=0.46)
titrationPlots[["CD14"]] <- titrationPlot("CD14", gate.PBMC=0.68, gate.Lung=0.86)
titrationPlots[["CD39"]] <- titrationPlot("CD39", gate.PBMC=0.54, gate.Lung=0.33)
# Category D
titrationPlots[["HLA-ABC"]] <- titrationPlot("HLA-ABC", y.axis=TRUE)
titrationPlots[["CD44"]] <- titrationPlot("CD44")
titrationPlots[["CD45"]] <- titrationPlot("CD45")
# Category E
titrationPlots[["CD56"]] <- titrationPlot("CD56", y.axis=TRUE)
titrationPlots[["CD117"]] <- titrationPlot("CD117")
titrationPlots[["CD127"]] <- titrationPlot("CD127")

## Make a panel row for each category
A <- cowplot::plot_grid(plotlist=titrationPlots[1:3], nrow=1)
B <- cowplot::plot_grid(plotlist=titrationPlots[4:6], nrow=1)
C <- cowplot::plot_grid(plotlist=titrationPlots[7:9], nrow=1)
D <- cowplot::plot_grid(plotlist=titrationPlots[10:12], nrow=1)
E <- cowplot::plot_grid(plotlist=titrationPlots[13:15], nrow=1)

# a bit of a hack to make celltype legend
p.legend <- cowplot::get_legend(ggplot(data.frame(supercluster=object$supercluster), 
                                           aes(color=supercluster,x=1,y=1)) + 
  geom_point(shape=15, size=1.5) + 
  scale_color_manual(values=color.supercluster) + 
  theme(legend.title=element_blank(), 
        legend.margin=margin(0,0,0,0), 
        legend.key.size = unit(0.15,"cm"),
        legend.position = c(0.98,1.1), 
        legend.justification=c(1,1), 
        legend.direction="horizontal"))

## Combine final figure
p.figure2 <- cowplot::plot_grid(NULL, p.legend, A, NULL, B, NULL, C, NULL, D, NULL, E, 
                                labels=c("","","A","","B","","C","","D","","E"), 
                                ncol=1, 
                                rel_heights= c(1.3,1.5,20,1.5,20,1.5,20,1.5,20,1.5,20),
                                vjust=-0.5, 
                                hjust=panel.label_hjust, 
                                label_size=panel.label_size)


## Make PNG of final figure
png(file=file.path(outdir,"Figure 2.png"), 
    units=figure.unit, 
    res=figure.resolution, 
    width=figure.width.full, 
    height=10, 
    antialias=figure.antialias)

p.figure2

dev.off()
```

    ## png 
    ##   2

``` r
p.figure2
```

![](Antibody-titration_files/figure-gfm/titrationExamples-1.png)<!-- -->

## Individual titration plots

For supplementary information.

``` r
plots.columns = 3

markers <- abpanel[rownames(object[["ADT.kallisto"]]),]
markers <- markers[order(markers$Category, markers$Marker),]

## Divide titration plots into categories
plotCat <- list("A"=list(),"B"=list(),"C"=list(),"D"=list(),"E"=list())

## Make individual plots for each marker
for(i in 1:nrow(markers)){
  curMarker <- markers[i,]
  curMarker.name <- curMarker$Marker
  y.axis <- ifelse(length(plotCat[[curMarker$Category]]) %in% c(0,3,6,9,12),TRUE,FALSE)
  plotCat[[curMarker$Category]][[curMarker.name]] <- titrationPlot(curMarker.name, y.axis=y.axis)
}

## Make a supplementary figure for each category
for(i in seq_along(plotCat)){
  category <- names(plotCat)[i]
  plots <- plotCat[[category]]
  
  plots.num <- length(plots)
  plots.rows <- ceiling(plots.num/plots.columns)
  
  curPlots <- cowplot::plot_grid(plotlist=plots,ncol=plots.columns, align="h", axis="tb")
  curPlots.layout <- cowplot::plot_grid(NULL, p.legend, curPlots, vjust=-0.5, hjust=panel.label_hjust, label_size=panel.label_size, ncol=1, rel_heights= c(0.5, 1.3, 70/5*plots.rows))
  
  png(file=file.path(outdir,paste0("Supplementary Figure 2",category,".png")), 
      units=figure.unit, 
      res=figure.resolution, 
      width=figure.width.full, 
      height=(2*plots.rows),
      antialias=figure.antialias)

  print(curPlots.layout)
  
  dev.off()
  
  print(curPlots.layout)
}
```

![](Antibody-titration_files/figure-gfm/suppFig1-1.png)<!-- -->![](Antibody-titration_files/figure-gfm/suppFig1-2.png)<!-- -->![](Antibody-titration_files/figure-gfm/suppFig1-3.png)<!-- -->![](Antibody-titration_files/figure-gfm/suppFig1-4.png)<!-- -->![](Antibody-titration_files/figure-gfm/suppFig1-5.png)<!-- -->
