library("Seurat")

set.seed(114)
require("Seurat")
require("tidyverse")
theme_set(theme_bw() + 
            theme(
              axis.text.x=element_text(angle=45, hjust=1), 
              panel.grid.minor = element_blank(), 
              strip.background=element_blank(), 
              strip.text=element_text(face="bold", size=10)))

data10XADTDir <- "F:/Projects/ECCITE-seq/TotalSeqC_TitrationA/cellranger_A1/outs/raw_feature_bc_matrix"
dataCSCADTDir <- "F:/Projects/ECCITE-seq/TotalSeqC_TitrationA/cite-seq-count/A1_S5_d1_ADT/umi_count"

## Load helper functions (ggplot themes, biexp transformation etc.)
source("R/Utilities.R")
## Load predefined color schemes
source("R/color.R")

outdir <- "C:/Users/Terkild/OneDrive - Københavns Universitet/Koralovlab/ECCITE-seq/20200106 Titration 1"
datadir <- "D:/Data/DSBdata"

object <- readRDS(file.path(datadir,"H1_day0_demultilexed_singlets.RDS"))
object.neg <- readRDS(file.path(datadir,"neg_control_object.rds"))

object@raw.data[1:5,1:5]

res_mat <- object@raw.data
dim(res_mat)
bc_rank <- DropletUtils::barcodeRanks(res_mat, lower = 10)
knee_plot(bc_rank)

  res_mat <- object.neg@assay$CITE@raw.data
dim(res_mat)
bc_rank <- DropletUtils::barcodeRanks(res_mat, lower = 10)
knee_plot(bc_rank)
