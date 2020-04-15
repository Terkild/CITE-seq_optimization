library("Seurat")

datadir <- "F:/Data/DSBdata"

object <- readRDS(file.path(datadir,"H1_day0_demultilexed_singlets.RDS"))
object.neg <- readRDS(file.path(datadir,"neg_control_object.rds"))

object@assay$CITE
