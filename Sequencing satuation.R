set.seed(114)
require("Seurat", quietly=T)
require("tidyverse", quietly=T)

data.Seurat <- "data/5P-CITE-seq_Titration.rds"

dataCSCADTDirReads <- file.path(data.drive,data.project.dir,"cite-seq-count/A1_S5_d1_ADT_nocorrect/read_count")
dataCSCADTDir <- file.path(data.drive,data.project.dir,"cite-seq-count/A1_S5_d1_ADT_nocorrect/umi_count")

CSC.ADT.uncorrected <- Read10X(data.dir=dataCSCADTDir, gene.column=1)
CSC.ADT.uncorrected <- CSC.ADT.uncorrected[rownames(CSC.ADT.uncorrected) != "unmapped",]

CSC.ADT.uncorrected.reads <- Read10X(data.dir=dataCSCADTDirReads, gene.column=1)
CSC.ADT.uncorrected.reads <- CSC.ADT.uncorrected.reads[rownames(CSC.ADT.uncorrected.reads) != "unmapped",]

object <- readRDS(file=data.Seurat)

## Show number of cells from each sample
table(object$group)

object <- subset(object, subset=volume == "50µl")
object

UMI <- Matrix::colSums(CSC.ADT.uncorrected)[colnames(object)]
reads <- Matrix::colSums(CSC.ADT.uncorrected.reads)[colnames(object)]

df <- data.frame(barcode=colnames(object),sample=object$group,UMI=UMI,reads=reads)

df %>% group_by(sample) %>% summarize(UMI=sum(UMI), reads=sum(reads)) %>% mutate(saturation=(1-(UMI/reads)))
