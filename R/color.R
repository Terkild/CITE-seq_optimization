superclusters.levels <- c("T/NK cells","MO/MØ/DC","B/Plasma cells","Other")
color.supercluster <- RColorBrewer::brewer.pal(4,"Dark2")
names(color.supercluster) <- superclusters.levels

## Make color schemes for each condition that is "overlapping" so that 
## the same sample gets the same color even in other comparisons.
color.dilution <- c("DF1"="#e6194b", "DF4"="#0082c8")
color.volume <- c("50µl"="#0082c8", "25µl"="#f58231")
color.cellsAtStaining <- c("1000k"="#f58231","200k"="#911eb4")

color.tissue <- c("PBMC"="red", "Lung"="darkblue")