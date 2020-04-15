superclusters.levels <- c("T/NK cells","MO/MØ/DC","B/Plasma cells","Other")
color.supercluster <- RColorBrewer::brewer.pal(4,"Dark2")
names(color.supercluster) <- superclusters.levels

color.dilution <- c("DF1"="darkblue", "DF4"="lightblue")
color.tissue <- c("PBMC"="red", "Lung"="darkblue")
color.volume <- c("50µl"="darkgreen", "25µl"="darkorange")
color.number <- c("1000k"="black","200k"="grey")