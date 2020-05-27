feature_rankplot_hist_custom <- function(data,marker,group=NULL,barcodeGroup=NULL,wrap=NULL,conc=NULL,title=NULL,histogram.colors=c("red","blue"),barcode.colors=color.supercluster,...){

  data <- FetchData(data, vars=c(marker,barcodeGroup,group,wrap), slot = "counts")
  colnames(data)[1:3] <- c("value","barcodeGroup","group")

  color.manual <- histogram.colors

  if(is.null(wrap)){
    curWrap <- NULL
  } else {
    colnames(data)[4] <- "wrap"
  }

  if(group == "dilution"){
    if(!is.null(conc)){
      data$conc <- conc
      data$conc[data$group == "DF4"] <- conc/4
      data$conc <- factor(data$conc, levels=rev(sort(unique(data$conc))))
      levels(data$conc) <- sprintf("%2.2fug/mL",as.double(levels(data$conc)))
    } else {
      data$conc <- data$group
    }
    
    data$group <- data$conc
    curWrap <- data$wrap
    names(color.manual) <- levels(data$group)
  }

    p <- feature_rankplot_hist(data=data$value, 
                             group=data$group, 
                             wrap=curWrap, 
                             barcodeGroup=data$barcodeGroup, 
                             title=ifelse(!is.null(title),title,marker),
                             barcode.refGroups=levels(data$group)[1],
                             histogram.colors=color.manual,
                             barcode.colors=barcode.colors,
                             ...)
}