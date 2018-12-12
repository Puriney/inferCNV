#!/usr/bin/env Rscript


options(error = function() { traceback(2); stop("Error encountered") } )

C_CHR <- "chr"
C_START <- "start"
C_STOP <- "stop"
C_HCLUST_METHODS <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
C_OUTPUT_FORMAT <- c("pdf", "png")

#' @importFrom grDevices col2rgb colorRampPalette dev.off pdf png rgb
#' @importFrom graphics abline axis boxplot hist image layout lines mtext par plot points rect text title legend
#' @importFrom stats as.dendrogram as.dist cutree density dist filter median order.dendrogram quantile reorder sd complete.cases cor t.test p.adjust predict rnorm runif smooth.spline var wilcox.test
#' @importFrom utils flush.console read.table write.table tail
#' @importFrom binhf ansc
#' @import futile.logger
#' @importFrom methods setClass new is
#' @importFrom gplots bluered
#' @importFrom ape write.tree as.phylo
#' @importFrom fastcluster hclust
#' @import RColorBrewer
#' @importFrom Matrix Matrix rowMeans colSums
#' @import coin
#' @importFrom dplyr %>% count

NULL

