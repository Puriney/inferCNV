#!/usr/bin/env Rscript
library(optparse)
library(ggplot2)
library(data.table)
library(gridExtra)
library(grid)
# Command line arguments
pargs <- optparse::OptionParser()
pargs <- optparse::add_option(pargs, c("--file_dir"),
                              type="character",
                              action="store",
                              dest="file_dir",
                              metavar="File_Directory",
                              help=paste("Path to the expression data files created by inferCNV.",
                                         "[Default %default][REQUIRED]"))
pargs <- optparse::add_option(pargs, c("--cluster_by_groups"),
                              type="logical",
                              action="store",
                              dest="cluster_by_groups",
                              default = "TRUE",
                              metavar="Cluster_By_Groups",
                              help=paste("Option to cluster the observed cells by unique type.",
                                         "[Default %default]"))
pargs <- optparse::add_option(pargs, c("--one_reference"),
                              type="logical",
                              action="store",
                              dest="one_reference",
                              default = "FALSE",
                              metavar="One_Reference",
                              help=paste("Option to treat all references as one.",
                                         "[Default %default]"))
pargs <- optparse::add_option(pargs, c("--genes"),
                              type="character",
                              action="store",
                              dest="genes",
                              metavar="Genes",
                              default = NULL,
                              help=paste("Desired Gene to plot.",
                                         "[Default %default]"))
pargs <- optparse::add_option(pargs, c("--contig"),
                              type="character",
                              action="store",
                              dest="contig",
                              default = NULL,
                              metavar="CHOOSE_CHROMOSOME",
                              help=paste("Chromosome in which the median density will be plotted.",
                                         "[Default %default]"))
args_parsed <- optparse::parse_args(pargs)




#' @title Plot Gene(s) or Chromosome Distributions Along InferCNV Steps, or Chromosome wide Statistical Plots. 
#' @description Takes in expression data created by inferCNV. If provided genes or a chromosome label,
#'              will plot their distribution along each step in inferCNV. If no gene or chromosome is 
#'              provided, will plot several statistical plots along all chromosomes. These statistics 
#'              include; expression a single cell line chosen from each grouping, Mean, Median Variance,
#'              MAD in expression, KS-values, t-tests, and linear regressions.
#'              Plots are output as a PDF file. 
#' @param files_dir (string) Path to where the expression data files produced by inferCNV are located.
#' @param cluster_by_groups (logical) Option to cluster the observed cells by unique type. Used for statistical plots.
#' @param one_reference (logical) Option to have references combined into one if TRUE or separated into different references if FALSE. Default value is FALSE.
#' @param genes (string or vector) Gene(s) of interest. Used for gene/contig plots.
#' @param contig (string) The label for the chromosome that will be plotted. Used for gene/contig plots.
#' 
#' @return 
#' PDF File with density plots showing gene expression levels along each step in inferCNV. 
#' @export
# Requires:
#   gridExtra, grid, data.table, ggplot2






pull_legend<-function(plot){
    plot_params <- ggplot_gtable(ggplot_build(plot))
    # get which grob in list is the legend called "guide=box" 
    legend_idx <- which(sapply(plot_params$grobs, function(x) x$name) == "guide-box")
    legend <- plot_params$grobs[[legend_idx]]
    return(legend)
}
# ----------------------Select Gene & Add Cell Types-----------------------------------------------------------------------
configure_stat_data <- function(infercnv_obj,
                           cluster_by_groups,
                           one_ref
) 
{
    # get the cell ID names in order, (reference, observed)
    row_order <- c(colnames(infercnv_obj@expr.data))
    
    cell_type <- replace(row_order, 1:length(row_order) %in% unlist(infercnv_obj@observation_grouped_cell_indices), paste("Observed"))
    
    ## Label the references based on index locations 
    if (one_ref == FALSE || length(infercnv_obj@reference_grouped_cell_indices) > 1) {
        for(i in 1:length(infercnv_obj@reference_grouped_cell_indices)){ 
            cell_type <- replace(cell_type, infercnv_obj@reference_grouped_cell_indices[[i]], paste("Reference",toString(i),sep = "")) 
        }
    } else {
        for(i in 1:length(infercnv_obj@reference_grouped_cell_indices)){ 
            cell_type <- replace(cell_type, 1:length(cell_type) %in% unlist(infercnv_obj@reference_grouped_cell_indices), paste("Reference"))
        }
    }
    plotting_data <- t(infercnv_obj@expr.data)
    # add the labled cell annotation information
    final_data <- data.table(plotting_data, cell_type = cell_type, check.names = FALSE, stringsAsFactors = FALSE)
    return(final_data)
}

configure_gene_data <- function(infercnv_obj,
                                gene,
                                cluster_by_groups,
                                one_reference,
                                contig){
    # Function to identify references and observed cell lines within expresion dataframe 
    #		Input: 
    #				infercnv_obj: (S4) infercnv object 
    #				gene: (string) Gene of interest.
    #               one_reference: (boolean) If True, will combine references into one 
    #       Output:
    #               Data frame with expression data and annotation information
    
    # get the cell ID names in order, (reference, observed)
    row_order <- c(colnames(infercnv_obj@expr.data))
    cell_type <- replace(row_order, 1:length(row_order) %in% unlist(infercnv_obj@observation_grouped_cell_indices), paste("Observed"))
    
    ## Label the references based on index locations 
    if (one_reference == FALSE || length(infercnv_obj@reference_grouped_cell_indices) > 1) {
        for(i in 1:length(infercnv_obj@reference_grouped_cell_indices)){ 
            cell_type <- replace(cell_type, infercnv_obj@reference_grouped_cell_indices[[i]], paste("Reference",toString(i),sep = "")) 
        }
    } else {
        for(i in 1:length(infercnv_obj@reference_grouped_cell_indices)){ 
            cell_type <- replace(cell_type, 1:length(cell_type) %in% unlist(infercnv_obj@reference_grouped_cell_indices), paste("Reference"))
        }
    }
    
    #filter the expression set to only include the gene of interest 
    if (!missing(gene)) {
        new_plot_data <- infercnv_obj@expr.data[gene,] 
        # add the labled cell annotation information
        final_data <- data.frame(new_plot_data, cell_type,
                                 check.rows = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
        colnames(final_data) <- c("exp","cell_type")
    } 
    if (!missing(contig)){
        # check if contig in data
        if (contig %in% infercnv_obj@gene_order[[1]]){
            final_data <- data.frame(t(infercnv_obj@expr.data), cell_type,
                                     check.rows = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
        } else {
            error_message <- paste("Contig selected is not found in the data.\n",
                                   "Check to make sure the correct contig label was selected.")
            stop(error_message)
        }
    }
    
    # check if all reference cells are in cell type 
    ref_cells <- which(cell_type != "Observed")
    reference_ids <- colnames(infercnv_obj@expr.data)[unlist(infercnv_obj@reference_grouped_cell_indices)]
    if (!(all(reference_ids %in% colnames(infercnv_obj@expr.data)[ref_cells]))){
        missing_refs <- reference_idx[which(!(reference_idx %in% names(cell_type)))]
        error_message <- paste("Error: Not all references are accounted for.",
                               "Make sure the reference names match the names in the data.\n",
                               "Check the following reference cell lines: ", 
                               paste(missing_refs, collapse = ","))
        stop(error_message)
    }
    # check if all observed cells are in cell type 
    observed_ids <-colnames(infercnv_obj@expr.data)[unlist(infercnv_obj@observation_grouped_cell_indices)]
    obs_cells <- which(cell_type == "Observed")
    if (!(all(observed_ids %in% colnames(infercnv_obj@expr.data)[obs_cells]))){
        missing_refs <- reference_idx[which(!(reference_idx %in% names(cell_type)))]
        error_message <- paste("Error: Not all observed cell lines are accounted for.",
                               "Make sure the cell id names match the names in the data.\n",
                               "Check the following oberved cell lines: ", 
                               paste(missing_refs, collapse = ","))
        stop(error_message)
    }
    return(final_data)
}

concat <- function(files){ 
    # Function to Create the titles for each graph:
    #   Gets the names of the files and capitalizes first letter of each word 
    #       input: the paths to file 
    #       output: vector of graph titles 
    sub_titles <- sapply(files, function(x){
        unlist(strsplit(basename(x),split="\\."))[1]
    })
    temp <- strsplit(sub_titles, "_")
    sapply(temp, function(x){
        tmp <- paste0(toupper(substring(x[-1], 1,1)), substring(x[-1], 2), collapse=" ")
        paste(x[1], tmp)
    })
}


# Function to generate the plots for the data:
StatisticalPlots <- function(files_dir,
                        cluster_by_groups = TRUE,
                        one_reference = FALSE){
    # ----------------------Load The InferCNV Object----------------------------------------------------------------------------------------
    env <- list.files(files_dir, pattern="*.infercnv_obj", full.names=TRUE)
    if (all(!(file.exists(env)))){
        path_error <- c(env)[which(!file.exists(env))]
        error_message <- paste("Error in argument pathway.",
                               "Please make sure the enviorments created by inferCNV are in the following pathway: ")
        stop(error_message, paste(path_error, collapse = ", "))
    } else {
        infercnv_obj<- get(load(file = env[length(env)]))
    }
    # ----------------------Load Data----------------------------------------------------------------------------------------
    dat <- configure_stat_data(infercnv_obj = infercnv_obj,
                          cluster_by_groups = cluster_by_groups,
                          one_ref = one_reference)
    # Create ouput directory
    dir.create(file.path(files_dir, "Plots"), showWarnings = FALSE) # will create directory if it doesnt exist 
    out_directory <- paste(files_dir, "Plots", sep = .Platform$file.sep)
    
    # ----------------------Create the plots for each gene of interest -----------------------------------------------------------------------
    #### get chromosome seperator locations
    ## get the correct order of the chromosomes
    ordering <- unique(infercnv_obj@gene_order[['chr']]) 
    ordered_locations <- table(infercnv_obj@gene_order[['chr']],useNA = "no")[ordering] 
    cumulativeCounts <- cumsum(c(1,ordered_locations[which(!is.na(names(ordered_locations)))]))
    chrIdx <- head(cumulativeCounts,n=-1L) # exclude the last number 
    geneOrder <- colnames(t(infercnv_obj@expr.data))
    # create PDF 
    pdf(paste(out_directory,"ExploratoryPlots.pdf", sep =.Platform$file.sep ), onefile = TRUE, height = 6, width = 9)
    
    # ----------------------plot 2 differnet cells  -----------------------------------------------------------------------
    message("Plotting a single reference cell against a observed cell.")
    #get index for one cell from each cell type
    single_idx <- sapply(unique(dat$cell_type), function(x){which(dat$cell_type==x)[1]})
    singleData <- dat[single_idx,]
    singleLongData <- melt(singleData, id = "cell_type")
    colnames(singleLongData) <- c("cell_type", "Gene", "Expression")
    singlePlot <- ggplot(singleLongData, aes(x = Gene, y = Expression, group = cell_type)) +
        geom_line(aes(color = cell_type)) +
        labs(title="Expression Of One Reference And One Observed Cell") + 
        theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 10),
              panel.background = element_blank(),
              panel.grid.major.x = element_line(colour = "grey"),
              legend.position = "bottom",
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
        scale_x_discrete(limits = geneOrder,
                         breaks = geneOrder[chrIdx],
                         labels = names(cumulativeCounts[-1]))
    print(singlePlot)
    
    
    # ----------------------PLOT MEAN VALUES-----------------------------------------------------------------------
    message("Plotting mean expression values.")
    geneMeans <- dat[,lapply(.SD,mean), by = cell_type]
    meanData <- melt(geneMeans, id = "cell_type")
    colnames(meanData) <- c("cell_type", "Gene", "Expression")
    meanPlot <- ggplot(meanData, aes(x = Gene, y = Expression, group = cell_type)) +
        geom_line(aes(color = cell_type)) +
        labs(title="Mean Expression") + 
        theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
              panel.background = element_blank(),
              panel.grid.major.x = element_line(colour = "grey"),
              legend.position = "bottom",
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
        scale_x_discrete(limits = geneOrder,
                         breaks = geneOrder[chrIdx],
                         labels = names(cumulativeCounts[-1]))
    print(meanPlot)
    
    # ----------------------MEDIAN VALUES -----------------------------------------------------------------------
    message("Plotting median expression values.")
    geneMedians <- dat[,lapply(.SD,median), by = cell_type]
    medianData <- melt(geneMedians, id = "cell_type")
    colnames(medianData) <- c("cell_type", "Gene", "Expression")
    medianPlot <- ggplot(medianData, aes(x = Gene, y = Expression, group = cell_type)) +
        geom_line(aes(color = cell_type)) +
        labs(title="Median Expression") +
        theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
              panel.background = element_blank(),
              panel.grid.major.x = element_line(colour = "grey"),
              legend.position = "bottom",
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
        scale_x_discrete(limits = geneOrder,
                         breaks = geneOrder[chrIdx],
                         labels = names(cumulativeCounts[-1]))
    print(medianPlot)
    
    # ----------------------VAR VALUES -----------------------------------------------------------------------
    message("Plotting Variance in expression within genes.")
    geneVar <- dat[,lapply(.SD,var), by = cell_type]
    varData <- melt(geneVar, id = "cell_type")
    colnames(varData) <- c("cell_type", "Gene", "Variance")
    varPlot <- ggplot(varData, aes(x = Gene, y = Variance, group = cell_type)) +
        geom_line(aes(color = cell_type)) +
        labs(title="Variance In Expression Within Genes") +
        theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
              panel.background = element_blank(),
              panel.grid.major.x = element_line(colour = "grey"),
              legend.position = "bottom",
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
        scale_x_discrete(limits = geneOrder,
                         breaks = geneOrder[chrIdx],
                         labels = names(cumulativeCounts[-1]))
    print(varPlot)
    
    # ----------------------MAD VALUES -----------------------------------------------------------------------
    message("Plotting MAD of gene expression.")
    geneMAD <- dat[,lapply(.SD,mad), by = cell_type]
    MADData <- melt(geneMAD, id = "cell_type")
    colnames(MADData) <- c("cell_type", "Gene", "MAD")
    madPlot <- ggplot(MADData, aes(x = Gene, y = MAD, group = cell_type)) +
        geom_line(aes(color = cell_type)) +
        labs(title="MAD Expression") +
        theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
              panel.background = element_blank(),
              panel.grid.major.x = element_line(colour = "grey"),
              legend.position = "bottom",
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
        scale_x_discrete(limits = geneOrder,
                         breaks = geneOrder[chrIdx],
                         labels = names(cumulativeCounts[-1]))
    print(madPlot)
    
    
    # ----------------------KS.test between genes -----------------------------------------------------------------------
    message("Plotting KS and P-values.")
    # seperate the observed and reference data 
    obs_dat <-  dat[unlist(infercnv_obj@observation_grouped_cell_indices),]
    ref_dat <- dat[-unlist(infercnv_obj@observation_grouped_cell_indices),]
    obs_dat <- as.data.frame(obs_dat)[,1:(ncol(obs_dat)-1)]
    ref_dat <- as.data.frame(ref_dat)[,1:(ncol(ref_dat)-1)] 
    
    ksValues <- lapply(geneOrder,function(i,d1,d2){
        ks_values <- ks.test(d1[,i],d2[,i])
        c(KS = ks_values$statistic, p.value = (-10*log10(ks_values$p.value)))}, obs_dat, ref_dat)
    names(ksValues) <- geneOrder
    ksData <- t(data.frame(ksValues, check.names = FALSE) )
    ksData <- transform(ksData, KS.D = as.numeric(KS.D),
                        p.value = as.numeric(p.value))
    ksData$Gene <- row.names(ksData)
    #plot KS values
    ksPlot <- ggplot(data = ksData, aes(x = Gene, group = 1)) +
        geom_line(aes(y = KS.D), color = "#FF9999") +
        ylab(label="KS Value") +
        theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
              panel.background = element_blank(),
              panel.grid.major.x = element_line(colour = "grey"),
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
        scale_x_discrete(limits = geneOrder,
                         breaks = geneOrder[chrIdx],
                         labels = names(cumulativeCounts[-1]))
    #plot p-values
    pValPlot <- ggplot(data = ksData, aes(x = Gene, group = 1)) +
        geom_line(aes(y = p.value), color = "skyblue2") +
        ylab(label="Phred Score (-10log10(P-Value))") +
        theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
              panel.background = element_blank(),
              panel.grid.major.x = element_line(colour = "grey"),
              legend.position = "none",
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
        scale_x_discrete(limits = geneOrder,
                         breaks = geneOrder[chrIdx],
                         labels = names(cumulativeCounts[-1]))
    
    # create the pdf and plot it on a grid 
    KSPLOT <- gridExtra::arrangeGrob( ksPlot, pValPlot,
                                      ncol=1, # use apply and get to get variables for each step instead
                                      top    = grid::textGrob("KS Values With Associated P-Values",
                                                              gp = grid::gpar(fontsize=20,font=2)),
                                      bottom = grid::textGrob("Gene",
                                                              gp = grid::gpar(fontsize=15)))
    gridExtra::grid.arrange(KSPLOT, ncol = 1, heights = c(6))
    
    # ----------------------t.test between genes -----------------------------------------------------------------------
    message("Plotting t-test.")
    # seperate the observed and reference data 
    plotTtest <- function(dat){
        obs_dat <-  dat[unlist(infercnv_obj@observation_grouped_cell_indices),]
        ref_dat <- dat[-unlist(infercnv_obj@observation_grouped_cell_indices),]
        obs_dat <- as.data.frame(obs_dat)[,1:(ncol(obs_dat)-1)]
        ref_dat <- as.data.frame(ref_dat)[,1:(ncol(ref_dat)-1)] 
        
        ttestValues <- data.frame(t(sapply(geneOrder,function(i,d1,d2){
            values <- (t.test(d1[,i],d2[,i]))
            ttest <- c(t.score = values$statistic, # t value the calculated difference represented in units of standard error (sigma)
                       t.pvalue = values$p.value,
                       confidence = values$conf.int) 
        }, obs_dat, ref_dat)))
        ttestValues$Gene <- row.names(ttestValues)
        ttestValues$t.fdr <- p.adjust(ttestValues$t.pvalue, method = "BH")
        ttestValues$t.pvalue <- -10*log10(ttestValues$t.pvalue)
        ttestValues$t.fdr <- -10*log10(ttestValues$t.fdr)
        
        ttest_score <- ggplot(data = ttestValues, aes(x = Gene)) +
            geom_line(aes(y = t.score.t, color = "T Value", group = 1)) +
            #geom_line(aes(y = confidence1, color = "blue", group = 1)) +
            labs(title="T-Value") +
            theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
                  panel.background = element_blank(),
                  panel.grid.major.x = element_line(colour = "grey"),
                  legend.position = "bottom",
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
            scale_x_discrete(limits = geneOrder,
                             breaks = geneOrder[chrIdx],
                             labels = names(cumulativeCounts[-1]))
        ttest_pval <- ggplot(data = ttestValues, aes(x = Gene)) +
            geom_line(aes(y = t.pvalue, color = "P-Value", group = 1)) +
            geom_line(aes(y = t.fdr, color = "Adjust p-Value", group = 1)) +
            labs(title="P-Value And Adjusted P-Value") +
            theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
                  panel.background = element_blank(),
                  panel.grid.major.x = element_line(colour = "grey"),
                  legend.position = "bottom",
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
            scale_x_discrete(limits = geneOrder,
                             breaks = geneOrder[chrIdx],
                             labels = names(cumulativeCounts[-1]))
        TTESTPLOT <- gridExtra::arrangeGrob( ttest_score, ttest_pval,
                                             ncol=1, # use apply and get to get variables for each step instead
                                             top    = grid::textGrob("t-Test Values With Associated P-Values",
                                                                     gp = grid::gpar(fontsize=20,font=2)))
        gridExtra::grid.arrange(TTESTPLOT, ncol = 1, heights = c(6))
        
    }
    my.t.test <- function(...) {
        temp<-try(plotTtest(...), silent=TRUE)
        if (is(temp, "try-error")) return(NA) else return(temp)
    }
    #print(plotTtest(dat))
    test <- my.t.test(dat)
    if (!is.na(test)){
        print(test)
    }
    
    # ----------------------linear regression between genes-----------------------------------------------------------------------
    linearTest <- function(dat, 
                           one_reference=one_reference){
        plotting <- function(dat, 
                             one_reference){
            # run linear models
            test <- lapply(1:(ncol(dat)-1), function(i){
                summary(lm(unlist(dat[,i,with =FALSE])~ dat$cell_type))
            })
            coef <- data.frame(t(sapply(1:length(test), function(x){
                temp <- test[[x]]$coefficients
                estimates <- temp[,1]
            })))
            pval <- data.frame(t(sapply(1:length(test), function(x){
                temp <- test[[x]]$coefficients
                estimates <- temp[,4]
            })))
            
            coef$Gene <- colnames(dat)[-ncol(dat)]
            coefData <- melt(coef, id = "Gene")
            colnames(coefData) <- c("Gene", "cell_type", "Coefficient")
            coefPlot <- ggplot(data = coefData, aes(x = Gene, y = Coefficient, colour = cell_type, fill = cell_type, group = cell_type)) +
                geom_line(alpha = .3) +
                scale_x_discrete(limits = geneOrder,
                                 breaks = geneOrder[chrIdx],
                                 labels = names(cumulativeCounts[-1])) +
                theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
                      panel.background = element_blank(),
                      panel.grid.major.x = element_line(colour = "grey"),
                      legend.position = "none",
                      axis.title.x = element_blank(),
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))
            
            pval$Gene <- colnames(dat)[-ncol(dat)]
            pvalData <- melt(pval, id = "Gene")
            colnames(pvalData) <- c("Gene", "cell_type", "pVal")
            lmPvalPlot <- ggplot(data = pvalData, aes(x = Gene, y = -10*log10(pVal), colour = cell_type, fill = cell_type, group = cell_type)) +
                geom_line(alpha = .3) +
                scale_x_discrete(limits = geneOrder,
                                 breaks = geneOrder[chrIdx],
                                 labels = names(cumulativeCounts[-1]))+
                theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 15),
                      panel.background = element_blank(),
                      panel.grid.major.x = element_line(colour = "grey"),
                      legend.position = "none",
                      axis.title.x = element_blank(),
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))
            lmPlot_list <- list(coefPlot,lmPvalPlot)
            return(lmPlot_list)
        }
        lmPlots <- c()
        if (length(infercnv_obj@reference_grouped_cell_indices) > 1 & one_reference == FALSE){
            dat$cell_type <- factor(dat$cell_type)
            types <- length(unique(dat$cell_type))
            refs <- sapply(1:length(infercnv_obj@reference_grouped_cell_indices), function(i){paste("Reference",toString(i),sep = "")}) 
            names(refs) <- refs
            
            plotrefs <- lapply(refs, function(i){
                contr.treatment(types, base = which(sort(unique(dat$cell_type)) == i))})
            for (i in 1:length(plotrefs)){
                contrasts(dat$cell_type) <- plotrefs[[i]]
                tempPlot <- plotting(dat, one_reference = one_reference)
                lmPlots <- c(lmPlots, tempPlot)}
            LMPLOT <- gridExtra::arrangeGrob( 
                grobs = lmPlots,
                ncol=1, # use apply and get to get variables for each step instead
                top    = grid::textGrob("Linear Regression Model Estimates and P-Values",
                                        gp = grid::gpar(fontsize=15,font=2)),
                bottom = grid::textGrob("Gene",
                                        gp = grid::gpar(fontsize=15)))
            plot_with_legend <- lmPlots[[1]] + theme(legend.position = "bottom") +
                scale_colour_discrete(labels=c("Observed","Reference","Reference")) +
                guides(fill = guide_legend(title = "Cell Type"), colour = guide_legend(title = "Cell Type"))
            plot_legend <- pull_legend(plot_with_legend)
            print(gridExtra::grid.arrange(LMPLOT,
                                          plot_legend, ncol = 1, heights = c(7,1)))
        } else {
            # need to make the reference the intercept 
            dat$cell_type <- factor(dat$cell_type)
            types <- length(unique(dat$cell_type))
            contrasts(dat$cell_type) <- contr.treatment(types, base = which(sort(unique(dat$cell_type)) == "Reference"))
            lmPlots <- plotting(dat, one_reference = one_reference)
            LMPLOT <- gridExtra::arrangeGrob( 
                grobs = lmPlots,
                ncol=1, # use apply and get to get variables for each step instead
                top    = grid::textGrob("Linear Regression Model Estimates and P-Values",
                                        gp = grid::gpar(fontsize=15,font=2)),
                bottom = grid::textGrob("Gene",
                                        gp = grid::gpar(fontsize=15)))
            plot_with_legend <- lmPlots[[1]] + 
                theme(legend.position = "bottom") +
                guides(fill = guide_legend(title = "Cell Type"), colour = guide_legend(title = "Cell Type"))
            plot_legend <- pull_legend(plot_with_legend)
            print(gridExtra::grid.arrange(LMPLOT, plot_legend, ncol = 1, heights = c(10,1)))
        }
    }
    linearTest(dat, 
               one_reference = one_reference)
    dev.off()
}


# ----------------------Function to plot Gene(s) or chromosome----------------------------------------------------------------------------------------
GeneContigPlot <- function (files_dir,
                         genes = NULL,
                         one_reference = FALSE,
                         contig = NULL){
    
    # ----------------------Check Pathways----------------------------------------------------------------------------------------
    # Error handling 
    ## check if pathway exists
    obj_files <- list.files(files_dir, pattern="*infercnv_obj", full.names=TRUE)
    
    if (any(!c(file.exists(files_dir, obj_files)))){
        path_error <- c(files_dir,obj_files)[which(!file.exists(files_dir,obj_files))]
        error_message <- paste("Error in argument pathway.",
                               "Please make sure the following pathways are correct: ")
        stop(error_message, paste(path_error, collapse = ", "))
    }
    # ----------------------Load infercnv objects----------------------------------------------------------------------------------------
    if ( length(obj_files) == 0 ) {
        stop("No enviornmental files found")
    }
    for (i in obj_files) {
        load(file = i)
    }
    
    obj_id <- lapply(obj_files, function(x){
        load(x)})
    final_obj <- get(tail(obj_id,n=1)[[1]])
    
    ref_groups = names(final_obj@reference_grouped_cell_indices)
    reference_idx = unlist(final_obj@reference_grouped_cell_indices)
    
    # ----------------------Create output directory----------------------------------------------------------------------------------------   
    dir.create(file.path(files_dir, "Plots"), showWarnings = FALSE) # will create directory if it doesnt exist 
    out_directory <- paste(files_dir, "Plots", sep = .Platform$file.sep)
    
    # ----------------------Create graph title----------------------------------------------------------------------------------------
    titles <- concat(obj_files)
    
    if (!is.null(genes)){
        if (is.character(genes)) {
            genes <- unlist(strsplit(genes, split = ","))
        }
        # check to see if the gene is in the data 
        logging::loginfo(paste0("Checking the following genes exists in data: ", paste(genes, collapse = ", ")))
        missing_gene <- lapply(obj_id, function(x){
            if (!(all(genes %in% row.names(get(x)@expr.data)))){
                genes[which(!(genes %in% row.names(get(x)@expr.data)))]
            }
        })
        if (!is.null(unlist(missing_gene))) {
            stop(paste0("Error: The following gene(s) are not found in expression data: ", paste(unique(unlist(missing_gene)), collapse = ", ")))
        }
        
        
        
        # ----------------------Create the plots for each gene of interest -----------------------------------------------------------------------
        for (gene in genes){
            logging::loginfo(paste0("Plotting ",gene))
            # run configure_gene_data function on each file to gather the expression data 
            data_list <- list()
            for (q in 1:length(obj_id)){
                q<-q
                data_list[[q]] <- configure_gene_data(infercnv_obj = get(obj_id[[q]]), 
                                                 gene = gene,
                                                 one_reference = one_reference)
            }
            # check that all cell id's match 
            if(length(unique(lapply(data_list, function(x) row.names(x))))!=1){
                stop("Cell id's don't match between plotting data.")
            }
            
            plot_list <- list()
            for (i in 1:length(data_list)) {
                plot_list[[i]] <-
                    assign(paste0("plot",i), 
                           ggplot(na.omit(data_list[[i]]), aes(x = exp, fill = cell_type, colour = cell_type)) +
                               geom_density(alpha = 0.3, adjust = 3/5) +
                               labs(title=titles[i]) +
                               theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 9),
                                     panel.background = element_blank(),
                                     legend.position = "none",
                                     axis.title.x = element_blank(),
                                     axis.title.y = element_blank())+
                               scale_x_continuous(breaks = pretty( data_list[[i]]$exp, n = 5)))
            }
            # Generate output file name with the desired output directory 
            ## File name is the name of the gene 
            output_file_name <- paste0(gene, ".pdf")
            pdf(file.path(out_directory,output_file_name), onefile = FALSE, height = 9, width = 9)
            temp <- plot_list[[1]] + 
                theme(legend.position = "bottom") + 
                guides(fill = guide_legend(title = "Cell Type"), colour = guide_legend(title = "Cell Type"))
            plot_legend <- pull_legend(temp)
            
            arrange_plots <- gridExtra::arrangeGrob( grobs = plot_list,
                                                     ncol=3, # use apply and get to get variables for each step instead
                                                     top    = grid::textGrob(gene,
                                                                             gp = grid::gpar(fontsize=20,font=2)),
                                                     left   = grid::textGrob("Density",
                                                                             gp = grid::gpar(fontsize=15), rot=90),
                                                     bottom = grid::textGrob("Expression Levels",
                                                                             gp = grid::gpar(fontsize=15)))
            # create the pdf and plot it on a grid 
            gridExtra::grid.arrange(arrange_plots, plot_legend, 
                                    nrow = 2, ncol = 1, 
                                    heights = c(10,1))
            dev.off()
        }
    }
    
    if (!is.null(contig)){
        data_list <- list()
        for (q in 1:length(obj_id)){
            q<-q
            data_list[[q]] <- configure_gene_data(infercnv_obj = get(obj_id[[q]]),
                                             one_reference = one_reference,
                                             contig = contig)
        }
        
        plot_list <- list()
        for (i in 1:length(data_list)) {
            current_order <- get(obj_id[[i]])@gene_order
            current_genes <- row.names(current_order[which(current_order[1] == contig),])
            
            geneMedians <- as.data.table(data_list[[i]])[,c("cell_type",current_genes), with = FALSE][,lapply(.SD,median), by = cell_type]
            df <- melt(geneMedians, id = "cell_type")
            colnames(df) <- c("cell_type", "Gene", "Median")
            plot_list[[i]] <- assign(paste0("plot",i), 
                                     ggplot(df, aes(x = Median, fill = cell_type, colour = cell_type)) + 
                                         geom_density(alpha = 0.3, adjust = 3/5) +
                                         labs(title=titles[i]) +
                                         theme(plot.title = element_text(hjust = 0.5, colour = "Black", size = 9),
                                               panel.background = element_blank(),
                                               legend.position = "none",
                                               axis.title.x = element_blank(),
                                               axis.title.y = element_blank())+
                                         scale_x_continuous(breaks = pretty(df$Median, n = 5)))
        }
        temp <- plot_list[[1]] + 
            theme(legend.position = "bottom") +
            guides(fill = guide_legend(title = "Cell Type"), colour = guide_legend(title = "Cell Type"))
        plot_legend <- pull_legend(temp)
        
        output_file_name <- paste0(contig,".pdf")
        pdf(file.path(out_directory,output_file_name), onefile = FALSE, height = 9, width = 9)
        # create the pdf and plot it on a grid 
        grid.arrange(arrangeGrob( grobs = plot_list,
                                  ncol=3, # use apply and get to get variables for each step instead
                                  top    = textGrob(contig,
                                                    gp = gpar(fontsize=20,font=2)),
                                  left   = textGrob("Density",
                                                    gp = gpar(fontsize=15), rot = 90),
                                  bottom = textGrob("Expression Levels",
                                                    gp = gpar(fontsize=15))),
                     plot_legend, nrow = 2, ncol = 1, heights = c(10,1))
    }
}

if (!is.null(args)) {
    if (is.null(args_parsed$genes) & is.null(args_parsed$contig)){
        StatisticalPlots(files_dir = args_parsed$file_dir,
                    cluster_by_groups = args_parsed$cluster_by_groups,
                    one_reference = args_parsed$one_reference)
    } else {
        GeneContigPlot(files_dir = args_parsed$file_dir, 
                    genes = args_parsed$genes,
                    one_reference = args_parsed$one_reference,
                    contig = args_parsed$contig)
    }
}
