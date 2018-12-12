#!/usr/bin/env Rscript

#' @title Create Next Generation Clustered Heat Map (NG-CHM)
#' @description  Create highly interactive heat maps for single cell expression data using 
#' Next Generation Clustered Heat Map (NG-CHM). NG-CHM was developed and 
#' maintained by MD Anderson Department of Bioinformatics and Computational 
#' Biology in collaboration with In Silico Solutions. 
#'
#' @param infercnv_obj (S4) InferCNV S4 object holding expression data, gene location data, annotation information.
#' @param path_to_shaidyMapGen (string) Path to the java application ShaidyMapGen.jar
#' @param out_dir (string) Path to where the infercnv.ngchm output file should be saved to 
#' @param title (string) Title that will be used for the heatmap 
#' @param gene_symbol (string) Specify the label type that is given to the gene needed to create linkouts, default is NULL
#' @param x.center (integer) Center expression value for heatmap coloring.
#' @param x.range (integer) Values for minimum and maximum thresholds for heatmap coloring. 
#'
#' @return
#' 
#' Exports a NGCHM file named infercnv.ngchm and saves it to the output directory given to infercnv. 
# Requires:
#	NGCHM, ape, RcolorBrewer

Create_NGCHM <- function(infercnv_obj,
                         path_to_shaidyMapGen,
                         out_dir,
                         title = NULL, 
                         gene_symbol = NULL,
                         x.center = NA,
                         x.range = NA) {
    
    # ----------------------Check/create Pathways-----------------------------------------------------------------------------------------------
    ## check out_dir
    if (file.exists(out_dir)){
        file_path <- paste(out_dir, "infercnv.ngchm", sep = .Platform$file.sep)
    }else{
        dir.create(file.path(out_dir))
        paste("Creating the following Directory: ", out_dir)
    }
    
    #----------------------Initialize Next Generation Clustered Heat Map------------------------------------------------------------------
    # transpose the expression data so columns are the cell lines and rows are genes 
    plot_data <- t(infercnv_obj@expr.data)
    # create color map for the heat map and save it as a new data layer 
    # cut_value is the value that represents cuts on the heatmap
    cut_value <- -2147483648.0 
    # if specific center value is not given, set to 1 
    if (any(is.na(x.center))) {
        x.center <- 1
    }
    # if the range values are not given, will set appropriate values 
    if (! any(is.na(x.range))) {
        ## if the range values are provided, use defined values
        low_threshold <- x.range[1]
        high_threshold  <- x.range[2]
        if (low_threshold > x.center | high_threshold < x.center | low_threshold >= high_threshold) {
            x.center <- 0
            if (low_threshold > x.center | high_threshold < x.center | low_threshold >= high_threshold) {
                stop(paste("Error, problem with relative values of x.range: ", x.range, ", and x.center: ", x.center))
            }
        }
    } else {
        ## else, if not given, set the values 
        bounds <- get_average_bounds(infercnv_obj)
        low_threshold <- as.numeric(bounds[1])
        high_threshold <- as.numeric(bounds[2])
    }
    colMap <- NGCHM::chmNewColorMap(values        = c(cut_value, low_threshold, x.center, high_threshold),
                                    colors        = c("grey45","darkblue","white","darkred"),
                                    missing.color = "white", 
                                    type          = "linear") 
    layer <- NGCHM::chmNewDataLayer("DATA", as.matrix(plot_data), colMap, summarizationMethod = "average")
    # create the heat map object with the name "inferCNV"
    if (is.null(title)){
        title = "inferCNV"
    }
    hm <- NGCHM::chmNew(title, layer)
    # set the column (gene) order 
    hm@colOrder <- colnames(plot_data)
    hm@colOrderMethod <- "User"
    
    # add linkouts to each gene (column) for more information 
    if (!is.null(gene_symbol)) {
        hm <- NGCHM::chmAddAxisType(hm, 'col', gene_symbol)
    }
    
    ## set variables 
    ref_index = infercnv_obj@reference_grouped_cell_indices
    reference_idx = row.names(plot_data[unlist(ref_index),])
    ref_groups = names(ref_index)
    
    # ---------------------- Import Dendrogram & Order Rows -----------------------------------------------------------------------------------
    # IF Cluster By Group is set to TRUE:
    # Get the order of the rows (cell lines) from the dendrogram created by infer_cnv 
    
    # read the file containing the groupings created by infer_cnv
    row_groups_path <- paste(out_dir, "observation_groupings.txt", sep=.Platform$file.sep)
    row_groups <- read.table(row_groups_path, header = TRUE, check.names = FALSE) # genes are the row names 
    obs_order <- rev(row.names(row_groups)) # Reveerse names to correct order 
    row_order <- c(as.vector(reference_idx), obs_order) # put the reference cells above the observed cells 
    ## check for correct dimensions of new row order 
    if (length(row_order) != nrow(plot_data)) {
        stop("Error: After ordering the rows, row length does not match original dimensions of the data.
             \n Difference in row length: Original ", nrow(plot_data), ", After ordering ", length(row_order))
    }
    ## set the row order for the heatmap
    hm@rowOrder <- row_order
    hm@rowOrderMethod <- "User"
    
    # ----------------------Add Divisions Between References And Chromosomes ------------------------------------------------------------------
    # Column Separation: separation between the chromosomes 
    
    ## get the correct order of the chromosomes
    ordering <- unique(infercnv_obj@gene_order[['chr']]) 
    ## get gene locations in correct order, then find frequency of each chromosome
    ## add locations to each gene 
    location_data <- infercnv_obj@gene_order
    location_data$Gene <- row.names(infercnv_obj@gene_order)
    gene_order = colnames(plot_data)
    gene_locations_merge <- merge(data.frame(Gene = colnames(plot_data), stringsAsFactors = FALSE), location_data, by.x = "Gene")
    gene_locations <- gene_locations_merge[match(gene_order,gene_locations_merge$Gene),]
    # 
    # ## check if the number of genes has changed 
    if (nrow(gene_locations) != length(colnames(plot_data))){
        warning(paste0("Number of similar genes between expression data and locations:", nrow(gene_locations),
                       "\n Total number of genes in expression data: ", length(colnames(plot_data)),
                       "\n Check to make sure all the genes are in the location file and the gene names are the same between files."))
    }
    ## put in order 
    ordered_locations <- table(gene_locations[['chr']])[ordering] 
    cumulative_len <- cumsum(ordered_locations) #cumulative sum, separation locations
    sep_col_idx <- cumulative_len[-1 * length(cumulative_len)] # drop the last index because we do not want to add a break at the very end
    sep_col_idx <- rep(1,length(sep_col_idx)) + sep_col_idx # add one to each index, want to be to the right of the last gene in that chr
    ## colCutLocations: locations where the cuts will occur 
    ## colCutWidth: the width of the cuts 
    hm@colCutLocations <- as.integer(sep_col_idx)
    hm@colCutWidth <- as.integer(30)
    
    # Row separation: separation between reference samples and observed samples 
    hm@rowCutLocations <- as.integer(length(reference_idx)+1)
    # make the size of the separation proportional to the size of the heat map 
    ## use ncol because plot_data has cell ID's as the columns 
    row_sep <- ceiling(nrow(plot_data)*.01)
    hm@rowCutWidth <- as.integer(row_sep)
    
    #----------------------Create Covariate Bar----------------------------------------------------------------------------------------------------------------------------------------
    # Returns the color palette for contigs.
    get_group_color_palette <- function(){
        return(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3")))
    }
    # # check to make sure all cell lines are included 
    if (!(all(obs_order %in% row_order))) {
        missing_ids <- row_groups[which(!(obs_order %in% row_order))]
        error_message <- paste("Groupings of cell line ID's in observation_groupings.txt \n",
                               "do not match the ID's in the expression data.\n",
                               "Check the following cell line ID's: ",
                               paste(missing_ids, collapse = ","))
    }
    #---------------------COLUMN Covariate bar----------------------------------------------------------------------------------------------------------------------
    # COLUMN Covariate bar
    ## map the genes to their chromosome 
    ## gene_locations: created earlier, Genes and their locations 
    chr_labels <- as.vector(ordering)
    ## get the chromosomes 
    chr <- as.character(gene_locations$chr)
    ## get the gene ID's 
    names(chr) <- gene_locations$Gene
    
    chr_palette <- get_group_color_palette()(length(unique(location_data$chr)))
    names(chr_palette) <- unique(location_data$chr)
    
    ## create color mapping
    colMap_chr <- NGCHM::chmNewColorMap(values        = as.vector(chr_labels),
                                        colors        = chr_palette,
                                        missing.color = "white")
    chr_cov <- NGCHM::chmNewCovariate(fullname        = 'Chromosome', 
                                      values           = chr, 
                                      value.properties = colMap_chr,
                                      type             = "discrete")
    hm <- NGCHM::chmAddCovariateBar(hm, "column", chr_cov, 
                                    display   = "visible", 
                                    thickness = as.integer(20))
    
    #---------------------ROW Covariate bar----------------------------------------------------------------------------------------------------------------------
    
    # create covariate bar from dendrogram groups
    ## row_groups is taken from the dendrogram created by inferCNV
    ## create better column names 
    colnames(row_groups) <- c("Dendrogram.Group", "Dendrogram.Color", "Annotation.Group", "Annotation.Color")
    dendrogram_col <- as.character(unlist(row_groups["Dendrogram.Color"]))# group colors
    dendrogram_group <- as.character(unlist(row_groups["Dendrogram.Group"]))# group number
    dendrogram_unique_group <- unique(dendrogram_group)
    cells <- row.names(row_groups) # cell line ID's
    names(dendrogram_col) <- cells
    names(dendrogram_group) <- cells
    dendrogram_palette <- get_group_color_palette()(length(unique(dendrogram_col)))
    ## create color mapping
    colMap_dendrogram <- NGCHM::chmNewColorMap(values        = as.vector(dendrogram_unique_group),
                                               colors        = dendrogram_palette,
                                               missing.color = "white")
    dendrogram_cov <- NGCHM::chmNewCovariate(fullname         = 'Dendrogram',
                                             values           = dendrogram_group,
                                             value.properties = colMap_dendrogram,
                                             type             = "discrete")
    hm <- NGCHM::chmAddCovariateBar(hm, "row", dendrogram_cov,
                                    display   = "visible",
                                    thickness = as.integer(20))
    
    # Covariate to identify Reference and Observed data
    annotation_col <- as.character(unlist(row_groups["Annotation.Color"])) # group colors
    annotation_group <- as.character(unlist(row_groups["Annotation.Group"]))# group number
    names(annotation_group) <- cells
    names(annotation_col) <- cells
    annotation_unique_group <- unique(annotation_group)
    
    len <-lengths(ref_index)
    ref_bar_labels <- unlist(sapply(1:length(len), function(x){ rep(ref_groups[x],len[x]) }))
    names(ref_bar_labels) <- reference_idx
    
    # if you want the exact coloring as the original inferCNV plots 
    #annotation_palette <- c(get_group_color_palette()(length(ref_index)), get_group_color_palette()(length(annotation_unique_group)))
    
    # combine reference and observed labels 
    annotation_group <- c(ref_bar_labels,annotation_group)
    
    # change the observed group names in bar to group namnes 
    observed_data <- infercnv_obj@observation_grouped_cell_indices
    lapply(1:length(observed_data), function(x) { 
        tmp <- names(observed_data[x])
        annotation_group <<- replace(annotation_group, observed_data[[x]], tmp) } )
    unique_group <- unique(annotation_group)
    annotation_palette <- get_group_color_palette()(length(unique_group))
    
    # check if all reference cells are included 
    if (!(all(reference_idx %in% names(annotation_group)))){
        missing_refs <- reference_idx[which(!(reference_idx %in% names(annotation_group)))]
        error_message <- paste("Error: Not all references are accounted for.",
                               "Make sure the reference names match the names in the data.\n",
                               "Check the following reference cell lines: ", 
                               paste(missing_refs, collapse = ","))
        stop(error_message)
    }
    # check if all observed cells are included
    observed_idx <- row.names(plot_data[unlist(infercnv_obj@observation_grouped_cell_indices),])
    if (!(all(observed_idx %in% names(annotation_group)))){
        missing_obs <- reference_idx[which(!(observed_idx %in% names(annotation_group)))]
        error_message <- paste("Error: Not all observed cell lines are accounted for.",
                               "Make sure the reference names match the names in the data.\n",
                               "Check the following reference cell lines: ", 
                               paste(missing_obs, collapse = ","))
        stop(error_message)
    }
    
    ## create color mapping
    colMap_annotation <- NGCHM::chmNewColorMap(values        = as.vector(unique_group), 
                                               colors        = annotation_palette,
                                               missing.color = "white")
    annotation_cov <- NGCHM::chmNewCovariate(fullname         = 'Annotation', 
                                             values           = annotation_group, 
                                             value.properties = colMap_annotation,
                                             type             = "discrete")
    hm <- NGCHM::chmAddCovariateBar(hm, "row", annotation_cov, 
                                    display   = "visible", 
                                    thickness = as.integer(20))
    
    #---------------------------------------Export the heat map-----------------------------------------------------------------------------------------------------------------------
    ## adjust the size of the heat map 
    #hm@width <- as.integer(500)
    #hm@height <- as.integer(500)
    ## adjust label display size 
    #hm@rowDisplayLength <- as.integer(10) 
    futile.logger::flog.info(paste("Saving new NGCHM object"))
    NGCHM::chmExportToFile(hm, file_path, overwrite = TRUE, shaidyMapGen = path_to_shaidyMapGen)
    }

