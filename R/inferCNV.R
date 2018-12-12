

#' @details
#' The main functions you will need to use are CreateInfercnvObject() and run(infercnv_object).
#' For additional details on running the analysis step by step, please refer to the example vignette.
#' @aliases infercnv-package
"_PACKAGE"


#' The infercnv Class
#'
#' An infercnv object encapsulates the expression data and gene chromosome ordering information
#' that is leveraged by infercnv for data exploration.  The infercnv object is passed among the
#' infercnv data processing and plotting routines.
#'
#' Slots in the infercnv object include:
#'
#' @slot expr.data  <matrix>  the count or expression data matrix, manipulated throughout infercnv ops
#'
#' @slot count.data <matrix>  retains the original count data, but shrinks along with expr.data when genes are removed.
#' 
#' @slot gene_order  <data.frame> chromosomal gene order
#'
#' @slot reference_grouped_cell_indices <list>  mapping [['group_name']] to c(cell column indices) for reference (normal) cells
#'
#' @slot observation_grouped_cell_indices <list> mapping [['group_name']] to c(cell column indices) for observation (tumor) cells
#'
#' @export
#'

infercnv <- methods::setClass(
                         "infercnv",
                         slots = c(
                             expr.data = "ANY",
                             count.data = "ANY",
                             gene_order= "data.frame",
                             reference_grouped_cell_indices = "list",
                             observation_grouped_cell_indices = "list") )




#' @title CreateInfercnvObject
#'
#' @param raw_counts_matrix  the matrix of genes (rows) vs. cells (columns) containing the raw counts
#'                           If a filename is given, it'll be read via read.table()
#'                           otherwise, if matrix or Matrix, will use the data directly.
#' 
#' @param gene_order_file data file containing the positions of each gene along each chromosome in the genome.
#'
#' @param annotations_file a description of the cells, indicating the cell type classifications
#'
#' @param ref_group_names a vector containing the classifications of the reference (normal) cells to use for infering cnv
#'
#' @param delim delimiter used in the input files
#'
#' @description Creation of an infercnv object. This requires the following inputs:
#' A more detailed description of each input is provided below:
#'
#' The raw_counts_matrix:
#'
#'           MGH54_P16_F12 MGH53_P5_C12 MGH54_P12_C10 MGH54_P16_F02 MGH54_P11_C11  ...
#' DDX11L1     0.0000000     0.000000      0.000000      0.000000     0.0000000
#' WASH7P      0.0000000     2.231939      7.186235      5.284944     0.9650009
#' FAM138A     0.1709991     0.000000      0.000000      0.000000     0.0000000
#' OR4F5       0.0000000     0.000000      0.000000      0.000000     0.0000000
#' OR4F29      0.0000000     0.000000      0.000000      0.000000     0.0000000
#' ...
#'
#' The gene_order_file, contains chromosome, start, and stop position for each gene, tab-delimited:
#'
#'          chr  start   stop
#' DDX11L1 chr1  11869  14412
#' WASH7P  chr1  14363  29806
#' FAM138A chr1  34554  36081
#' OR4F5   chr1  69091  70008
#' OR4F29  chr1 367640 368634
#' OR4F16  chr1 621059 622053
#' ...
#' 
#' The annotations_file, containing the cell name and the cell type classification, tab-delimited.
#'
#'             V1                   V2
#' 1 MGH54_P2_C12 Microglia/Macrophage
#' 2 MGH36_P6_F03 Microglia/Macrophage
#' 3 MGH53_P4_H08 Microglia/Macrophage
#' 4 MGH53_P2_E09 Microglia/Macrophage
#' 5 MGH36_P5_E12 Oligodendrocytes (non-malignant)
#' 6 MGH54_P2_H07 Oligodendrocytes (non-malignant)
#' ...
#' 179  93_P9_H03 malignant
#' 180 93_P10_D04 malignant
#' 181  93_P8_G09 malignant
#' 182 93_P10_B10 malignant
#' 183  93_P9_C07 malignant
#' 184  93_P8_A12 malignant
#' ...
#'
#'
#' and the ref_group_names vector might look like so:  c("Microglia/Macrophage","Oligodendrocytes (non-malignant)")
#'
#' @return infercnv
#'
#' @export
#'

CreateInfercnvObject <- function(raw_counts_matrix, gene_order_file, annotations_file, ref_group_names, delim="\t") {

    # input expression data
    if (class(raw_counts_matrix) == "character") {
        raw.data <- read.table(raw_counts_matrix, sep=delim, header=TRUE, row.names=1, check.names=FALSE)    
        raw.data = as.matrix(raw.data)
    } else if (class(raw_counts_matrix) %in% c("dgCMatrix", "matrix", "data.frame")) {
        # use as is:
        raw.data <- raw_counts_matrix
    } else {
        stop("CreateInfercnvObject:: Error, raw_counts_matrix isn't recognized as a matrix, data.frame, or filename")
    }
    
    # get gene order info
    gene_order <- read.table(gene_order_file, header=FALSE, row.names=1, sep="\t")
    names(gene_order) <- c(C_CHR, C_START, C_STOP)

    # read annotations file
    input_classifications <- read.table(annotations_file, header=FALSE, row.names=1, sep=delim, stringsAsFactors=FALSE)

    # just in case the first line is a default header, remove it:
    if (rownames(input_classifications)[1] == "V1") {
        input_classifications = input_classifications[-1, , drop=F]
    }
    
    # make sure all reference samples are accounted for:
    if (! all( rownames(input_classifications) %in% colnames(raw.data)) ) {
        
        missing_cells <- rownames(input_classifications)[ ! ( rownames(input_classifications) %in% colnames(raw.data) ) ]

        error_message <- paste("Please make sure that all the annotated cell ",
                               "names match a sample in your data matrix. ",
                               "Attention to: ",
                               paste(missing_cells, collapse=","))
        stop(error_message)
    }

    # restrict expression data to the annotated cells.
    raw.data <- raw.data[,colnames(raw.data) %in% rownames(input_classifications)]
        
    # reorder cell classifications according to expression matrix column names
    input_classifications <- input_classifications[order(match(row.names(input_classifications), colnames(raw.data))), , drop=FALSE]

    
    # get indices for reference cells
    ref_group_cell_indices = list()
    for (name_group in ref_group_names) {
        cell_indices = which(input_classifications[,1] == name_group)
        if (length(cell_indices) == 0 ) {
            stop(sprintf("Error, not identifying cells with classification %s", name_group))
        }
        cell_names =  rownames(input_classifications)[cell_indices]
        ref_group_cell_indices[[ name_group ]] <- cell_indices
    }
    
    # rest of the cells are the 'observed' set.
    all_group_names <- unique(input_classifications[,1])
    obs_group_names <- setdiff(all_group_names, ref_group_names)
        
    # extract the genes indicated in the gene ordering file:
    order_ret <- .order_reduce(data=raw.data, genomic_position=gene_order)

    num_genes_removed = dim(raw.data)[1] - dim(order_ret$exp)[1]

    if (num_genes_removed > 0) {
        flog.info(paste("num genes removed taking into account provided gene ordering list: ",
                        num_genes_removed,
                        " = ",
                        num_genes_removed / dim(raw.data)[1] * 100,
                        "% removed.", sep=""))
    }
    
    raw.data <- order_ret$expr
    input_gene_order <- order_ret$order
    
    if(is.null(raw.data)) {
        error_message <- paste("None of the genes in the expression data",
                               "matched the genes in the reference genomic",
                               "position file. Analysis Stopped.")
        stop(error_message)
    }
    
    # define groupings according to the observation annotation names
    

    obs_group_cell_indices = list()
    for (name_group in obs_group_names) {
        cell_indices = which(input_classifications[,1] == name_group)
        cell_names = rownames(input_classifications)[cell_indices]
        obs_group_cell_indices[[ name_group ]] <- cell_indices
    }

        
    object <- new(
        Class = "infercnv",
        expr.data = raw.data, 
        count.data = raw.data,
        gene_order = input_gene_order,
        reference_grouped_cell_indices = ref_group_cell_indices,
        observation_grouped_cell_indices = obs_group_cell_indices)


    validate_infercnv_obj(object)
    
    return(object)
}


# Order the data and subset the data to data in the genomic position file.
#
# Args:
# @param data Data (expression) matrix where the row names should be in
#                 the row names of the genomic_position file.
# @param genomic_position Data frame read in from the genomic position file
#
# @return Returns a matrix of expression in the order of the
#            genomic_position file. NULL is returned if the genes in both
#            data parameters do not match.
#

.order_reduce <- function(data, genomic_position){
    flog.info(paste("::order_reduce:Start.", sep=""))
    ret_results <- list(expr=NULL, order=NULL, chr_order=NULL)
    if (is.null(data) || is.null(genomic_position)){
        return(ret_results)
    }

    # Drop pos_gen entries that are position 0
    remove_by_position <- -1 * which(genomic_position[2] + genomic_position[3] == 0)
    if (length(remove_by_position)) {
        flog.debug(paste("::process_data:order_reduce: removing genes specified by pos == 0, count: ",
                                length(remove_by_position), sep=""))

        genomic_position <- genomic_position[remove_by_position, , drop=FALSE]
    }

    # Reduce to genes in pos file

    flog.debug(paste("::process_data:order_reduce: gene identifers in expression matrix: ",
                            row.names(data), collapse="\n", sep=""))
    flog.debug(paste("::process_data:order_reduce: gene identifers in genomic position table: ",
                            row.names(data), collapse="\n", sep=""))



    keep_genes <- intersect(row.names(data), row.names(genomic_position))
        
    flog.debug(paste("::process_data:order_reduce: keep_genes size: ", length(keep_genes),
                     sep=""))
    
    # Keep genes found in position file
    if (length(keep_genes)) {
        ret_results$expr <- data[match(keep_genes, rownames(data)), , drop=FALSE]
        ret_results$order <- genomic_position[match(keep_genes, rownames(genomic_position)), , drop=FALSE]
    } else {
        flog.info(paste("::process_data:order_reduce:The position file ",
                        "and the expression file row (gene) names do not match."))
        return(list(expr=NULL, order=NULL, chr_order=NULL))
    }
    
    ## ensure expr and order match up!
    if (isTRUE(all.equal(rownames(ret_results$expr), rownames(ret_results$order)))) {
        flog.info(".order_reduce(): expr and order match.")
    }
    else {
        stop("Error, .order_reduce(): expr and order don't match! must debug")
    }
        
    # Set the chr to factor so the order can be arbitrarily set and sorted.
    chr_levels <- unique(genomic_position[[C_CHR]])
    ret_results$order[[C_CHR]] <- factor(ret_results$order[[C_CHR]],
                                   levels=chr_levels)
    
    # Sort genomic position file and expression file to genomic position file
    # Order genes by genomic region
    order_names <- row.names(ret_results$order)[with(ret_results$order, order(chr,start,stop))]
    
    ret_results$expr <- ret_results$expr[match(order_names, rownames(ret_results$expr)), , drop=FALSE] #na.omit is to rid teh duplicate gene entries (ie. Y_RNA, snoU13, ...) if they should exist.
    
    # This is the contig order, will be used in visualization.
    # Get the contig order in the same order as the genes.
    ret_results$order <- ret_results$order[match(order_names, rownames(ret_results$order)), , drop=FALSE]
    ret_results$chr_order <- ret_results$order[1]
    
    # Remove any gene without position information
    # Genes may be sorted correctly by not have position information
    # Here they are removed.
    flog.info(paste("::process_data:order_reduce:Reduction from positional ",
                           "data, new dimensions (r,c) = ",
                           paste(dim(data), collapse=","),
                           " Total=", sum(data),
                           " Min=", min(data),
                           " Max=", max(data),
                           ".", sep=""))
    flog.debug(paste("::process_data:order_reduce end."))
    return(ret_results)
}


#' @title remove_genes()
#'
#' @description infercnv obj accessor method to remove genes from the matrices
#'
#' @param infercnv_obj infercnv object
#' 
#' @param gene_indices_to_remove matrix indices for genes to remove
#'
#' @return infercnv_obj
#'
#' @export
#'

remove_genes <- function(infercnv_obj, gene_indices_to_remove) {

    infercnv_obj@expr.data <- infercnv_obj@expr.data[ -1 * gene_indices_to_remove, , drop=F]
    
    infercnv_obj@count.data <- infercnv_obj@count.data[ -1 * gene_indices_to_remove, , drop=F]

    infercnv_obj@gene_order <- infercnv_obj@gene_order[ -1 * gene_indices_to_remove, , drop=F] 

    validate_infercnv_obj(infercnv_obj)
    
    return(infercnv_obj)
}


#' @title validate_infercnv_obj()
#'
#' @description validate an infercnv_obj
#' ensures that order of genes in the @gene_order slot match up perfectly with the gene rows in the @expr.data matrix.
#' Otherwise, throws an error and stops execution.
#'
#' @param infercnv_obj infercnv_object
#'
#' @return none
#'

validate_infercnv_obj <- function(infercnv_obj) {

    flog.info("validating infercnv_obj")

    if (isTRUE(all.equal(rownames(infercnv_obj@expr.data), rownames(infercnv_obj@gene_order)))) {
            # all good.
                    return();

    }
        else {

        flog.error("hmm.... rownames(infercnv_obj@expr.data != rownames(infercnv_obj@gene_order))")
                broken.infercnv_obj = infercnv_obj
                        save('broken.infercnv_obj', file="broken.infercnv_obj")

    }


    genes = setdiff(rownames(infercnv_obj@expr.data), rownames(infercnv_obj@gene_order))
        if (length(genes) != 0) {
                flog.error(paste("The following genes are in infercnv_obj@expr.data and not @gene_order:", paste(genes, collapse=","),
                                         sep=" "))

    }

    genes = setdiff(rownames(infercnv_obj@gene_order), rownames(infercnv_obj@expr.data))
        if (length(genes) != 0) {
                flog.error(paste("The following genes are in @gene_order and not infercnv_obj@expr.data:", paste(genes, collapse=","),
                                         sep=" "))

    }

    stop("Problem detected w/ infercnv_obj")

}


get_cell_name_by_grouping <- function(infercnv_obj) {

    cell_name_groupings = list()
    
    groupings = c(infercnv_obj@reference_grouped_cell_indices, infercnv_obj@observation_grouped_cell_indices)

    for (group_name in names(groupings)) {

        cell_names = colnames(infercnv_obj@expr.data[, groupings[[ group_name ]] ] )

        cell_name_groupings[[ group_name ]] = cell_names

    }

    return(cell_name_groupings)
}



