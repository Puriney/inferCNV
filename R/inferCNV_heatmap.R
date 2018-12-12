#!/usr/bin/env Rscript


# Returns the color palette for contigs.
#
# Returns:
# Color Palette
get_group_color_palette <- function(){
    return(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3")))
}


#' @description Formats the data and sends it for plotting.
#'
#' @title Plot the matrix as a heatmap, with cells as rows and genes as columns, ordered according to chromosome
#'
#' @param infercnv_obj infercnv object
#' @param out_dir Directory in which to save pdf and other output.
#' @param title Plot title.
#' @param obs_title Title for the observations matrix.
#' @param ref_title Title for the reference matrix.
#' @param cluster_by_groups Whether to cluster observations by their annotations or not. Using this ignores k_obs_groups.
#' @param k_obs_groups Number of groups to break observation into.
#' @param contig_cex Contig text size. 
#' @param x.center Value on which to center expression.
#' @param x.range vector containing the extreme values in the heatmap (ie. c(-3,4) )
#' @param hclust_method Clustering method to use for hclust.
#' @param color_safe_pal Logical indication of using a color blindness safe palette.
#' @param output_filename Filename to save the figure to.
#' @param output_format format for heatmap image file (default: 'png'), options('png', 'pdf', NA)
#'                      If set to NA, will print graphics natively
#' @param ref_contig If given, will focus cluster on only genes in this contig.
#' @param write_expr_matrix Includes writing a matrix file containing the expression data that is plotted in the heatmap.
#' 
#' @return
#' No return, void.
#'
#' @export
#'


plot_cnv <- function(infercnv_obj,
                     out_dir=".",
                     title="inferCNV",
                     obs_title="Observations (Cells)",
                     ref_title="References (Cells)",
                     cluster_by_groups=TRUE,
                     k_obs_groups = 3,
                     contig_cex=1,
                     x.center=0,
                     x.range=NA,
                     hclust_method='ward.D',
                     color_safe_pal=TRUE,
                     output_filename="infercnv",
                     output_format="png", #pdf, png, NA
                     ref_contig = NULL,
                     write_expr_matrix=FALSE) {


    # arg validations
    if (! hclust_method %in% C_HCLUST_METHODS) {
        stop(sprintf("Error, hclust_method: %s is not supported", hclust_method))
    }
    if ( (! is.na(output_format) ) & (!  output_format %in% C_OUTPUT_FORMAT) )  {
        stop(sprintf("Error, output_format: %s is not supported", output_format) )
    }
    
    if(!file.exists(out_dir)){
        dir.create(out_dir)
    }
    
    plot_data = infercnv_obj@expr.data
    
    flog.info(paste("::plot_cnv:Start", sep=""))
    flog.info(paste("::plot_cnv:Current data dimensions (r,c)=",
                           paste(dim(plot_data), collapse=","),
                           " Total=", sum(plot_data, na.rm=T),
                           " Min=", min(plot_data, na.rm=T),
                           " Max=", max(plot_data, na.rm=T),
                           ".", sep=""))
    flog.info(paste("::plot_cnv:Depending on the size of the matrix",
                           " this may take a moment.",
                           sep=""))


    
    if (write_expr_matrix) {
        expr_dat_file <- paste(out_dir, paste("expr.", output_filename, ".dat", sep=""), sep="/")

        if (class(plot_data) %in% c("matrix", "data.frame")) {
            write.table(plot_data, file=expr_dat_file, quote=F, sep="\t")
        }
        
    }
    
    if (! any(is.na(x.range))) {


        if ( (length(x.range) == 1) & (x.range[1] == "auto") ) {

            # examine distribution of data that's off-center, since much of the center could
            # correspond to a mass of data that has been wiped out during noise reduction
            quantiles = quantile(plot_data[plot_data != x.center], c(0.05, 0.95))

            # determine max distance from the center.
            delta = max( abs( c(x.center - quantiles[1],  quantiles[2] - x.center) ) )
            low_threshold = x.center - delta
            high_threshold = x.center + delta
            
            flog.info(sprintf("plot_cnv(): auto thresholding at: (%f , %f)", low_threshold, high_threshold))
            
        } else {
        
            # use defined values
            low_threshold = x.range[1]
            high_threshold = x.range[2]
            
            if (low_threshold > x.center | high_threshold < x.center | low_threshold >= high_threshold) {
                stop(paste("Error, problem with relative values of x.range: ", x.range, ", and x.center: ", x.center))
            }
        }

        plot_data[plot_data < low_threshold] <- low_threshold
        plot_data[plot_data > high_threshold] <- high_threshold
        
        infercnv_obj@expr.data <- plot_data  #because used again below...
        
    }
    
    
    # Contigs
    contigs = infercnv_obj@gene_order[[C_CHR]]
    unique_contigs <- unique(contigs)
    n_contig <- length(unique_contigs)
    ct.colors <- get_group_color_palette()(n_contig)
    names(ct.colors) <- unique_contigs

    # Select color palette
    custom_pal <- color.palette(c("purple3", "white", "darkorange2"),
                                c(2, 2))
    if (color_safe_pal == FALSE){
        custom_pal <- color.palette(c("darkblue", "white", "darkred"),
                                    c(2, 2))
    }
    
    # Row separation based on reference
    ref_idx <- unlist(infercnv_obj@reference_grouped_cell_indices)
    ref_idx = ref_idx[order(ref_idx)]
        
    # Column seperation by contig and label axes with only one instance of name
    contig_tbl <- table(contigs)[unique_contigs]
    col_sep <- cumsum(contig_tbl)
    col_sep <- col_sep[-1 * length(col_sep)]   ## FIXME:  removing last entry?
    # These labels are axes labels, indicating contigs at the first column only
    # and leaving the rest blank.
    contig_labels <- c()
    contig_names <-c()
    for (contig_name in names(contig_tbl)){
        contig_labels <- c(contig_labels,
                           contig_name,
                           rep("", contig_tbl[contig_name] - 1))
        contig_names <- c(contig_names,rep(contig_name,contig_tbl[contig_name]))
    }

    # Calculate how many rows will be made for the number of columns in the grouping key
    grouping_key_coln <- c()
    obs_annotations_names <- names(infercnv_obj@observation_grouped_cell_indices)

    
    # obs_annotations_groups: integer vec named by cells, set to index according to category name vec above.
    obs_annotations_groups = rep(-1, length(colnames(infercnv_obj@expr.data))) # init
    names(obs_annotations_groups) = colnames(infercnv_obj@expr.data)
    obs_index_groupings = infercnv_obj@observation_grouped_cell_indices
    counter <- 1
    for (obs_index_group in obs_index_groupings) {
        obs_annotations_groups[ obs_index_group ] <- counter
        counter <- counter + 1
    }
    # restrict to just the obs indices

    obs_annotations_groups <- obs_annotations_groups[-ref_idx]
    
    grouping_key_coln[1] <- floor(123/(max(nchar(obs_annotations_names)) + 4))  ## 123 is the max width in number of characters, 4 is the space taken by the color box itself and the spacing around it
    if (grouping_key_coln[1] < 1) {
        grouping_key_coln[1] <- 1
    }

    name_ref_groups = names(infercnv_obj@reference_grouped_cell_indices)
    grouping_key_coln[2] <- floor(123/(max(nchar(name_ref_groups)) + 4))  ## 123 is the max width in number of characters, 4 is the space taken by the color box itself and the spacing around it
    if (grouping_key_coln[2] < 1) {
        grouping_key_coln[2] <- 1
    }
    
    grouping_key_rown <- c()
    grouping_key_rown[1] <- ceiling(length(obs_annotations_names)/grouping_key_coln[1])
    grouping_key_rown[2] <- ceiling(length(name_ref_groups)/grouping_key_coln[2])
    # Calculate how much bigger the output needs to be to accodomate for the grouping key
    grouping_key_height <- c((grouping_key_rown[2] + 2) * 0.175, (grouping_key_rown[1] + 3) * 0.175)

    # Rows observations, Columns CHR
    if (! is.na(output_format)) {
        if (output_format == "pdf") {
            pdf(paste(out_dir, paste(output_filename, ".pdf", sep=""), sep="/"),
                useDingbats=FALSE,
                width=10,
                height=(8.13 + sum(grouping_key_height)),
                paper="USr")
        } else if (output_format == "png") {
            png(paste(out_dir, paste(output_filename, ".png", sep=""), sep="/"),
                width=10,
                height=(8.13 + sum(grouping_key_height)),
                units="in",
                res=600)
        }
    }
    
    # Plot observations
    ## Make Observation Samples
    ## Remove observation col names, too many to plot
    ## Will try and keep the reference names
    ## They are more informative anyway
    obs_data <- infercnv_obj@expr.data
    if (!is.null(ref_idx)){
        obs_data <- plot_data[, -ref_idx, drop=FALSE]
        if (ncol(obs_data) == 1) {
            # hack for dealing with single entries
            plot_data <- cbind(obs_data, obs_data)
            names(obs_data) <- c("", names(obs_data)[1])
        }
    }
    
    obs_data_t <- t(obs_data)

    # Subsample the data to only the references and update the ref_group indexes to match their new indexes
    # ref_data_t <- plot_data[, ref_idx, drop=FALSE]
    ref_data_t <- NULL
    updated_ref_groups <- list()
    current_ref_count <- 1
    current_grp_idx <- 1
    plot_data <-infercnv_obj@expr.data
    ref_groups = infercnv_obj@reference_grouped_cell_indices
    for (ref_grp in ref_groups) {
        ref_data_t <- cbind(ref_data_t, plot_data[, ref_grp, drop=FALSE])
        updated_ref_groups[[current_grp_idx]] = seq(current_ref_count, current_ref_count + length(ref_grp) - 1)
        current_ref_count <- current_ref_count + length(ref_grp)
        current_grp_idx <- current_grp_idx + 1
    }
    ref_groups <- updated_ref_groups
    
    nb_breaks <- 16
    breaksList_t <-
        seq(min(min(obs_data_t, na.rm=TRUE), min(ref_data_t, na.rm=TRUE)),
        max(max(obs_data_t,na.rm=TRUE), max(ref_data_t, na.rm=TRUE)),
        length.out=nb_breaks)


    # Create file base for plotting output
    force_layout <- .plot_observations_layout(grouping_key_height=grouping_key_height)
    .plot_cnv_observations(obs_data=obs_data_t,
                          file_base_name=out_dir,
                          output_filename_prefix=output_filename,
                          cluster_contig=ref_contig,
                          contig_colors=ct.colors[contigs],
                          contig_labels=contig_labels,
                          contig_names=contig_names,
                          col_pal=custom_pal,
                          contig_seps=col_sep,
                          num_obs_groups=k_obs_groups,
                          obs_annotations_groups=obs_annotations_groups,
                          obs_annotations_names=obs_annotations_names,
                          grouping_key_coln=grouping_key_coln[1],
                          cluster_by_groups=cluster_by_groups,
                          cnv_title=title,
                          cnv_obs_title=obs_title,
                          contig_lab_size=contig_cex,
                          breaksList=breaksList_t,
                          x.center=x.center,
                          hclust_method=hclust_method,
                          layout_lmat=force_layout[["lmat"]],
                          layout_lhei=force_layout[["lhei"]],
                          layout_lwid=force_layout[["lwid"]])
    obs_data <- NULL

    if(!is.null(ref_idx)){
        .plot_cnv_references(ref_data=ref_data_t,
                            ref_groups=ref_groups,
                            name_ref_groups=name_ref_groups,
                            grouping_key_coln=grouping_key_coln[2],
                            col_pal=custom_pal,
                            contig_seps=col_sep,
                            file_base_name=out_dir,
                            output_filename_prefix=output_filename,
                            cnv_ref_title=ref_title,
                            breaksList=breaksList_t,
                            x.center=x.center,
                            layout_add=TRUE)
    }
    if (! is.na(output_format)) {
        dev.off()
    }
}

# TODO Tested, test make files so turned off but can turn on and should pass.
# Plot the observational samples
#
# Args:
#' @param obs_data Data to plot as observations. Rows = Cells, Col = Genes.
#' @param col_pal The color palette to use.
#' @param contig_colors The colors for the contig bar.
#' @param contig_labels The labels for the contigs.
#' @param contig_names Names of the contigs.
#' @param contig_seps Indices for line seperators of contigs.
#' @param num_obs_groups Number of groups of observations to create.
#' @param file_base_name Base of the file to used to make output file names.
#' @param cnv_title Title of the plot.
#' @param cnv_obs_title Title for the observation matrix.
#' @param contig_lab_size Text size for contigs.
#' @param cluster_contig A value directs cluster to only genes on this contig.
#' @param obs_annotations_groups Vector of observations annotations group assignation (as indices).
#' @param cluster_by_groups Whether to cluster observations by their annotations or not. Using this ignores num_obs_groups.
#' @param breaksList List of values used as splitters on coloring range.
#' @param hclust_method Method to use for hclust.
#' @param testing If TRUE, does not plot anything.
#' @param layout_lmat lmat values to use in layout.
#' @param layout_lhei lhei values to use in layout.
#' @param layout_lwid lwid values to use in layout.
#
#' @return Void.
#'
#' @keywords internal
#' @noRd
#'

.plot_cnv_observations <- function(obs_data,
                                  col_pal,
                                  contig_colors,
                                  contig_labels,
                                  contig_names,
                                  contig_seps,
                                  num_obs_groups,
                                  file_base_name,
                                  output_filename_prefix,
                                  cnv_title,
                                  cnv_obs_title,
                                  contig_lab_size=1,
                                  cluster_contig=NULL,
                                  obs_annotations_groups,
                                  obs_annotations_names,
                                  grouping_key_coln,
                                  cluster_by_groups,
                                  breaksList,
                                  x.center,
                                  hclust_method="ward.D",
                                  testing=FALSE,
                                  layout_lmat=NULL,
                                  layout_lhei=NULL,
                                  layout_lwid=NULL){

    flog.info("plot_cnv_observation:Start")
    flog.info(paste("Observation data size: Cells=",
                    nrow(obs_data),
                    "Genes=",
                    ncol(obs_data),
                    sep=" "))
    observation_file_base <- paste(file_base_name, sprintf("%s.observations.txt", output_filename_prefix), sep=.Platform$file.sep)
    
    # Output dendrogram representation as Newick
    # Need to precompute the dendrogram so we can manipulate
    # it before the heatmap plot
    ## Optionally cluster by a specific contig
    hcl_desc <- "General"
    hcl_group_indices <- 1:ncol(obs_data)
    
    if(!is.null(cluster_contig)){
        # restricting to single contig
        hcl_contig_indices <- which(contig_names == cluster_contig)
        if(length(hcl_group_indices) > 0 ) {
            hcl_group_indices <- hcl_contig_indices
            hcl_desc <- cluster_contig
            flog.info(paste("plot_cnv_observation:Clustering only by contig ", cluster_contig))
        } else {
           flog.warn(paste("plot_cnv_observations: Not able to cluster by",
                                     cluster_contig,
                                     "Clustering by all genomic locations.",
                                     "To cluster by local genomic location next time",
                                     "select from:",
                                     unique(contig_names),
                                     collapse=",",
                                     sep=" "))
        }
    }
    
    # HCL with a inversely weighted euclidean distance.
    flog.info(paste("clustering observations via method: ", hclust_method, sep=""))
    # obs_hcl <- NULL
    obs_dendrogram <- list()
    ordered_names <- NULL
    isfirst <- TRUE
    hcl_obs_annotations_groups <- vector()
    obs_seps <- c()

    
    if (cluster_by_groups) {

        ## Clustering separately by groups (ie. patients)

        for (i in seq(1, max(obs_annotations_groups))) {
            gene_indices_in_group <- which(obs_annotations_groups == i)
            num_genes_in_group <- length(gene_indices_in_group)
            flog.info(sprintf("Number of genes in group(%d) is %d", i, num_genes_in_group))

            if (num_genes_in_group < 2) {
                flog.info(sprintf("Skipping group: %d, since less than 2 entries", i))
                next
            }
            
            data_to_cluster <- obs_data[gene_indices_in_group, hcl_group_indices, drop=F]
            flog.info(paste("group size being clustered: ", paste(dim(data_to_cluster), collapse=","), sep=" "))
            group_obs_hcl <- hclust(dist(data_to_cluster), method=hclust_method)
            ordered_names <- c(ordered_names, row.names(obs_data[which(obs_annotations_groups == i), hcl_group_indices])[group_obs_hcl$order])

                           
            if (isfirst) {
                write.tree(as.phylo(group_obs_hcl),
                           file=paste(file_base_name, sprintf("%s.observations_dendrogram.txt", output_filename_prefix), sep=.Platform$file.sep))
                isfirst <- FALSE
            }
            else {
                write.tree(as.phylo(group_obs_hcl),
                           file=paste(file_base_name, sprintf("%s.observations_dendrogram.txt", output_filename_prefix), sep=.Platform$file.sep), append=TRUE)
            }
            
            group_obs_dend <- as.dendrogram(group_obs_hcl)
            obs_dendrogram[[length(obs_dendrogram) + 1]] <- group_obs_dend
            hcl_obs_annotations_groups <- c(hcl_obs_annotations_groups, rep(i, length(which(obs_annotations_groups == i))))
            obs_seps <- c(obs_seps, length(ordered_names))
        }
        if (length(obs_dendrogram) > 1) {
            # merge separate dendrograms into a single dendrogram
            obs_dendrogram <- do.call(merge, obs_dendrogram)
        } else {
            obs_dendrogram <- obs_dendrogram[[1]]
        }
        split_groups <- rep(1, dim(obs_data)[1])
        names(split_groups) <- ordered_names
    }
    else {
        # clustering all groups together
        obs_hcl <- hclust(dist(obs_data[,hcl_group_indices]), method=hclust_method)
                                        
        write.tree(as.phylo(obs_hcl),
                   file=paste(file_base_name, sprintf("%s.observations_dendrogram.txt", output_filename_prefix), sep=.Platform$file.sep))
        
        obs_dendrogram <- as.dendrogram(obs_hcl)
        ordered_names <- row.names(obs_data)[obs_hcl$order]
        split_groups <- cutree(obs_hcl, k=num_obs_groups)
        split_groups <- split_groups[ordered_names]
        hcl_obs_annotations_groups <- obs_annotations_groups[obs_hcl$order]
        
        # Make a file of members of each group
        flog.info("plot_cnv_observation:Writing observations by grouping.")
        for (cut_group in unique(split_groups)){
          group_memb <- names(split_groups)[which(split_groups == cut_group)]
          # Write group to file
          memb_file <- file(paste(file_base_name,
                                  paste(hcl_desc,"HCL",cut_group,"members.txt",sep="_"),
                                  sep=.Platform$file.sep))
          write.table(obs_data[group_memb,], memb_file)
          # Record seperation
          ordered_memb <- which(ordered_names %in% group_memb)
          if (is.null(obs_seps)) {
            obs_seps <- c(length(ordered_memb))
          }
          else {
            obs_seps <- c(obs_seps, (obs_seps[length(obs_seps)] + length(ordered_memb)))
          }
        }
        obs_seps <- c(obs_seps, length(ordered_names))
    }

    if (length(obs_seps) > 1) {
        obs_seps <- obs_seps[length(obs_seps)] - obs_seps[(length(obs_seps) - 1):1]
    }

    # Output HCL group membership.
    # Record locations of seperations

    # Make colors based on groupings
    row_groupings <- get_group_color_palette()(length(table(split_groups)))[split_groups]
    row_groupings <- cbind(row_groupings, get_group_color_palette()(length(table(hcl_obs_annotations_groups)))[hcl_obs_annotations_groups])
    annotations_legend <- cbind(obs_annotations_names, get_group_color_palette()(length(table(hcl_obs_annotations_groups))))

    # Make a file of coloring and groupings
    flog.info("plot_cnv_observation:Writing observation groupings/color.")
    groups_file_name <- file.path(file_base_name, sprintf("%s.observation_groupings.txt", output_filename_prefix))
    # file_groups <- rbind(split_groups,row_groupings)
    file_groups <- cbind(split_groups, row_groupings[,1], hcl_obs_annotations_groups, row_groupings[,2])
    #row.names(file_groups) <- c("Group","Color")
    colnames(file_groups) <- c("Dendrogram Group", "Dendrogram Color", "Annotation Group", "Annotation Color")
    # write.table(t(file_groups), groups_file_name)
    write.table(file_groups, groups_file_name)

    # Generate the Sep list for heatmap.3
    contigSepList <- create_sep_list(row_count=nrow(obs_data),
                                     col_count=ncol(obs_data),
                                     row_seps=obs_seps,
                                     col_seps=contig_seps)

    # reorder expression matrix based on dendrogram ordering
    obs_data <- obs_data[ordered_names, ]

    # Remove row/col labels, too cluttered
    # and print.
    orig_row_names <- row.names(obs_data)
    row.names(obs_data) <- rep("", nrow(obs_data))

    data_observations <- heatmap.cnv(obs_data,
                                        Rowv=obs_dendrogram,
                                        Colv=FALSE,
                                        cluster.by.row=TRUE,
                                        cluster.by.col=FALSE,
                                        main=cnv_title,
                                        ylab=cnv_obs_title,
                                        margin.for.labCol=2,
                                        xlab="Genomic Region",
                                        key=TRUE,
                                        labCol=contig_labels,
                                        cexCol=contig_lab_size,
                                        notecol="black",
                                        density.info="histogram",
                                        denscol="blue",
                                        trace="none",
                                        dendrogram="row",
                                        cexRow=0.8,
                                        breaks=breaksList,
                                        scale="none",
                                        x.center=x.center,
                                        color.FUN=col_pal,
                                        if.plot=!testing,
                                        # Seperate by contigs
                                        sepList=contigSepList,
                                        sep.color=c("black","black"),
                                        sep.lty=1,
                                        sep.lwd=1,
                                        # Color rows by user defined cut
                                        RowIndividualColors=row_groupings,
                                        annotations_legend=annotations_legend,
                                        grouping_key_coln=grouping_key_coln,
                                        # Color by contigs
                                        ColIndividualColors=contig_colors,
                                        # Legend
                                        key.title="Distribution of Expression",
                                        key.xlab="Modified Expression",
                                        key.ylab="Count",
                                        # Layout
                                        force_lmat=layout_lmat,
                                        force_lwid=layout_lwid,
                                        force_lhei=layout_lhei)

    # Write data to file.
    if (class(obs_data) %in% c("matrix", "data.frame")) {
        flog.info(paste("plot_cnv_references:Writing observation data to",
                        observation_file_base,
                        sep=" "))
        row.names(obs_data) <- orig_row_names
        write.table(t(obs_data[data_observations$rowInd,data_observations$colInd]),
                    file=observation_file_base)
    }
}


#' Not Testing, params ok.
#' Create the layout for the plot
#' This is a modification of the original
#' layout from the GMD heatmap.3 function
#'
#' Returns:
#' list with slots "lmat" (layout matrix),
#'                             "lhei" (height, numerix vector),
#'                             and "lwid" (widths, numeric vector)
#'
#' @keywords internal
#' @noRd
#'
.plot_observations_layout <- function(grouping_key_height)
{
    ## Plot observational samples
    obs_lmat <- c(0,  0,  0,  0,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,
                  # 8,  0, 10,  0,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,
                  8, 11, 10,  0,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  # 9 = reference heatmap
                  0,  0,  0,  0,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
                  5,  2,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  5,  2,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  5,  2,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  # 1 observations heatmap
                  5,  2,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  5,  2,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  5,  2,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  5,  2,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  5,  2,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  5,  2,  3,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
                  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,  6,
                  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
    obs_lmat <- matrix(obs_lmat,ncol=14,byrow=TRUE)

    # vector of heights for each of the above rows.
    obs_lhei <- c(1.125, 2.215, .15,
                   .5, .5, .5,
                   .5, .5, .5,
                   .5, .5, .5,
                  0.03, grouping_key_height[1]+0.04, grouping_key_height[2]+0.04, 0.03)

    return(list(lmat=obs_lmat,
           lhei=obs_lhei,
           lwid=NULL))
}

# TODO Tested, test make files so turned off but can turn on and should pass.
# Plot the reference samples
#
#' Args:
#' ref_data Data to plot as references. Rows = Cells, Col = Genes
#' ref_groups Groups of references to plot together.
#' col_pal The color palette to use.
#' contig_seps Indices for line seperators of contigs.
#' file_base_name Base of the file to used to make output file names.
#' cnv_ref_title Title for reference matrix.
#' layout_lmat lmat values to use in the layout.
#' layout_lwid lwid values to use in the layout.
#' layout_lhei lhei values to use in the layout.
#' layout_add Indicates the ref image shoudl be added to the previous plot.
#' testing: Turns off plotting when true.
#'
#' Returns:
#' Void
#'
#' @keywords internal
#' @noRd
#'
.plot_cnv_references <- function(ref_data,
                                ref_groups,
                                name_ref_groups,
                                grouping_key_coln,
                                col_pal,
                                contig_seps,
                                file_base_name,
                                output_filename_prefix,
                                cnv_ref_title,
                                breaksList,
                                x.center=x.center,
                                layout_lmat=NULL,
                                layout_lwid=NULL,
                                layout_lhei=NULL,
                                layout_add=FALSE,
                                testing=FALSE){

    flog.info("plot_cnv_references:Start")
    flog.info(paste("Reference data size: Cells=",
                           ncol(ref_data),
                           "Genes=",
                           nrow(ref_data),
                           sep=" "))
    number_references <- ncol(ref_data)
    reference_ylab <- NA
    reference_data_file <- paste(file_base_name, sprintf("%s.references.txt", output_filename_prefix), sep=.Platform$file.sep)
    
    ref_seps <- c()
    # Handle only one reference
    # heatmap3 requires a 2 x 2 matrix, so with one reference
    # I just duplicate the row and hid the second name so it
    # visually looks like it is just taking up the full realestate.
    if(number_references == 1){
        ref_data <- cbind(ref_data, ref_data)
        names(ref_data) <- c("",names(ref_data)[1])
        colnames(ref_data) <- c("", colnames(ref_data)[1])
    }

    # Handle reference groups
    # If there is more than one reference group, visually break
    # up the groups with a row seperator. Also plot the rows in
    # order so the current groups are shown and seperated.
    if(length(ref_groups) > 1){
        i_cur_idx <- 0
        order_idx <- c()
        grouping <- 0
        for (ref_grp in ref_groups){
            i_cur_idx <- i_cur_idx + length(ref_grp)
            ref_seps <- c(ref_seps, i_cur_idx)
            order_idx <- c(order_idx, ref_grp)
        }
        ref_seps <- ref_seps[1:(length(ref_seps) - 1)]
        ref_data <- ref_data[, order_idx, drop=FALSE]
    }

    split_groups <- c()
    i_cur_idx <- 1
    for (ref_grp in ref_groups){
        split_groups <- c(split_groups, rep(i_cur_idx, length(ref_grp)))
        i_cur_idx <- i_cur_idx + 1
    }

    # Make row color column colors from groupings
    flog.info(paste("plot_cnv_references:Number reference groups=",
                           length(ref_groups)),
                           sep=" ")

    # Transpose data.
    ref_data <- t(ref_data)
    # Remove labels if too many.
    ref_orig_names <- row.names(ref_data)
    if (number_references > 20){
        # The reference labs can become clustered
        # Dynamically change labels given a certain number of labels.
        reference_ylab <- cnv_ref_title
        row.names(ref_data) <- rep("", number_references)
    }

    row_groupings <- as.matrix(get_group_color_palette()(length(table(split_groups)))[split_groups])
    annotations_legend <- cbind(name_ref_groups, get_group_color_palette()(length(name_ref_groups)))

    # Generate the Sep list for heatmap.3
    contigSepList <- create_sep_list(row_count=nrow(ref_data),
                                     col_count=ncol(ref_data),
                                     row_seps=ref_seps,
                                     col_seps=contig_seps)

    # Print controls
    flog.info("plot_cnv_references:Plotting heatmap.")
    data_references <- heatmap.cnv(ref_data,
                                   main=NULL, #NA,
                                   ylab=reference_ylab,
                                   xlab=NULL, #NA,
                                   key=FALSE,
                                   labCol=rep("", nrow(ref_data)),
                                   notecol="black",
                                   trace="none",
                                   dendrogram="none",
                                   Colv=FALSE,
                                   Rowv=FALSE,
                                   cexRow=0.4,
                                   breaks=breaksList,
                                   scale="none",
                                   x.center=x.center,
                                   color.FUN=col_pal,
                                   sepList=contigSepList,
                                   # Row colors and legend
                                   RowIndividualColors=row_groupings,
                                   annotations_legend=annotations_legend,
                                   grouping_key_coln=grouping_key_coln,
                                   # Seperate by contigs
                                   sep.color=c("black","black"),
                                   sep.lty=1,
                                   sep.lwd=1,
                                   if.plot=!testing,
                                   # Layout
                                   force_lmat=layout_lmat,
                                   force_lwid=layout_lwid,
                                   force_lhei=layout_lhei,
                                   force_add=layout_add)

                                        # Write data to file
    if (class(ref_data) %in% c("matrix", "data.frame")) {

        ## TODO: write files for dgcMatrix too.
        row.names(ref_data) <- ref_orig_names
        flog.info(paste("plot_cnv_references:Writing reference data to",
                        reference_data_file,
                        sep=" "))
        write.table(t(ref_data[data_references$rowInd,data_references$colInd]),
                    file=reference_data_file)
    }
}




#########################################
# Plotting from GMD
###################
###
# This package relies on heatmap.3 from the library GMD
# This code is a modified version of heatmap.3, some changes
# were required for our plotting.
# The heatmap.cnv function should be considered a modification
# of th GMD library function heatmap.3, all credit goes to
# their authors.
## Please note this code is from the library GMD
## All credit for this code goes to GMD's authors.
## I do not recommend using this version of the code, which
## has been poorly modified for our use but recommend using
## the official version from the package GMD
## https://cran.r-project.org/web/packages/GMD/index.html
## A copy of gtools::invalid
##
## see \code{invalid} in package:gtools for details
## Test if a value is missing, empty, or contains only NA or NULL values
## param: x value to be tested
.invalid <-
  function(x)
{
  if (missing(x) || is.null(x) || length(x) == 0)
    return(TRUE)
  if (is.list(x))
    return(all(sapply(x, .invalid)))
  else if (is.vector(x))
    return(all(is.na(x)))
  else return(FALSE)
}


## Please note this code is from the library GMD
## All credit for this code goes to GMD's authors.
## I do not recommend using this version of the code, which
## has been poorly modified for our use but recommend using
## the official version from the package GMD
## https://cran.r-project.org/web/packages/GMD/index.html
.is.grouped <-
  function(x)
{
  x <- as.character(x)
  x.f <- factor(x,levels=unique(x),ordered=TRUE)
  identical(as.character(sort(x.f)),x)
}


## Please note this code is from the library GMD
## All credit for this code goes to GMD's authors.
## I do not recommend using this version of the code, which
## has been poorly modified for our use but recommend using
## the official version from the package GMD
## https://cran.r-project.org/web/packages/GMD/index.html
## Call a function with arguments
##
## Call a function with arguments
## FUN function or function name
## ... unnameed function arguments
## MoreArgs named (or unnameed) function arguments
.call.FUN <-
  function(FUN,...,MoreArgs)
{
  FUN <- match.fun(FUN)
  tmp.MoreArgs <- list(...)
  if (!.invalid(MoreArgs)){
    if (length(MoreArgs)>=1){
      for (i in 1:length(MoreArgs)) tmp.MoreArgs[[names(MoreArgs)[i]]] <- MoreArgs[[i]]
    }
  }
  ret <- do.call(FUN, tmp.MoreArgs)
  if ("call" %in% names(ret)){
    ret$call <- match.call()
  }
  if ("call" %in% names(attributes(ret))){
    attr(ret,"call") <- match.call()
  }
  return(ret)
}

## Please note this code is from the library GMD
## All credit for this code goes to GMD's authors.
## I do not recommend using this version of the code, which
## has been poorly modified for our use but recommend using
## the official version from the package GMD
## https://cran.r-project.org/web/packages/GMD/index.html
## Scale values to make them follow Standard Normal Distribution
##
## Scale values to make them follow Standard Normal Distribution
## param x numeric
## param scale character, indicating the type to scale.
## param na.rm logical
## return an object with the same dimention of `x'.
.scale.data <-
  function(x,scale,na.rm=TRUE)
{
  if(scale=="row"){
    x <- sweep(x,1,rowMeans(x,na.rm=na.rm),FUN="-")
    sx <- apply(x,1,sd,na.rm=na.rm)
    x <- sweep(x,1,sx,FUN="/")
  } else if(scale=="column"){
    x <- sweep(x,2,colMeans(x,na.rm=na.rm),FUN="-")
    sx <- apply(x,2,sd,na.rm=na.rm)
    x <- sweep(x,2,sx,,FUN="/")
  }
  x
}

## Please note this code is from the library GMD
## All credit for this code goes to GMD's author's.
## I do not recommend using this version of the code, which
## has been poorly modified for our use but recommend using
## the official version from the package GMD
## https://cran.r-project.org/web/packages/GMD/index.html
## Scale values to a new range: c(low, high)
## x numeric
## low numeric, lower bound of target values
## high numeric, higher bound of target values
## return an object with the same dimention of `x'.
.scale.x <-
  function(x,low=0,high=1,na.rm=TRUE)
{
  if(identical(max(x,na.rm=na.rm),min(x,na.rm=na.rm))) NA
  a <- 1/(max(x)-min(x))
  b <- -min(x)/(max(x)-min(x))
  a*x+b
}

## Please note this code is from the library GMD
## All credit for this code goes to GMD's author's.
## I do not recommend using this version of the code, which
## has been poorly modified for our use but recommend using
## the official version from the package GMD
## https://cran.r-project.org/web/packages/GMD/index.html
## Plot text
##
## Plot text
## x character, text to plot
## cex
## forecolor color of foreground
## bg color of background
## bordercolor color of border
## axes as in \code{graphics:::plot}
## ... additional arguments for \code{graphics:::text}
.plot.text <- function(x,xlim=c(0,1),ylim=c(0,1),cex=1,forecolor=par("fg"),bg=par("bg"),bordercolor=NA,axes=FALSE,...){
  if (.invalid(x)){
    x <- NULL
  }
  if (is.null(x)){
    x <- ""
  } else if (is.na(x)){
    x <- 'NA'
  }

  plot(xlim,ylim,type="n",ylab="",xlab="",xaxt="n",yaxt="n",axes=axes)
  rect(xleft=0, ybottom=0, xright=1, ytop=1, col=bg, border=bordercolor)
  text(0.5,0.5,x,cex=cex,...)
}

## Please note this code is from the library GMD
## All credit for this code goes to GMD's author's.
## I do not recommend using this version of the code, which
## has been poorly modified for our use but recommend using
## the official version from the package GMD
## https://cran.r-project.org/web/packages/GMD/index.html
## This was originally heatmap.3.
heatmap.cnv <-
  function(x,

           ## whether a dissimilarity matrix
           diss=inherits(x,"dist"),

           ## dendrogram control
           Rowv=TRUE,
           Colv=TRUE,
           dendrogram=c("both","row","column","none"),

           ## dist object
           dist.row,
           dist.col,
           dist.FUN=gdist,
           dist.FUN.MoreArgs=list(method="euclidean"),

           ## hclust object
           hclust.row,
           hclust.col,
           hclust.FUN=hclust,
           hclust.FUN.MoreArgs=list(method="ward.D"),

           ## data scaling
           scale=c("none","row","column"),
           na.rm=TRUE,

           ## clustering control
           cluster.by.row=FALSE,
           cluster.by.col=FALSE,
           kr=NA,
           kc=NA,
           row.clusters=NA,
           col.clusters=NA,

           ## image plot
           revR=FALSE,
           revC=FALSE,
           add.expr,

           ## mapping data to colors
           breaks,
           ## centering colors to a value
           x.center,
           ## colors
           color.FUN=gplots::bluered,
           ##
           ## block sepration
           sepList=list(NULL,NULL),
           sep.color=c("gray45","gray45"),
           sep.lty=1,
           sep.lwd=2,

           ## cell labeling
           cellnote,
           cex.note=1.0,
           notecol="cyan",
           na.color=par("bg"),

           ## level trace
           trace=c("none","column","row","both"),
           tracecol="cyan",
           hline,
           vline,
           linecol=tracecol,

           ## Row/Column Labeling
           labRow=TRUE, ## shown by default
           labCol=TRUE, ## shown by default
           srtRow=NULL,
           srtCol=NULL,
           sideRow=4,
           sideCol=1,
           ##
           margin.for.labRow,
           margin.for.labCol,
           ColIndividualColors,
           RowIndividualColors,
           annotations_legend,
           grouping_key_coln,
           cexRow,
           cexCol,
           labRow.by.group=FALSE,
           labCol.by.group=FALSE,


           ## plot color key + density info
           key=TRUE,
           key.title="Color Key",
           key.xlab="Value",
           key.ylab="Count",

           keysize=1.5,
           mapsize=9,
           mapratio=4/3,
           sidesize=3,
           cex.key.main=0.75,
           cex.key.xlab=0.75,
           cex.key.ylab=0.75,
           density.info=c("histogram","density","none"),
           denscol=tracecol,
           densadj=0.25,

           ## plot titles/labels
           main="Heatmap",
           sub="",
           xlab="",
           ylab="",
           cex.main=2,
           cex.sub=1.5,
           font.main=2,
           font.sub=3,
           adj.main=0.5,
           mgp.main=c(1.5,0.5,0),
           mar.main=3,
           mar.sub=3,
           ## plot ##
           if.plot=TRUE,

           ## plot of partition (left/top of heatmap)
           plot.row.partition=FALSE,
           plot.col.partition=FALSE,
           cex.partition=1.25,
           color.partition.box="gray45",
           color.partition.border="#FFFFFF",

           ## plot of summary (right/bottom of heatmap)
           plot.row.individuals=FALSE,
           plot.col.individuals=FALSE,
           plot.row.clusters=FALSE,
           plot.col.clusters=FALSE,
           plot.row.clustering=FALSE,
           plot.col.clustering=FALSE,

           ##
           plot.row.individuals.list=FALSE,
           plot.col.individuals.list=FALSE,
           plot.row.clusters.list=FALSE,
           plot.col.clusters.list=FALSE,
           plot.row.clustering.list=FALSE,
           plot.col.clustering.list=FALSE,

           ## for plot of clusters - row
           row.data=FALSE,
           ## for plot of clusters - col
           col.data=FALSE,

           ##
           if.plot.info=FALSE,
           text.box,
           cex.text=1.0,

           ## Force in layout info
           force_lmat=NULL,
           force_lwid=NULL,
           force_lhei=NULL,
           force_add=FALSE,

           ## extras
           ...
           )
{
  ## check input - take1 ##
  if (! is.matrix(x)) {
    x <- as.matrix(x)
  }
  x.ori <- x

  if(!inherits(x,"dist") & !is.matrix(x)){
    stop("`x' should either be a matrix, a data.frame or a `dist' object.")
  }

  if (! sideRow %in% c(2,4)){
    stop('sideRow must be either 2 or 4.')
  }

  if (! sideCol %in% c(1,3)){
    stop('sideCol must be either 1 or 3.')
  }

  ## store input
  Rowv.ori <- Rowv
  Colv.ori <- Colv

  ## check
  dendrogram <- match.arg(dendrogram)
  if ( (dendrogram %in% c("both","row")) & !inherits(Rowv,"dendrogram") ){
    warning("Discrepancy: row dendrogram is asked;  Rowv is set to `TRUE'.")
    Rowv <- TRUE
  }

  if ( (dendrogram %in% c("both","col")) & !inherits(Colv,"dendrogram") ){
    warning("Discrepancy: col dendrogram is asked;  Colv is set to `TRUE'.")
    Colv <- TRUE
  }

  if (identical(Rowv, FALSE) | missing(Rowv)){
    if(!identical(cluster.by.row,FALSE)){
      warning("Discrepancy: No row dendrogram is asked; cluster.by.row is set to `FALSE'.")
      cluster.by.row <- FALSE
    }
  } else {
    if(!identical(cluster.by.row,TRUE)){
      warning("Discrepancy: row dendrogram is asked; cluster.by.row is set to `TRUE'.")
      cluster.by.row <- TRUE
    }
  }

  if (identical(Colv, FALSE) | .invalid(Colv)){
    if(!identical(cluster.by.col,FALSE)){
      warning("Discrepancy: No col dendrogram is asked; cluster.by.col is set to `FALSE'.")
      cluster.by.col <- FALSE
    }
  } else {
    if(!identical(cluster.by.col,TRUE)){
      warning("Discrepancy: col dendrogram is asked; cluster.by.col is set to `TRUE'.")
      cluster.by.col <- TRUE
    }
  }

  if (!.invalid(kr)){
    if (is.numeric(kr)){
      if(!plot.row.partition){
        warning("Discrepancy: kr is set, therefore plot.row.partition is set to `TRUE'.")
        plot.row.partition <- TRUE
      }
    }
  }

  if (!.invalid(kc)){
    if (is.numeric(kc)){
      if(!plot.col.partition){
        warning("Discrepancy: kc is set, therefore plot.col.partition is set to `TRUE'.")
        plot.col.partition <- TRUE
      }
    }
  }

  ## generate dist.obj - row/col ##
  if (inherits(x,"dist")){
    dist.row <- dist.col <- x ## dist.obj
    x <- as.matrix(x)
    mat.row <- mat.col <- x ## mat.obj
    symm <- TRUE
  } else if (is.matrix(x)){
    symm <- isSymmetric(x)
    if (diss){
      if (!symm){
        stop("Dissimilary matrix should be symmetric. Please set `diss' to FALSE if `x' is not dissimilary matrix.")
      } else {
        flush.console()
      }
      mat.row <- mat.col <- x
      dist.row <- dist.col <- as.dist(x)
    } else{
      if (cluster.by.row) {
        if (.invalid(dist.row)){
          dist.row <- .call.FUN(dist.FUN,x,MoreArgs=dist.FUN.MoreArgs)
        }
        mat.row <- as.matrix(dist.row)
      } else {
        dist.row <- NULL
        mat.row <- NULL
      }
      if (cluster.by.col) {
        if (.invalid(dist.col)){
          dist.col <- .call.FUN(dist.FUN,t(x),MoreArgs=dist.FUN.MoreArgs)
        }
        mat.col <- as.matrix(dist.col)
      } else {
        dist.col <- NULL
        mat.col <- NULL
      }
    }
  }
  ## check input - take2: di ##
  di <- dim(x)

  if(length(di)!=2 || !is.numeric(x)){
    stop("`x' should only contain `numeric' values and can be converted to a 2-D matrix.")
  }

  ## parse param ##
  scale <- if(symm && .invalid(scale)) "none" else match.arg(scale) ## no scale on symmetric matrix
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  dist.FUN <- match.fun(dist.FUN)
  hclust.FUN <- match.fun(hclust.FUN)
  color.FUN <- match.fun(color.FUN)

  ## NG if both breaks and scale are specified ##
  if(!.invalid(breaks) & (scale!="none")){
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.",
            "Please consider using only one or the other.")
  }

  ## nr and nc ##
  nr <- di[1]
  nc <- di[2]

  ## check input - take3: nr,nc ##
  if(nr <=1 || nc <=1)
    stop("`x' must have at least 2 rows and 2 columns")

  ## font size of row/col labels ##
  cexRow0 <- 0.2+1/log10(nr)
  cexCol0 <- 0.2+1/log10(nc)

  if (.invalid(cexRow)) {
    cexRow <- cexRow0
  } else {
    # message('heatmap.3 | From GMD 0.3.3, please use relative values for cexRow.')
    cexRow <- cexRow0*cexRow
  }
  if (.invalid(cexCol)) {
    cexCol <- cexCol0
  } else {
    # message('heatmap.3 | From GMD 0.3.3, please use relative values for cexCol.')
    cexCol <- cexCol0*cexCol
  }

  ## cellnote ##
  ## ##if(.invalid(cellnote)) cellnote <- matrix("",ncol=ncol(x),nrow=nrow(x))

  ## ------------------------------------------------------------------------
  ## parse dendrogram ##
  ## ------------------------------------------------------------------------

  if (missing(Rowv)) Rowv <- FALSE
  if (.invalid(Colv)) Colv <- if(symm) Rowv else FALSE
  if (Colv=="Rowv") {
    if ((!isTRUE(Rowv) | !symm) ){
      Colv <- FALSE
      warning("`Colv' is specified to use \"Rowv\", but either `Rowv' is invalid or `x' is not symmetric; Colv is suppressed.")
    } else{
      Colv <- Rowv
    }
  }
  ## ------------------------------------------------------------------------
  ## generate hclust.obj - row/col
  ## ------------------------------------------------------------------------
  flush.console()

  if ( (!inherits(Rowv,"dendrogram") & !identical(Rowv,FALSE)) | (cluster.by.row & .invalid(row.clusters))){
    if (.invalid(hclust.row)){
      hclust.row <- .call.FUN(hclust.FUN,dist.row,MoreArgs=hclust.FUN.MoreArgs)
    } else {
      if (length(hclust.row$order) != nr){
        stop("`hclust.row' should have equal size as the rows.")
      }
    }
  } else{
    hclust.row <- NULL
  }

  if(symm){
    hclust.col <- hclust.row
  }

  if ( !inherits(Colv,"dendrogram") & !identical(Colv,FALSE) | (cluster.by.col & .invalid(col.clusters))){
    if (.invalid(hclust.col)){
      hclust.col <- .call.FUN(hclust.FUN,dist.col,MoreArgs=hclust.FUN.MoreArgs)
    } else {
      if (length(hclust.col$order) != nc){
        stop("`hclust.col' should have equal size as the cols.")
      }
    }
  } else {
    hclust.col <- NULL
  }
  ## ------------------------------------------------------------------------
  ## generate hclust.obj - row/col
  ## ------------------------------------------------------------------------
  ddr <- ddc <- NULL
  ## get the dendrograms and reordering row/column indices - row ##
  if(inherits(Rowv,"dendrogram")){
    if (attr(Rowv,"members") != nr){
      stop("`Rowv' should contain equal size of members as the rows.")
    }
    ddr <- Rowv ## use Rowv 'as-is',when it is dendrogram
    # rowInd <- order.dendrogram(ddr)
    rowInd <- 1:nr  ### data has already been sorted in the proper order, and coloring vectors too
  } else if (is.integer(Rowv)){ ## compute dendrogram and do reordering based on given vector
    ddr <- as.dendrogram(hclust.row)
    ddr <-  reorder(ddr,Rowv)
    rowInd <- order.dendrogram(ddr)
    if(nr != length(rowInd)){
      stop("`rowInd' is of wrong length.")
    }
  } else if (isTRUE(Rowv)){ ## if TRUE,compute dendrogram and do reordering based on rowMeans

    Rowv <- rowMeans(x,na.rm=TRUE)
    ddr <- as.dendrogram(hclust.row)
    ddr <- reorder(ddr,Rowv)
    rowInd <- order.dendrogram(ddr)
    if(nr !=length(rowInd)){
      stop("`rowInd' is of wrong length.")
    }
  } else{
    rowInd <- nr:1 ## from bottom.
  }
  ## get the dendrograms and reordering row/column indices - col ##
  if(inherits(Colv,"dendrogram")){
    if (attr(Colv,"members") != nc){
      stop("`Colv' should contain equal size of members as the cols.")
    }
    ddc <- Colv ## use Colv 'as-is',when it is dendrogram
    colInd <- order.dendrogram(ddc)
  } else if(identical(Colv,"Rowv")) {
    if(exists("ddr")){
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    } else{
      colInd <- rowInd
    }
  } else if(is.integer(Colv)){## compute dendrogram and do reordering based on given vector
    ddc <- as.dendrogram(hclust.col)
    ddc <- reorder(ddc,Colv)
    colInd <- order.dendrogram(ddc)
    if(nc != length(colInd))
      stop("`colInd' is of wrong length.")
  } else if (isTRUE(Colv)){## if TRUE,compute dendrogram and do reordering based on rowMeans
    Colv <- colMeans(x,na.rm=TRUE)
    ddc <- as.dendrogram(hclust.col)
    ddc <- reorder(ddc,Colv)
    colInd <- order.dendrogram(ddc)
    if(nc !=length(colInd))
      stop("`colInd' is of wrong length.")
  } else{
    colInd <- 1:nc ## from left
  }
  ## ------------------------------------------------------------------------
  ## check consistency
  ## ------------------------------------------------------------------------

  ## Xmisc::logme(dendrogram)
  ## Xmisc::logme(Colv)
  ## Xmisc::logme(Rowv)

  ## dendrogram - check consistency: Rowv ##
  if ( is.null(ddr) & (dendrogram %in% c("both","row"))){
    warning("Discrepancy: Rowv is invalid or FALSE, while dendrogram is `",
            dendrogram,"'. Omitting row dendogram.")
    if (is.logical(Colv) & (Colv.ori) & dendrogram=="both")
      dendrogram <- "column"
    else
      dendrogram <- "none"
  }

  ## dendrogram - check consistency: Colv ##
  if ( is.null(ddc) & (dendrogram %in% c("both","column"))){
    warning("Discrepancy: Colv is invalid or FALSE, while dendrogram is `",
            dendrogram,"'. Omitting column dendogram.")
    if (is.logical(Rowv) & (identical(Rowv.ori,TRUE) | is.numeric(Rowv.ori) | inherits(Rowv.ori,"dendrogram")) & dendrogram=="both")
      dendrogram <- "row"
    else
      dendrogram <- "none"
  }

  ## check consistency
  if (is.null(ddr)){
    if(isTRUE(cluster.by.row) | isTRUE(plot.row.partition) | isTRUE(plot.row.clusters) | isTRUE(plot.row.clustering) ){
      warning("Using invalid `Rowv' while allowing",
              "`cluster.by.row' or `plot.row.partition' or `plot.row.clusters' or `plot.row.clustering'",
              "can produce unpredictable results; Forced to be disabled.")
    }
  }
  if (is.null(ddc)){
    if(isTRUE(cluster.by.col) | isTRUE(plot.col.partition) | isTRUE(plot.col.clusters) | isTRUE(plot.col.clustering) ){
      warning("Using invalid `Colv' while allowing",
              "`cluster.by.col' or `plot.col.partition' or `plot.col.clusters' or `plot.col.clustering'",
              "can produce unpredictable results; Forced to be disabled.")
    }
  }

  if (is.null(ddr)) cluster.by.row <- plot.row.partition <- plot.row.clusters <- plot.row.clustering <- FALSE
  if (is.null(ddc)) cluster.by.col <- plot.col.partition <- plot.col.clusters <- plot.col.clustering <- FALSE
  ## ------------------------------------------------------------------------
  ## Reordering
  ## ------------------------------------------------------------------------
  flush.console()
  ## reorder x and cellnote ##
  x <- x[rowInd,colInd]

  if (!.invalid(cellnote)) cellnote <- cellnote[rowInd,colInd]

  ## reorder labels - row ##
  if(identical(labRow,TRUE)){ ## Note: x is already reorderred
    labRow <- if (is.null(rownames(x))) (1:nr)[rowInd] else rownames(x)
  } else if(identical(labRow,FALSE) | .invalid(labRow)){
    labRow <- rep("",nrow(x))
  } else if(is.character(labRow)){
    labRow <- labRow[rowInd]
  }

  ## reorder cellnote/labels - col ##
  if (identical(labCol,TRUE)){
    labCol <- if(is.null(colnames(x))) (1:nc)[colInd] else colnames(x)
  } else if(identical(labCol,FALSE) | .invalid(labCol)){
    labCol <- rep("",ncol(x))
  } else if(is.character(labCol)){
    labCol <- labCol[colInd]
  }

  ## ------------------------------------------------------------------------
  ## scale
  ## center to 0 and scale to 1 in row or col but not both! ##
  ## ------------------------------------------------------------------------
  flush.console()
  x <- .scale.data(x,scale,na.rm)
  ## ------------------------------------------------------------------------
  ## labels for observations/clusters/
  ## ------------------------------------------------------------------------
  ## margin for labels

  margin.for.labRow0 <- max(nchar(labRow, keepNA = FALSE))*0.75+0.2
  margin.for.labCol0 <- max(nchar(labCol, keepNA = FALSE))*0.75+0.2  ## Change compared to GMD to work with R 3.3+ that changed the default value of keepNA to "NA"

  if (.invalid(margin.for.labRow)){
    margin.for.labRow <- margin.for.labRow0
  }
  #} else {
  #  message('heatmap.3 | From GMD 0.3.3, please use relative values for margin.for.labRow.')
  #  margin.for.labRow <- margin.for.labRow0*margin.for.labRow
  #}
  if (margin.for.labRow < 2) {
      margin.for.labRow <- 2
  }

  # if (.invalid(margin.for.labCol)){
  margin.for.labCol <- margin.for.labCol0
  # }
  #  else {
  #   # message('heatmap.3 | From GMD 0.3.3, please use relative values for margin.for.labCol.')
  #   margin.for.labCol <- margin.for.labCol0*margin.for.labCol
  # }

  if (margin.for.labCol < 2) {
    margin.for.labCol <- 2
  }

  ## group unique labels - row ## ##??check
  if (!.invalid(labRow.by.group) & !identical(labRow.by.group,FALSE)){
    group.value <- unique(labRow)
    group.index <- sapply(group.value,function(x,y) min(which(y==x)),y=labRow)
    labRow <- rep("",length(labRow))
    labRow[group.index] <- group.value
  }

  ## group unique labels - col ## ##??check
  if (!.invalid(labCol.by.group) & !identical(labCol.by.group,FALSE)){
    group.value <- unique(labCol)
    group.index <- sapply(group.value,function(x,y) min(which(y==x)),y=labCol)
    labCol <- rep("",length(labCol))
    labCol[group.index] <- group.value
  }
  ## ------------------------------------------------------------------------
  ## color breaks
  ## ------------------------------------------------------------------------
  flush.console()

  ## set breaks for binning x into colors ##
  if(.invalid(breaks)){
    breaks <- 16
  } else {
      flog.debug(paste("inferCNV::heatmap.cnv, breaks parameter set to: [", paste(breaks, collapse=","), "]", sep=""))
  }

  ## get x.range according to the value of x.center ##
  if (!.invalid(x.center)){ ## enhanced
    if (is.numeric(x.center)){
      x.range.old <- range(x,na.rm=TRUE)
      if (length(breaks) > 1) {
          # important, use specified breakpoint info here if set by user
          x.range.old = range(breaks)
      }
      dist.to.x.center <- max(abs(x.range.old-x.center))
      x.range <- c(x.center-dist.to.x.center,x.center+dist.to.x.center)
      if (length(breaks) > 1) {
          # re-set the breaks according to the new x.range
          breaks=seq(x.range[1], x.range[2], length=16)
          flog.debug(paste("inferCNV::heatmap.cnv, resetting breaks to adjusted x.range: [",
                                  paste(breaks, collapse=","), "]", sep=""))
      }

    } else {
      stop("`x.center' should be numeric.")
    }
  } else{
      x.range <- range(x,na.rm=TRUE)
      if (length(breaks) > 1) {
          # important, use specified breakpoint info here if set by user
          x.range = range(breaks)
      }

  }
  flog.debug( paste("inferCNV::heatmap.cnv x range set to: ",
                          paste(x.range, collapse=",")), sep="" )

  ## set breaks for centering colors to the value of x.center ##
  if(length(breaks)==1){
    breaks <-
      seq(min(min(x,na.rm=TRUE),x.range[1]),
          max(max(x,na.rm=TRUE),x.range[2]),
          length.out=breaks)
  }

  ## count of breaks and colors ##
  nbr <- length(breaks)
  ncolor <- length(breaks)-1

  ## generate colors ##
  colors <- color.FUN(ncolor)
  ## set up breaks and force values outside the range into the endmost bins ##
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[] <- ifelse(x<min.breaks,min.breaks,x)
  x[] <- ifelse(x>max.breaks,max.breaks,x)
  ## ------------------------------------------------------------------------
  ## check if it is sufficient to draw side plots ##
  ## ------------------------------------------------------------------------
  if (cluster.by.row){
    if (!.invalid(row.clusters)) {## suppress kr
      if(!is.numeric(row.clusters) | length(row.clusters)!=nr | !(.is.grouped(row.clusters))){
        warning("`row.clusters' is not a grouped numeric vector of length nrow(x); cluster.by.row is set to FALSE.")
        cluster.by.row <- FALSE
      } else{
        row.clusters <- row.clusters[rowInd]
        kr <- length(unique(row.clusters))
      }
    } else {
      if (.invalid(kr)) kr <- 2
      if (is.numeric(kr) & length(kr)==1){
        row.clusters <- cutree(hclust.row,k=kr)
        row.clusters <- row.clusters[rowInd]
      } else {
        warning("`kr' should be numeric of length one; cluster.by.row is set to FALSE.")
        cluster.by.row <- FALSE
      }
    }
  }

  if (cluster.by.col){
    if (!.invalid(col.clusters)) {## suppress kc
      if(!is.numeric(col.clusters) | length(col.clusters)!=nc | !(.is.grouped(col.clusters))){
        warning("`col.clusters' is not a grouped numeric vector of length ncol(x); cluster.by.col is set to FALSE.")
        cluster.by.col <- FALSE
      } else{
        col.clusters <- col.clusters[colInd]
        kc <- length(unique(col.clusters))
        if(revC){ ## x columns reversed
          col.clusters <- rev(col.clusters)
        }

      }
    } else {
      if (.invalid(kc)) kc <- 2
      if (is.numeric(kc) & length(kc)==1){
        col.clusters <- cutree(hclust.col,k=kc)
        col.clusters <- col.clusters[colInd]
        if(revC){ ## x columns reversed
          col.clusters <- rev(col.clusters)
        }

      } else {
        warning("`kc' should be numeric of length one; cluster.by.col is set to FALSE.")
        cluster.by.col <- FALSE
      }
    }
  }
  ## ------------------------------------------------------------------------
  ## Plotting
  ## ------------------------------------------------------------------------
  if (if.plot){

    ir <- length(plot.row.individuals.list)
    ic <- length(plot.col.individuals.list)
    cr <- length(plot.row.clustering.list)
    cc <- length(plot.col.clustering.list)

    flush.console()
    if(mapratio<=1){
      sr <- 12
      sc <- sr*mapratio
    } else {
      sc <- 12
      sr <- sc/mapratio
    }

    ###############################
    ## calculate the plot layout ##
    ###############################

    ## 1) for heatmap
    lmat <- matrix(1,nrow=sr,ncol=sc)
    lwid <- c(rep(mapsize/sc,sc))
    lhei <- c(rep(mapsize/mapratio/sr,sr))

    ## 2) row.clusters
    if (plot.row.partition | plot.row.clusters){
      lmat <- cbind(max(lmat,na.rm=TRUE)+1,lmat)
      lwid <- c(0.3,lwid)
    } else {
      lmat <- cbind(NA,lmat)
      lwid <- c(0.02,lwid)

    }
    ## 3) col.clusters
    if (plot.col.partition | plot.col.clusters){
      lmat <- rbind(c(NA,rep(max(lmat,na.rm=TRUE)+1,sc)),lmat)
      lhei <- c(0.3/mapratio,lhei)
    } else {
      lmat <- rbind(NA,lmat)
      lhei <- c(0.02/mapratio,lhei)
    }

    if(!.invalid(RowIndividualColors)) { ## 4) add middle column to layout for vertical sidebar ##??check
      # if(!is.character(RowIndividualColors) || length(RowIndividualColors) !=nr)
      if(!is.character(RowIndividualColors) || dim(RowIndividualColors)[1] !=nr)
        stop("'RowIndividualColors' must be a character vector of length nrow(x)")
      lmat <- cbind(c(rep(NA,nrow(lmat)-sr),rep(max(lmat,na.rm=TRUE)+1,sr)),lmat)
      lwid <- c(0.2,lwid)
      lmat <- cbind(c(rep(NA,nrow(lmat)-sr),rep(max(lmat,na.rm=TRUE)+1,sr)),lmat)  # not really useful as this if is only true in plot_cnv_observations which forces layout_lmat and layout_lhei
      lwid <- c(0.2,lwid) # layout_lwid is however not forced, so need to actually update it here
    } else {
      lmat <- cbind(NA,lmat)
      lwid <- c(0.02,lwid)
    }

    if(!.invalid(ColIndividualColors)) { ## 5) add middle row to layout for horizontal sidebar ##??check
      if(!is.character(ColIndividualColors) || length(ColIndividualColors) !=nc){
        stop("'ColIndividualColors' must be a character vector of length ncol(x)")
      }
      lmat <- rbind(c(rep(NA,ncol(lmat)-sc),rep(max(lmat,na.rm=TRUE)+1,sc)),lmat)
      lhei <- c(0.2/mapratio,lhei)
    } else {
      lmat <- rbind(NA,lmat)
      lhei <- c(0.02/mapratio,lhei)
    }
    ## 6) for row-dend
    lmat <- cbind(c(rep(NA,nrow(lmat)-sr),
                    rep(max(lmat,na.rm=TRUE)+1,sr)),
                  lmat
                  )
    lwid <- c(keysize,lwid)

    ## 7) for col-dend, 8) for key
    lmat <- rbind(c(
                    max(lmat,na.rm=TRUE)+2,
                    rep(NA,ncol(lmat)-sc-1),
                    rep(max(lmat,na.rm=TRUE)+1,sc)
                    ),
                  lmat
                  )
    lhei <- c(keysize/mapratio,lhei)
    ## text.box##
    ## numbered 999 ##
    ## 9) for RowPlot (from bottom)
    if(.invalid(text.box)){
      text.box <- "made by\nFunction: heatmap.3\nPackage: GMD\nin R"
    }
    if(plot.row.individuals) { ## enhanced: add right column to layout for plots
      lmat <- cbind(lmat,
                    c(rep((1+max(lmat,na.rm=TRUE)),nrow(lmat)-sr),# text
                      rep((ir:1)+max(lmat,na.rm=TRUE)+(1),each=floor(sr/ir)),rep(NA,sr%%ir)
                      )
                    )
      lwid <- c(lwid,sidesize)
    } else {
      lmat <- cbind(lmat,c(rep(NA,nrow(lmat))))
      lwid <- c(lwid,0.01)
    }

    ## 10) for ColPlot from right
    if(plot.col.individuals) { ## enhanced: add bottom row to layout for plots
      lmat <- rbind(lmat,
                    c(rep((1+max(lmat,na.rm=TRUE)),ncol(lmat)-sc-1),# text
                      rep((1:ic)+max(lmat,na.rm=TRUE)+(1),each=floor(sc/ic)),rep(NA,sc%%ic),
                      ##NA # change to numeric if text.box
                      999
                      )
                    )
      lhei <- c(lhei,sidesize/mapratio)
    } else {
      lmat <- rbind(lmat,c(rep(NA,ncol(lmat))))
      lhei <- c(lhei,0.01/mapratio)
    }
    ## 11) for RowPlot (from bottom)
    if(plot.row.clusters) { ## enhanced: add right column to layout for plots
      lmat <- cbind(lmat,
                    c(rep((1+max(lmat[lmat!=999],na.rm=TRUE)),nrow(lmat)-sr-1), # text
                      rep((kr:1)+max(lmat[lmat!=999],na.rm=TRUE)+(1),each=floor(sr/kr)),rep(NA,sr%%kr),
                      ##NA
                      999
                      )
                    )
      lwid <- c(lwid,sidesize)
    } else {
      lmat <- cbind(lmat,c(rep(NA,nrow(lmat))))
      lwid <- c(lwid,0.01)
    }

    ## 12) for ColPlot from right
    if(plot.col.clusters) { ## enhanced: add bottom row to layout for plots
      lmat <- rbind(lmat,
                    c(rep((1+max(lmat[lmat!=999],na.rm=TRUE)),ncol(lmat)-sc-2),# text
                      rep((1:kc)+max(lmat[lmat!=999],na.rm=TRUE)+(1),each=floor(sc/kc)),rep(NA,sc%%kc),
                      ##NA,NA # change to numeric if text.box
                      999,999
                      )
                    )
      lhei <- c(lhei,sidesize/mapratio)
    } else {
      lmat <- rbind(lmat,c(rep(NA,ncol(lmat))))
      lhei <- c(lhei,0.01/mapratio)
    }

    ## 13) for RowPlot (from bottom)
    if(plot.row.clustering) { ## enhanced: add right column to layout for plots
      lmat <- cbind(lmat,
                    c(rep((1+max(lmat[lmat!=999],na.rm=TRUE)),nrow(lmat)-sr-2), # text
                      rep(c((cr:1)+max(lmat[lmat!=999],na.rm=TRUE)+(1)),each=floor(sr/cr)),rep(NA,sr%%cr),
                      ##NA,NA
                      999,999
                      )
                    )
      lwid <- c(lwid,sidesize)
    } else {
      lmat <- cbind(lmat,c(rep(NA,nrow(lmat))))
      lwid <- c(lwid,0.01)
    }


    ## 14) for ColPlot from right
    if(plot.col.clustering) { ## enhanced: add bottom row to layout for plots
      lmat <- rbind(lmat,
                    c(rep((1+max(lmat[lmat!=999],na.rm=TRUE)),ncol(lmat)-sc-3),# text
                      rep((1:cc)+max(lmat[lmat!=999],na.rm=TRUE)+(1),each=floor(sc/cc)),rep(NA,sc%%cc),
                      ##NA,NA,NA # change to numeric if text.box
                      999,999,999
                      )
                    )
      lhei <- c(lhei,sidesize/mapratio)
    } else {
      lmat <- rbind(lmat,c(rep(NA,ncol(lmat))))
      lhei <- c(lhei,0.01/mapratio)
    }

    lmat[is.na(lmat)] <- 0
    if (any(lmat==999)) flag.text <- TRUE else flag.text <- FALSE
    lmat[lmat==999] <- max(lmat[lmat!=999])+1

    ##-------------------------------  Overriding the lmat settings if already given.
    ## Graphics `output' ##
    ## layout
    if(!is.null(force_lmat)){
        lmat <- force_lmat
    }
    if(!is.null(force_lwid)){
        lwid <- force_lwid
    }
    if(!is.null(force_lhei)){
        lhei <- force_lhei
    }
    if(!force_add){
        layout(lmat,widths=lwid,heights=lhei,respect=FALSE)
    }
    ##---------------------------------
      
    ## reverse columns
    if(revC){ ## x columns reversed
      iy <- nr:1
      ddc <- rev(ddc)
      x <- x[iy,]
      if (!.invalid(cellnote)) cellnote <- cellnote[iy,]
    } else {
      iy <- 1:nr
    }

    ## reverse rows
    if(revR){ ## x columns reversed
      ix <- nc:1
      ddr <- rev(ddr)
      x <- x[,ix]
      if (!.invalid(cellnote)) cellnote <- cellnote[,ix]
    } else {
      ix <- 1:nc
    }

    ## 1) draw the main carpet/heatmap
    margins <- c(margin.for.labCol,0,0,margin.for.labRow)
    mgp <- c(3,1,0)
    par(mar=margins,mgp=mgp);outer=FALSE

    x.save <- x
    if(!symm || scale !="none"){ ##??
      x <- t(x)
      if (!.invalid(cellnote)) cellnote <- t(cellnote)
    }
    image(1:nc,1:nr,
          x,
          xlim=0.5+c(0,nc),ylim=0.5+c(0,nr),
          axes=FALSE,xlab="",ylab="",col=colors,breaks=breaks,
          ...)
    flog.info(paste("Colors for breaks: ", paste(colors, collapse=","), sep=" "))
    flog.info(paste("Quantiles of plotted data range:", paste(quantile(x), collapse=","), sep=" "))
    
    ## plot/color NAs
    if(!.invalid(na.color) & any(is.na(x))){
      mmat <- ifelse(is.na(x),1,NA)
      image(1:nc,1:nr,mmat,axes=FALSE,xlab="",ylab="",
            col=na.color,add=TRUE)
    }

    ##
    ## labCol (?)
    if ((dendrogram %in% c("both","col")) & sideCol==3) {
      warning("Discrepancy: col dendrogram is asked; srtCol is set to 1.")
      sideCol <- 1
    }
    if (!length(srtCol)) {
      axis(sideCol,1:nc,labels=labCol,las=2,line=-0.5,tick=0,cex.axis=cexCol,outer=outer)
    } else {
      if (sideCol==1){
        if (sideCol==1) .sideCol <- par("usr")[3]-0.5*srtCol/90 else .sideCol <- par("usr")[4]+0.5*srtCol/90
        text(1:nc,.sideCol,labels=labCol,srt=srtCol,pos=1,xpd=TRUE,cex=cexCol)
      }
    }

    if(!.invalid(xlab)) mtext(xlab,side=1,line=margins[1]-1.25)

    ## labRow (?)
    if ((dendrogram %in% c("both","row")) & sideRow==2) {
      warning("Discrepancy: row dendrogram is asked; sideRow is set to 4.")
      sideRow <- 4
    }
    if (!length(srtRow)) {
      axis(sideRow,iy,labels=labRow,las=2,line=-0.5,tick=0,cex.axis=cexRow,outer=outer)
    } else {
      if (sideRow==4){
        if (sideRow==4) .sideRow <- par("usr")[2]+0.5*srtRow/90 else .sideRow <- par("usr")[1]-0.5*srtRow/90
        text(.sideRow,iy,labels=labRow,srt=srtRow,pos=1,xpd=TRUE,cex=cexRow)
      }
    }

    if(!.invalid(ylab)) mtext(ylab,side=4,line=margins[4]-1.25)

    if (!.invalid(add.expr))
      eval(substitute(add.expr))
    ## Enhanced: add 'sep.color' colored spaces to visually separate sections
    if (plot.row.partition | plot.row.clusters){ ##??
      plot.row.partitionList <- get.sep(clusters=row.clusters,type="row")
    } else {
      plot.row.partitionList <- NULL
    }
    if (plot.col.partition | plot.col.clusters){ ##??
      plot.col.partitionList <- get.sep(clusters=col.clusters,type="column")
    } else {
      plot.col.partitionList <- NULL
    }

    row.sepList <- sepList[[1]]
    if (!.invalid(row.sepList)){
      for (i in 1:length(row.sepList)){
        i.sep <- row.sepList[[i]]
        rect(
             xleft=i.sep[1]+0.5,
             ybottom=i.sep[2]+0.5,
             xright=i.sep[3]+0.5,
             ytop=i.sep[4]+0.5,
             lty=sep.lty,
             lwd=sep.lwd,
             col=FALSE,
             border=sep.color[1]
             )
      }
    }
    col.sepList <- sepList[[2]]
    if (!.invalid(col.sepList)){
      for (i in 1:length(col.sepList)){
        i.sep <- col.sepList[[i]]
        rect(
             xleft=i.sep[1]+0.5,
             ybottom=i.sep[2]+0.5,
             xright=i.sep[3]+0.5,
             ytop=i.sep[4]+0.5,
             lty=sep.lty,
             lwd=sep.lwd,
             col=FALSE,
             border=sep.color[2]
             )
      }
    }
    ## show traces
    min.scale <- min(breaks)
    max.scale <- max(breaks)

    x.scaled  <- .scale.x(t(x),min.scale,max.scale)

    if(.invalid(hline)) hline=median(breaks)
    if(.invalid(vline)) vline=median(breaks)

    if(trace %in% c("both","column")){
      for( i in colInd ){
        if(!.invalid(vline)){
          vline.vals <- .scale.x(vline,min.scale,max.scale)
          abline(v=i-0.5+vline.vals,col=linecol,lty=2)
        }
        xv <- rep(i,nrow(x.scaled))+x.scaled[,i]-0.5
        xv <- c(xv[1],xv)
        yv <- 1:length(xv)-0.5
        lines(x=xv,y=yv,lwd=1,col=tracecol,type="s")
      }
    }

    if(trace %in% c("both","row")){
      for( i in rowInd ){
        if(!.invalid(hline)){
          hline.vals <- .scale.x(hline,min.scale,max.scale)
          abline(h=i+hline,col=linecol,lty=2)
        }
        yv <- rep(i,ncol(x.scaled))+x.scaled[i,]-0.5
        yv <- rev(c(yv[1],yv))
        xv <- length(yv):1-0.5
        lines(x=xv,y=yv,lwd=1,col=tracecol,type="s")
      }
    }
    ## cellnote
    if(!.invalid(cellnote)){
      text(x=c(row(cellnote)),
           y=c(col(cellnote)),
           labels=c(cellnote),
           col=notecol,
           cex=cex.note)
    }

    ## 2) plot.row.partition
    if(plot.row.partition |plot.row.clusters) { ##row.clusters
      par(mar=c(margins[1],0.5,0,0.1))

      row.clusters.unique <- unique(row.clusters)
      row.clusters.unique <- row.clusters.unique[!is.na(row.clusters.unique)]

      image(rbind(1:nr),
            xlim=0.5+c(0,1),ylim=0.5+c(0,nr),
            col=par("bg"),
            axes=FALSE,add=force_add)
      if (!.invalid(plot.row.partitionList)){
        for (i in 1:length(plot.row.partitionList)){
          i.sep <- plot.row.partitionList[[i]]
          rect(
               xleft=0+0.5,
               ybottom=i.sep[2]+0.5,
               xright=1+0.5,
               ytop=i.sep[4]+0.5,
               lty=sep.lty,
               lwd=sep.lwd,
               col=color.partition.box,
               border=color.partition.border
               )
          g <- row.clusters.unique[i]
          s <- g
          text(x=1,y=(i.sep[2]+0.5+i.sep[4]+0.5)/2,labels=s,col=color.partition.border,
               cex=cex.partition,
               srt=90
               )
        }
      }

    }
    ## 3) plot.col.partition
    if(plot.col.partition | plot.col.clusters) {
      par(mar=c(0.1,0,0,margins[4]))
      col.clusters.unique <- unique(col.clusters)
      col.clusters.unique <- col.clusters.unique[!is.na(col.clusters.unique)]

      image(cbind(1:nc),
            xlim=0.5+c(0,nc),ylim=0.5+c(0,1),
            col=par("bg"),
            axes=FALSE,add=force_add)

      if (!.invalid(plot.col.partitionList)){
        for (i in 1:length(plot.col.partitionList)){
          i.sep <- plot.col.partitionList[[i]]
          rect(
               xleft=i.sep[1]+0.5,
               ybottom=0+0.5,
               xright=i.sep[3]+0.5,
               ytop=1+0.5,
               lty=sep.lty,
               lwd=sep.lwd,
               col=color.partition.box,
               border=color.partition.border
               )
          g <- col.clusters.unique[i]
          s <- g
          text(x=(i.sep[1]+0.5+i.sep[3]+0.5)/2,y=1,labels=s,col=color.partition.border,
               cex=cex.partition,
               srt=0
               )
        }
      }

    }
    ## 4) draw the side color bars - for row
    if(!.invalid(RowIndividualColors)) {
      par(mar=c(margins[1],0,0,0.5))
      # image(rbind(1:nr),col=RowIndividualColors[rowInd, 1],axes=FALSE,add=force_add)
      image(rbind(1:nr),col=RowIndividualColors[rowInd, 1],axes=FALSE,add=FALSE)

        if (dim(RowIndividualColors)[2] > 1) {
        par(mar=c(margins[1],0,0,0.5))
            image(rbind(1:nr),col=RowIndividualColors[rowInd, 2], axes=FALSE, add=force_add)
        }
    }

    ## 5) draw the side color bars - for col
    if(!.invalid(ColIndividualColors)) {
      par(mar=c(0.5,0,0,margins[4]))
      image(cbind(1:nc),col=ColIndividualColors[colInd],axes=FALSE,add=force_add)
    }

    ## 6) row-dend
    par(mar=c(margins[1],0,0,0))
    if(dendrogram %in% c("both","row")){
      plot(ddr,horiz=TRUE,axes=FALSE,yaxs="i",leaflab="none")
    }else{
      .plot.text(ylim=range(iy))
      if (sideRow==2){
        .sideRow <- par("usr")[2]-0.5*srtCol/90
        text(.sideRow,iy,labels=labRow,srt=srtRow,pos=1,xpd=TRUE,cex=cexRow)
      }
    }


    par(mar=c(0,0,0,0))
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
    legend(x=0.6,y=1.2 ,legend=annotations_legend[,1],
        cex=1.2,
        fill=annotations_legend[,2],
        ncol=grouping_key_coln)


    ## 7) col-dend and title
    mar3 <- (if(!is.null(main)) mar.main else 0) +
      (if(!is.null(sub)) mar.sub else 0)
    par(mar=c(0,0,mar3,margins[4]))

    if(dendrogram %in% c("both","column")){
      plot(ddc,axes=FALSE,xaxs="i",leaflab="none")
    } else{
      if(key) {
        .plot.text(xlim=range(1:nc))
      }
      if (sideCol==3){
        .sideCol <- par("usr")[3]+0.5*srtCol/90
        text(1:nc,.sideCol,labels=labCol,srt=srtCol,pos=1,xpd=TRUE,cex=cexCol)
      }
    }


    ## title
    if (is.null(sub)) main.line <- 1 else main.line <- 3
    if(!is.null(main)) title(main,cex.main=cex.main,adj=adj.main,mgp=mgp.main,font.main=font.main,line=main.line)
    if(!is.null(sub)) title(sub,cex.main=cex.sub,adj=adj.main,mgp=mgp.main,font.main=font.sub,line=0)
    ##if(!is.null(main)) title(main,cex.main=1.5*op[["cex.main"]])
    ## 8) plot the color-key
    if(key){
      cex.key <- 0.75
      op.ori <- par()

      par(mar=c(2,1.5,0.75,1)*keysize,cex=cex.key,mgp=c(0.75,0,0),tcl=-0.05)
      z <- seq(x.range[1],x.range[2],length=length(colors))
      flog.debug(paste("::inferCNV::heatmap.cnv colorkey z range: ", paste(z, collapse=","), sep=""))
      flog.debug(paste("::inferCNV::heatmap.cnv colorkey breaks range: ", paste(breaks, collapse=","), sep=""))
      flog.debug(paste("::inferCNV::heatmap.cnv colorkey colors range: ", paste(colors, collapse=","), sep=""))

      image(z=matrix(z,ncol=1),
            col=colors,
            breaks=breaks,
            xaxt="n",
            yaxt="n",
            xlab=key.xlab,
            ylab="",
            main="",add=force_add
            )

      par(usr=c(0,1,0,1))
      lv <- pretty(breaks)
      xv <- .scale.x(as.numeric(lv),x.range[1],x.range[2])
      axis(1,at=xv,labels=lv,cex.axis=cex.key*1)
      if(density.info=="density"){
        ## Experimental : also plot density of data
        dens <- density(x,adjust=densadj,na.rm=TRUE)
        omit <- dens$x < min(breaks) | dens$x > max(breaks)
        dens$x <- dens$x[-omit]
        dens$y <- dens$y[-omit]
        dens$x <- .scale.x(dens$x,x.range[1],x.range[2])
        lines(dens$x,dens$y / max(dens$y) * 0.95,col=denscol,lwd=1)
        axis(2,at=pretty(dens$y)/max(dens$y) * 0.95,pretty(dens$y),cex.axis=cex.key*1)
        ##title("Color Key and Density",cex.lab=cex.key*0.25)
        title(key.title,cex.main=cex.key,font.main=1)
        mtext(side=2,"Density",line=0.75,cex=cex.key)
      } else if(density.info=="histogram"){
        h <- hist(x,plot=FALSE,breaks=breaks)
        hx <- .scale.x(breaks,x.range[1],x.range[2])
        hy <- c(h$counts,h$counts[length(h$counts)])
        lines(hx,hy/max(hy)*0.95,lwd=1,type="s",col=denscol)
        axis(2,at=pretty(hy)/max(hy)*0.95,pretty(hy),cex.axis=cex.key*1)
        ##title("Color Key and Histogram",cex.main=cex.key*0.25)
        title(key.title,cex.main=cex.key,font.main=1)
        mtext(side=2,key.ylab,line=0.75,cex=cex.key)
      } else{
        title(key.title,cex.main=cex.key,font.main=1)
      }
    } else{
      if(!force_add){
      .plot.text()
      }
    }


    ## 9)
    if(plot.row.individuals) {
      .plot.text("Row\nIndividuals",cex=cex.text,bg="white")
      for (i in 1:ir) {
        ##.plot.text()
        tmp <- plot.row.individuals.list[[i]]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
      }
    }

    ## 10)
    if(plot.col.individuals) {
      .plot.text("Column\nIndividuals",cex=cex.text,bg="white",srt=90)
      for (i in 1:ic) {
        ##.plot.text()
        tmp <- plot.col.individuals.list[[i]]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
      }
    }

    ## 11) for RowPlot from bottom
    if (plot.row.clusters){
      .plot.text("Row\nClusters",cex=cex.text,bg="white")

      tmp <- plot.row.clusters.list[[1]]
      row.data <- row.data[rowInd]
      for (i in unique(row.clusters)){
        i.x <- row.data[row.clusters==i]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
        i.main <- sprintf("Row group %s (n=%s)",i,length(i.x))
        title(i.main,cex.main=1,font.main=1)
      }
    }
    ## 12) for ColPlot from left
    if (plot.col.clusters){
      .plot.text("Col\nClusters",cex=cex.text,bg="white",srt=90)

      tmp <- plot.col.clusters.list[[1]]
      col.data <- if(revC) col.data[rev(colInd)] else col.data[colInd]
      for (i in unique(col.clusters)){
        i.x <- col.data[col.clusters==i]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
        i.main <- sprintf("Col group %s (n=%s)",i,length(i.x))
        title(i.main,cex.main=1,font.main=1)
      }
    }

    ## 13)
    if(plot.row.clustering) {
      .plot.text("Row\nClustering",cex=cex.text,bg="white")

      for (i in 1:cr) {
        ##.plot.text()
        tmp <- plot.row.clustering.list[[i]]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
      }
    }
    ## 14)
    if(plot.col.clustering) {
      .plot.text("Column\nClustering",cex=cex.text,bg="white",srt=90)

      for (i in 1:cc) {
        ##.plot.text()
        tmp <- plot.col.clustering.list[[i]]
        for(j in 1:length(tmp)){
          eval(tmp[[j]])
        }
      }
    }

    ## 15) text
    if (!.invalid(text.box) & if.plot.info){
      .plot.text(text.box,cex=cex.text,bg="gray75")
    } else{
      if (flag.text){
        .plot.text()
      }
    }

  } else {
    x.save <- x
  }

  ret <-
    list(x.ori=x.ori,
         x=x.save,
         rowInd=rowInd,colInd=colInd,
         row.clusters=row.clusters,col.clusters=col.clusters,
         dist.row=dist.row,dist.col=dist.col,
         hclust.row=hclust.row,hclust.col=hclust.col,
         kr=kr,kc=kc
         )
  class(ret) <- c("hclustering",class(ret))
  invisible(ret)
}

## Please note this code is from the library GMD
## All credit for this code goes to GMD's authors.
## The official version from the package GMD is now archived
## https://cran.r-project.org/web/packages/GMD/index.html
gdist <-
  function(x,
           method="euclidean",
           MoreArgs=NULL,
           diag=FALSE,
           upper=FALSE
           )
{
  if(method %in% c("correlation","correlation.of.observations")){
    FUN <- function(x,...){
      as.dist(1-cor(t(x),y=NULL,...),diag=diag,upper=upper)}
    if (.invalid(MoreArgs)) MoreArgs=list(method="pearson",use="everything")
  } else if(method %in% c("correlation.of.variables")){
    FUN <- function(x,...){
      as.dist(1-cor(x,y=NULL,...),diag=diag,upper=upper)}
    if (.invalid(MoreArgs)) MoreArgs=list(method="pearson",use="everything")
  }

  COMMON_METHODS <-
    c("euclidean","maximum",
      "manhattan","canberra",
      "binary","minkowski"
      )

  if(method %in% COMMON_METHODS){
    d <- dist(x=x,method=method,diag=diag,upper=upper,p=MoreArgs$p)
  } else if (method %in% c("correlation","correlation.of.observations","correlation.of.variables")){
    ##d <- .call.FUN(FUN,x,MoreArgs)
    d <- FUN(x,method=MoreArgs$method,use=MoreArgs$use)
    attr(d,"method") <- method
  } else {
    FUN <- match.fun(method)
    MoreArgs[["diag"]] <- diag
    MoreArgs[["upper"]] <- upper
    d <- .call.FUN(FUN,x,MoreArgs)

    ## check attributes of the dist object ##
    if(is.null(attr(d,"method"))){
      attr(d,"method") <- method
    }
    if(is.null(attr(d,"call"))){
      attr(d,"call") <- match.call()
    }
    if(is.null(attr(d,"Size"))){
      warning(sprintf("The `dist' object returned by %s does not contain a specified attribute of `Size'.",attr(d,"method")))
    }
    if(is.null(attr(d,"Labels"))){
      warning(sprintf("The `dist' object returned by %s does not contain a specified attribute of `Labels'.",attr(d,"method")))
    }
  }
  attr(d,"Diag") <- diag
  attr(d,"Upper") <- upper
  class(d) <- "dist"
  return(d)
}


## Please note this code is from the library GMD
## All credit for this code goes to GMD's authors.
## I do not recommend using this version of the code, which
## has been poorly modified for our use but recommend using
## the official version from the package GMD
## https://cran.r-project.org/web/packages/GMD/index.html
## Get row or column lines of separation for \code{heatmap.3} according to clusters
## param clusters a numerical vector, indicating the cluster labels of observations.
## param type string, one of the following: \code{c("row","column","both")}
get.sep <-
  function(clusters,type=c("row","column","both"))
{
  ##   if(!is.numeric(clusters) | !(.is.grouped(clusters))){
  ## stop("`clusters' should be a grouped numeric vector.")
  ##   }
  tmp.whichna <- which(is.na(clusters))
  tmp.which <- which(!duplicated(clusters))

  tmp.sep <- data.frame(start=tmp.which,end=c(tmp.which[-1],length(clusters)+1)-1)
  tmp.sep2 <- tmp.sep[tmp.sep$start<=tmp.sep$end,]

  ## lines of separation
  sepList <- list()
  if (type=="row"){
    xmax <- length(clusters)
    for(i.s in 1:nrow(tmp.sep2)){
      sepList[[i.s]] <- c(0,tmp.sep2[i.s,1]-1,xmax,tmp.sep2[i.s,2])
    }
  } else if (type=="column"){
    ymax <- length(clusters)
    for(i.s in 1:nrow(tmp.sep2)){
      sepList[[i.s]] <- c(tmp.sep2[i.s,1]-1,0,tmp.sep2[i.s,2],ymax)
    }
  } else if (type=="both"){
    for(i.s in 1:nrow(tmp.sep2)){
      sepList[[i.s]] <- c(tmp.sep2[i.s,1]-1,tmp.sep2[i.s,1]-1,tmp.sep2[i.s,2],tmp.sep2[i.s,2])
    }
  }
  sepList
}


