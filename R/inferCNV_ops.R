
#' Function doing the actual analysis before calling the plotting functions.
#'
#' @title run() : Invokes a routine inferCNV analysis to Infer CNV changes given a matrix of RNASeq counts. 
#'
#' @param infercnv_obj An infercnv object populated with raw count data
#' 
#' @param cutoff Cut-off for the min average read counts per gene among reference cells. (default: 1)
#'
#' @param min_cells_per_gene minimum number of reference cells requiring expression measurements to include the corresponding gene. 
#'                           default: 3     
#' 
#' @param out_dir path to directory to deposit outputs (default: '.')
#'
#' @param normalize_factor scaling factor for total sum of counts (default: NA, in which case
#'                         will be set = 10^round(log10(mean(colSums))), typically setting to 1e5
#'
#' @param window_length Length of the window for the moving average
#'                          (smoothing). Should be an odd integer. (default: 101)#'
#'
#' @param num_ref_groups The number of reference groups or a list of
#'                           indices for each group of reference indices in
#'                           relation to reference_obs. (default: NULL)
#'
#' @param max_centered_threshold The maximum value a value can have after
#'                                   centering. Also sets a lower bound of
#'                                   -1 * this value.
#'
#' @param noise_filter  Values +- from the reference cell mean will be set to zero (whitening effect)
#'                      default(NA, instead will use sd_amplifier below.
#'
#' @param sd_amplifier  Noise is defined as mean(reference_cells) +- sdev(reference_cells) * sd_amplifier
#'                      default: 1.0
#' 
#' @param cluster_by_groups   If observations are defined according to groups (ie. patients), each group
#'                            of cells will be clustered separately. (default=FALSE, instead will use k_obs_groups setting)
#' 
#'
#' @param k_obs_groups Number of groups in which to break the observations. (default: 1)
#'
#' @param outlier_method_bound Method to use for bounding outlier values. (default: "average_bound")
#'                             Will preferentially use outlier_lower_bounda and outlier_upper_bound if set.
#' @param outlier_lower_bound  Outliers below this lower bound will be set to this value. 
#' @param outlier_upper_bound  Outliers above this upper bound will be set to this value.
#'
#'
#' @param hclust_method Method used for hierarchical clustering of cells. Valid choices are:
#' "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid".
#'
#' @param anscombe_normalize  Perform anscombe normalization on normalized counts before log transformation.
#' 
#' @param use_zscores If true, converts log(expression) data to zscores based on reference cell expr distribution.
#' 
#' @param remove_genes_at_chr_ends If true, removes the window_length/2 genes at both ends of the chromosome.
#'
#' @param mask_nonDE_genes If true, sets genes not significantly differentially expressed between tumor/normal to
#'                          the mean value for the complete data set
#'
#' @param mask_nonDE_pval  p-value threshold for defining statistically significant DE genes between tumor/normal
#
#' @param test.use statistical test to use.  (default: "wilcoxon") alternatives include 'perm' or 't'.'
#' 
#' @return infercnv_obj containing filtered and transformed data
#'
#' @param plot_steps If true, saves infercnv objects and plots data at the intermediate steps.
#'
#' @param include.spike  If true, introduces an artificial spike-in of data at ~0x and 2x for scaling residuals between 0-2. (default: F)
#'
#' @param spike_in_chrs  vector listing of chr names to use for modeling spike-ins (default: NULL - uses the two largest chrs.  ex.  c('chr1', 'chr2') )
#'
#' @param spike_in_multiplier vector of weights matching spike_in_chrs (default: c(0.01, 2.0) for modeling loss/gain of both chrs)
#'
#' @param pseudocount  Number of counts to add to each gene of each cell post-filtering of genes and cells and pre-total sum count normalization. (default: 0)
#'
#' @param debug If true, output debug level logging.
#' 
#' @export
#'

run <- function(infercnv_obj,

                # gene filtering settings
                cutoff=1,
                min_cells_per_gene=3,
                
                out_dir=".",
                normalize_factor=NA,
                window_length=101,
                
                num_ref_groups=NULL,

                max_centered_threshold=NA,

                # noise settings
                noise_filter=NA,
                sd_amplifier = 1.5,

                # observation cell clustering settings
                cluster_by_groups=FALSE,
                k_obs_groups=1,


                # outlier adjustment settings
                outlier_method_bound="average_bound",
                outlier_lower_bound=NA,
                outlier_upper_bound=NA,

                hclust_method='complete',

                anscombe_normalize=TRUE,
                use_zscores=FALSE,
                remove_genes_at_chr_ends=FALSE,

                mask_nonDE_genes=FALSE,
                mask_nonDE_pval=0.05,
                test.use='wilcoxon',
                
                plot_steps=FALSE,

                debug=FALSE, #for debug level logging

                include.spike = FALSE,

                #  must specify both below if to be used, and must match in vec length
                spike_in_chrs =  NULL, # use defaults
                spike_in_multiplier_vec = NULL, # use defaults
                
                pseudocount = 0
                
                ) {


    if (debug) {
        flog.threshold(DEBUG)
    }
    
    flog.info(paste("::process_data:Start", sep=""))

    if(out_dir != "." & !file.exists(out_dir)){
        dir.create(out_dir)
    }

    step_count = 0; 

    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %d: incoming data\n", step_count))
    
    # Plot incremental steps.
    if (plot_steps) {        

        infercnv_obj_incoming_data <- infercnv_obj
        save('infercnv_obj_incoming_data', file=file.path(out_dir, sprintf("%02d_incoming_data.infercnv_obj", step_count)))
        
    }

           
    ###################################################
    ## Step: removing insufficiently expressed genes
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: Removing lowly expressed genes\n", step_count))
    
    # Remove genes that aren't sufficiently expressed, according to min mean count cutoff.
    # Examines the original (non-log-transformed) data, gets mean for each gene, and removes genes
    #  with mean values below cutoff.

    infercnv_obj <- require_above_min_mean_expr_cutoff(infercnv_obj, cutoff)
    
    ## require each gene to be present in a min number of cells for ref sets

    infercnv_obj <- require_above_min_cells_ref(infercnv_obj, min_cells_per_gene=min_cells_per_gene)
    

    if (plot_steps){
        
        infercnv_obj_low_expr_genes_pruned <- infercnv_obj
        
        save('infercnv_obj_low_expr_genes_pruned', file=file.path(out_dir, sprintf("%02d_reduced_by_cutoff.infercnv_obj",step_count)))
        
    }


    ####################################################
    ##  Adding pseudocounts
    
    if (pseudocount != 0) {
        flog.info(sprintf("Adding pseudocount: %g", pseudocount))
        infercnv_obj <- add_pseudocount(infercnv_obj, pseudocount)
    }
        
    
    ###########################################
    ### STEP: normalization by sequencing depth
    
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: normalization by sequencing depth\n", step_count))
    
    infercnv_obj <- normalize_counts_by_seq_depth(infercnv_obj, normalize_factor=normalize_factor)

    if (plot_steps){
        
        infercnv_obj_normalize_by_depth <- infercnv_obj
        save('infercnv_obj_normalize_by_depth', file=file.path(out_dir, sprintf("%02d_normalized_by_depth.infercnv_obj", step_count)))
        
    }
    

    ##################################################
    ## spike-in
    
    if (include.spike) {
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: Spiking in genes with variation added for tracking\n", step_count))

        if (! (is.null(spike_in_chrs) && is.null(spike_in_multiplier_vec)) ) {
            infercnv_obj <- spike_in_variation_chrs(infercnv_obj, spike_in_chrs, spike_in_multiplier_vec)
        } else {
            infercnv_obj <- spike_in_variation_chrs(infercnv_obj)
        }
        
        # Plot incremental steps.
        if (plot_steps){

            infercnv_obj_spiked <- infercnv_obj
            save('infercnv_obj_spiked', file=file.path(out_dir, sprintf("%02d_spiked.infercnv_obj", step_count)))
            
            plot_cnv(infercnv_obj=infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     out_dir=out_dir,
                     color_safe_pal=FALSE,
                     x.center=mean(infercnv_obj@expr.data),
                     x.range="auto",
                     title=sprintf("%02d_spike_added",step_count),
                     obs_title="Observations (Cells)",
                     ref_title="References (Cells)",
                     output_filename=sprintf("infercnv.%02d_spike_added",step_count),
                     write_expr_matrix=TRUE
                     )
        }
        
    }

    
    ##################################
    ##### STEP: anscombe normalization

    if (anscombe_normalize) {
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: anscombe normalization\n", step_count))    

        infercnv_obj <- anscombe_transform(infercnv_obj)

        if (plot_steps) {
            infercnv_obj_anscombe_norm <- infercnv_obj
            save('infercnv_obj_anscombe_norm', file=file.path(out_dir, sprintf("%02d_anscombe_normalization.infercnv_obj", step_count)))

            plot_cnv(infercnv_obj=infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     out_dir=out_dir,
                     color_safe_pal=FALSE,
                     x.center=mean(infercnv_obj@expr.data),
                     x.range="auto",
                     title=sprintf("%02d_anscombe_norm",step_count),
                     obs_title="Observations (Cells)",
                     ref_title="References (Cells)",
                     output_filename=sprintf("infercnv.%02d_anscombe_norm",step_count),
                     write_expr_matrix=TRUE
                     )
            
        }
        
    }

    ###########################
    ## Step: log transformation
    
    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: log transformation of data\n", step_count))    
    
    infercnv_obj <- log2xplus1(infercnv_obj)
    
    # Plot incremental steps.
    if (plot_steps){
        
        infercnv_obj_log_transformed <- infercnv_obj
        save('infercnv_obj_log_transformed', file=file.path(out_dir, sprintf("%02d_logtransformed.infercnv_obj", step_count)))
        
        plot_cnv(infercnv_obj=infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_dir,
                 color_safe_pal=FALSE,
                 x.center=mean(infercnv_obj@expr.data),
                 x.range="auto",
                 title=sprintf("%02d_log_transformed_data",step_count),
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename=sprintf("infercnv.%02d_log_transformed",step_count),
                 write_expr_matrix=TRUE
                 )
    }
    
    

    ###############################
    ### STEP: ZScore transformation

    if (use_zscores) {
        
        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: Z-score transformation of data\n", step_count))    
        
        infercnv_obj <- transform_to_reference_based_Zscores(infercnv_obj)
        
        
        if (plot_steps){
            
            infercnv_obj_zscores <- infercnv_obj
            
            save('infercnv_obj_zscores', file=file.path(out_dir, sprintf("%02d_Z-scores.infercnv_obj", step_count)))
            
            plot_cnv(infercnv_obj=infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     out_dir=out_dir,
                     color_safe_pal=FALSE,
                     x.center=0,
                     x.range="auto",
                     title=sprintf("%02d_centering_gene_expr",step_count),
                     obs_title="Observations (Cells)",
                     ref_title="References (Cells)",
                     output_filename=sprintf("infercnv.%02d_centering_gene_expr",step_count),
                     write_expr_matrix=TRUE)
            
        }
    }
        
        
    #######################################################
    ## Apply maximum centered expression thresholds to data
    # Cap values between threshold and -threshold, retaining earlier center

    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: apply max centered expression threshold\n", step_count))
    
    threshold = max_centered_threshold
    if (is.na(max_centered_threshold)) {
        threshold = mean(abs(get_average_bounds(infercnv_obj)))
    }
    
    infercnv_obj <- apply_max_threshold_bounds(infercnv_obj, threshold=threshold)
    
    # Plot incremental steps.
    if (plot_steps){

        infercnv_obj_max_centered_expr <- infercnv_obj
        
        save('infercnv_obj_max_centered_expr', file=file.path(out_dir, sprintf("%02d_apply_max_centered_expr_threshold.infercnv_obj", step_count)))

        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_dir,
                 color_safe_pal=FALSE,
                 x.center=mean(infercnv_obj@expr.data),
                 x.range="auto",
                 title=sprintf("%02d_apply_max_centered_expr_threshold",step_count),
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename=sprintf("infercnv.%02d_apply_max_centred_expr_threshold",step_count),
                 write_expr_matrix=TRUE)
        
    }
    
    
    ###########################################################################
    # Step: For each cell, smooth the data along chromosome with gene windows

    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: Smoothing data per cell by chromosome\n", step_count))
    
    infercnv_obj <- smooth_by_chromosome(infercnv_obj, window_length=window_length, smooth_ends=TRUE)
    
    
    # Plot incremental steps.
    if (plot_steps){

        infercnv_obj_smoothed_by_chr <- infercnv_obj
        save('infercnv_obj_smoothed_by_chr', file=file.path(out_dir, sprintf("%02d_smoothed_by_chr.infercnv_obj", step_count)))
        
        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_dir,
                 color_safe_pal=FALSE,
                 x.center=mean(infercnv_obj@expr.data),
                 x.range="auto",
                 title=sprintf("%02d_smoothed_by_chr",step_count),
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename=sprintf("infercnv.%02d_smoothed_by_chr", step_count))
    }
    
    

    ## 
    # Step: 
    # Center cells/observations after smoothing. This helps reduce the
    # effect of complexity.

    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: re-centering data across chromosome after smoothing\n", step_count))
    
    infercnv_obj <- center_cell_expr_across_chromosome(infercnv_obj, method="median")
    
    
    # Plot incremental steps.
    if (plot_steps) {

        infercnv_obj_cell_centered <- infercnv_obj
        
        save('infercnv_obj_cell_centered', file=file.path(out_dir, sprintf("%02d_recentered_cells_by_chr.infercnv_obj", step_count)))
        
        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_dir,
                 color_safe_pal=FALSE,
                 x.center=mean(infercnv_obj@expr.data),
                 x.range="auto",
                 title=sprintf("%02d_centering_of_smoothed",step_count),
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename=sprintf("infercnv.%02d_centering_of_smoothed", step_count))
        
    }
    
    ###################################################
    # Step: Split the reference data into groups if requested
    
    if (!is.null(num_ref_groups)) {

        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: splitting reference data into %d clusters\n", step_count, num_ref_groups))
        
        infercnv_obj <- split_references(infercnv_obj,
                                         num_groups=num_ref_groups,
                                         hclust_method=hclust_method)
        
        
        
        
    }
    
    
    ####################################
    ## Step: Subtract average reference
    ## Since we're in log space, this now becomes log(fold_change)

    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: removing average of reference data\n", step_count))
    
    infercnv_obj <- subtract_ref_expr_from_obs(infercnv_obj, inv_log=TRUE)
    
    
    # Plot incremental steps.
    if (plot_steps){
                
        infercnv_obj_subtract_ref <- infercnv_obj
        
        save('infercnv_obj_subtract_ref', file=file.path(out_dir, sprintf("%02d_remove_ref_avg_from_obs.infercnv_obj", step_count)))

        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_dir,
                 color_safe_pal=FALSE,
                 x.center=0,
                 x.range="auto",
                 title=sprintf("%02d_remove_average",step_count),
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename=sprintf("infercnv.%02d_remove_average", step_count))
        
    }

    
    ## Step 08:
    # Remove Ends

    if (remove_genes_at_chr_ends == TRUE) {

        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: removing genes at chr ends\n", step_count))
        
        infercnv_obj <- remove_genes_at_ends_of_chromosomes(infercnv_obj, window_length)
        
        # Plot incremental steps.
        if (plot_steps){
            
            infercnv_obj_remove_chr_end_genes <- infercnv_obj
            
            save('infercnv_obj_remove_chr_end_genes', file=file.path(out_dir, sprintf("%02d_remove_gene_at_chr_ends.infercnv_obj", step_count)))
            
            plot_cnv(infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     out_dir=out_dir,
                     color_safe_pal=FALSE,
                     x.center=0,
                     x.range="auto",
                     title=sprintf("%02d_remove_genes_at_chr_ends",step_count),
                     obs_title="Observations (Cells)",
                     ref_title="References (Cells)",
                     output_filename=sprintf("infercnv.%02d_remove_genes_at_chr_ends",step_count),
                     write_expr_matrix=TRUE)
            
        }
    }
    
    
    #############################
    # Step: invert log transform  (convert from log(FC) to FC)

    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: invert log2(FC) to FC\n", step_count))
    
    infercnv_obj <- invert_log2(infercnv_obj)

    # Plot incremental steps.
    if (plot_steps) {
        
        infercnv_obj_invert_log_transform <- infercnv_obj
        
        save('infercnv_obj_invert_log_transform', file=file.path(out_dir, sprintf("%02d_invert_log_transform.infercnv_obj", step_count)))
        
        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_dir,
                 color_safe_pal=FALSE,
                 x.center=1,
                 x.range="auto",
                 title=sprintf("%02d_invert_log_transform log(FC)->FC",step_count),
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename=sprintf("infercnv.%02d_invert_log_FC",step_count),
                 write_expr_matrix=TRUE)
        
    }
    
    
    ################################
    # Step: de-noising 

    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: Denoising\n", step_count))
    
    if (! is.na(noise_filter)) {

        if (noise_filter > 0) {
            flog.info(paste("::process_data:Remove noise, noise threshold at: ", noise_filter))
            infercnv_obj <- clear_noise(infercnv_obj,
                                        threshold=noise_filter)
        }
        else {
                                        # noise == 0 or negative...
                                        # don't remove noise.
        }
        
    }
    else {
        # default, use quantiles, if NA 
        flog.info(paste("::process_data:Remove noise, noise threshold defined via ref mean sd_amplifier: ", sd_amplifier))
        infercnv_obj <- clear_noise_via_ref_mean_sd(infercnv_obj,
                                                    sd_amplifier = sd_amplifier)
    }



    
    if (plot_steps){
        
        infercnv_obj_denoised <- infercnv_obj
        
        save('infercnv_obj_denoised', file=file.path(out_dir, sprintf("%02d_denoise.infercnv_obj", step_count)))
        
        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_dir,
                 color_safe_pal=FALSE,
                 x.center=1,
                 x.range="auto",
                 title=sprintf("%02d_denoised", step_count),
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename=sprintf("infercnv.%02d_denoised", step_count))
        
    }
    
    ##################################
    # STEP: Remove outliers for viz

    step_count = step_count + 1
    flog.info(sprintf("\n\n\tSTEP %02d: Removing outliers\n", step_count))

    
    infercnv_obj = remove_outliers_norm(infercnv_obj,
                                        out_method=outlier_method_bound,
                                        lower_bound=outlier_lower_bound,
                                        upper_bound=outlier_upper_bound)
    
    
    # Plot incremental steps.
    if (plot_steps) {

        infercnv_obj_remove_outliers <- infercnv_obj

        save('infercnv_obj_remove_outliers', file=file.path(out_dir, sprintf("%02d_remove_outlier.infercnv_obj", step_count)))
        
        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 out_dir=out_dir,
                 color_safe_pal=FALSE,
                 x.center=1,
                 x.range="auto",
                 title=sprintf("%02d_removed_outliers",step_count),
                 obs_title="Observations (Cells)",
                 ref_title="References (Cells)",
                 output_filename=sprintf("infercnv.%02d_removed_outliers", step_count))
    }



    # define heatmap thresholds for final plots.
    plot_data = infercnv_obj@expr.data
    high_threshold = max(abs(quantile(plot_data[plot_data != 0], c(0.05, 0.95))))
    low_threshold = -1 * high_threshold
    
    if (include.spike) {

        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: Scaling according to spike\n", step_count))
                
        # normalize by spike
        infercnv_obj <- scale_cnv_by_spike(infercnv_obj)

        # now thresholds should be between 0 and 2 after spike-based scaling
        low_threshold = 0
        high_threshold = 2
        
        if (plot_steps) {
            infercnv_obj_scaled_by_spike <- infercnv_obj
            save('infercnv_obj_scaled_by_spike', file=file.path(out_dir, sprintf("%02d_scaled_by_spike.infercnv_obj", step_count)))
            
            plot_cnv(infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     out_dir=out_dir,
                     color_safe_pal=FALSE,
                     x.center=1,
                     x.range=c(low_threshold, high_threshold),
                     title=sprintf("%02d_scaled_by_spike",step_count),
                     obs_title="Observations (Cells)",
                     ref_title="References (Cells)",
                     output_filename=sprintf("infercnv.%02d_scaled_by_spike", step_count))
        }
        
    }
    
    
    
    ## Step: Filtering significantly DE genes
    if (mask_nonDE_genes) {

        step_count = step_count + 1
        flog.info(sprintf("\n\n\tSTEP %02d: Identify and mask non-DE genes\n", step_count))
        

        infercnv_obj <- mask_non_DE_genes_basic(infercnv_obj, test.use = test.use, center_val=mean(plot_data))


        # Plot incremental steps.
        if (plot_steps) {
            
            infercnv_obj_mask_nonDE <- infercnv_obj
            
            save('infercnv_obj_mask_nonDE', file=file.path(out_dir, sprintf("%02d_mask_nonDE.infercnv_obj", step_count)))
            
            plot_cnv(infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     out_dir=out_dir,
                     color_safe_pal=FALSE,
                     x.center=1,
                     x.range=c(low_threshold,high_threshold),
                     title=sprintf("%02d_mask_nonDE",step_count),
                     obs_title="Observations (Cells)",
                     ref_title="References (Cells)",
                     output_filename=sprintf("infercnv.%02d_mask_nonDE", step_count))
            
            
        }
    }

    if  (include.spike) {
        # remove the spike before making the final plot.
        infercnv_obj <- remove_spike(infercnv_obj)
    }
    
    
    save('infercnv_obj', file=file.path(out_dir, "run.final.infercnv_obj"))
    
    flog.info("Making the final infercnv heatmap")
    plot_cnv(infercnv_obj,
             k_obs_groups=k_obs_groups,
             cluster_by_groups=cluster_by_groups,
             out_dir=out_dir,
             color_safe_pal=FALSE,
             x.center=1,
             x.range=c(low_threshold,high_threshold),
             title="inferCNV",
             obs_title="Observations (Cells)",
             ref_title="References (Cells)",
             output_filename="infercnv")
    
        
    return(infercnv_obj)

}

#' Function for Generating a next-generation heatmap
#'
#' @title ngchm() : generates next gen heatmap
#'
#' @param infercnv_obj An infercnv object
#'
#' @param out_dir  output directory (default: '.')
#'
#' @param title title of the interactive heatmap (default: "NGCHM")
#'
#' @param gene_symbol ##TODO  (default: NULL)
#'
#' @param path_to_shaidyMapGen path to the shaidyMapGen jar file (default: NULL)
#' 
#' @param x.center (integer) Center expression value for heatmap coloring.
#' 
#' @param x.range (integer) Values for minimum and maximum thresholds for heatmap coloring. 
#' 
#' @export
#'

ngchm <- function(infercnv_obj, 
                       out_dir=".", 
                       title="NGCHM", 
                       gene_symbol=NULL, 
                       path_to_shaidyMapGen=NULL,
                       x.range = NA,
                       x.center = NA) {
    
    if (!is.null(path_to_shaidyMapGen)) {
        shaidy.path <- unlist(strsplit(path_to_shaidyMapGen, split = .Platform$file.sep))
        if (!file.exists(path_to_shaidyMapGen) || tail(shaidy.path, n = 1L) != "ShaidyMapGen.jar"){
            error_message <- paste("Cannot find the file ShaidyMapGen.jar using the parameter \"path_to_shaidyMapGen\".",
                                   "Check that the correct pathway is being used.")
            flog.error(error_message)
            stop(error_message)
        }
    } else {
        path_to_shaidyMapGen <- Sys.getenv("SHAIDYMAPGEN")
        if (!file.exists(path_to_shaidyMapGen)){ ## check if envionrmental variable is passed
            error_message <- paste("Cannot find the file ShaidyMapGen.jar using SHAIDYMAPGEN.",
                                   "Check that the correct pathway is being used.")
            flog.error(error_message)
            stop(error_message)
        }
    }
    
    if (!requireNamespace("NGCHM", quietly=TRUE)) {
        stop("The \"NGCHM\" library is required to use \"-ngchm=TRUE\" but it is not available.", .call=FALSE)
    }
    
    flog.info("Creating NGCHM as infercnv.ngchm")
    Create_NGCHM(infercnv_obj = infercnv_obj,
                 path_to_shaidyMapGen = path_to_shaidyMapGen,
                 out_dir = out_dir,
                 title = title,
                 gene_symbol = gene_symbol,
                 x.range = x.range,
                 x.center = x.center)
}



#' Subtracting the mean of the reference expr distributions from the observed cells.
#'
#' @title subtract_ref_expr_from_obs()
#'
#' @description Remove the average of the genes of the reference observations from all
#' observations' expression. In the case there are multiple reference groupings,
#' the averages are computed separately for each reference grouping, and the min|max
#' of the averages are subtracted from the observation expression levels. Any values within the range
#' of the min,max of the group are set to zero.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param inv_log mean values will be determined based on (2^x -1)
#'
#' @return infercnv_obj containing the reference subtracted values.
#'
#' @export
#'

subtract_ref_expr_from_obs <- function(infercnv_obj, inv_log=FALSE) {
                                            # r = genes, c = cells
    flog.info(paste("::subtract_ref_expr_from_obs:Start", sep=""))
            
    ref_groups = infercnv_obj@reference_grouped_cell_indices
    
    subtr_data <- .subtract_expr(infercnv_obj@expr.data, ref_groups, inv_log)
    
    colnames(subtr_data) <- colnames(infercnv_obj@expr.data)

    infercnv_obj@expr.data <- subtr_data
    
    flog.info("subtracting mean(normal) per gene per cell across all data")
                
    return(infercnv_obj)

}

.subtract_expr <- function(expr_matrix, ref_groups, inv_log=FALSE) {
    
    subtract_normal_expr_fun <- function(x) {
                
        grp_min = NA
        grp_max = NA
        
        for (ref_group in ref_groups) {
            
            if (inv_log) {
                ref_grp_mean = log2(mean(2^x[ref_group] - 1) + 1)
            } else {
                ref_grp_mean = mean(x[ref_group])
            }
            
            grp_min = min(grp_min, ref_grp_mean, na.rm=T)
            grp_max = max(grp_max, ref_grp_mean, na.rm=T)
            
        }
        
        row_init = rep(0, length(x))
        
        above_max = which(x>grp_max)
        below_min = which(x<grp_min)
        
        row_init[above_max] <- x[above_max] - grp_max
        row_init[below_min] <- x[below_min] - grp_min
        
        return(row_init)
    }
    
    subtr_data <- t(apply(expr_matrix, 1, subtract_normal_expr_fun))

    return(subtr_data)

}

    
# Not testing, parameters ok.
# Helper function allowing greater control over the steps in a color palette.
# Source:http://menugget.blogspot.com/2011/11/define-color-steps-for-
#               colorramppalette.html#more

# Args:
# steps Vector of colors to change use in the palette
# between: Steps where gradients change
#
# Returns:
# Color palette
color.palette <- function(steps,
                          between=NULL, ...){

    if (is.null(between)){
        between <- rep(0, (length(steps) - 1))
    }
    if (length(between) != length(steps) - 1){
        stop("Must have one less \"between\" value than steps")
    }

    fill.steps <- cumsum(rep(1, length(steps)) + c(0, between))
    RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
    RGB[, fill.steps] <- col2rgb(steps)
    for(i in which(between > 0)){
        col.start <- RGB[, fill.steps[i]]
        col.end <- RGB[, fill.steps[i + 1]]
        for(j in seq(3)){
            vals <- seq(col.start[j],
                        col.end[j],
                        length.out=between[i] + 2)
            vals <- vals[2:(2 + between[i] - 1)]
            RGB[j, (fill.steps[i] + 1):(fill.steps[i + 1] - 1)] <- vals
        }
    }
    new.steps <- rgb(RGB[1, ], RGB[2, ], RGB[3, ], maxColorValue=255)
    pal <- colorRampPalette(new.steps, ...)
    return(pal)
}

# Create a sepList forthe heatmap.3 plotting function given integer vectors
# of rows and columns where speration should take place.
# The expected input to the heatmap function is a list of 2 lists.
# The first list are column based rectangles, and the second row.
# To define a rectagle the index of the row or column where the line of the rectagle
# should be placed is done with a vector of integers, left, bottom, right and top line.
# Ie. list(list(c(1,0,3,10), c(5, 0, 10,10)), list(c(1,2,3,4)))
#
# Args:
# row_count Total number of rows
# col_count Total number of columns
# row_seps Vector of integers indices for row breaks
# col_seps Vector of integer indices for column breaks
#
# Returns
# List of lists of vectors
create_sep_list <- function(row_count,
                            col_count,
                            row_seps=NULL,
                            col_seps=NULL){
    sepList <- list()
    # Heatmap.3 wants a list of boxes for seperating columns
    # Column data
    if(!is.null(col_seps) &&
       !is.na(col_seps) &&
       (length(col_seps)>0) &&
       col_count > 0){
        colList <- list()
        for(sep in 1:length(col_seps)){
            colList[[sep]] <- c(col_seps[sep],0,col_seps[sep],row_count)
        }
        sepList[[1]] <- colList
    } else {
        sepList[[1]] <- list()
        sepList[[1]][[1]] <- c(0,0,0,0)
    }

    # Row data
    # This is measured from bottom to top
    # So you have to adjust the values of the data
    row_seps <- row_count-row_seps
    if(!is.null(row_seps) &&
       !is.na(row_seps) &&
       (length(row_seps)>0) &&
       row_count > 0){
        rowList <- list()
        for(sep in 1:length(row_seps)){
            rowList[[sep]] <- c(0,row_seps[sep],col_count,row_seps[sep])
        }
        sepList[[2]] <- rowList
    } else {
        sepList[[2]] <- list()
        sepList[[2]][[1]] <- c(0,0,0,0)
    }
    return(sepList)
}


#' @title split_references()
#' 
#' @description Split up reference observations in to k groups based on hierarchical clustering.
#'
#' 
#' @param infercnv_obj infercnv_object
#'
#' @param num_groups (default: 2)
#'
#' @param hclust_method clustering method to use (default: 'complete')
#'
#' @return infercnv_obj
#' 
#' @export
#'

split_references <- function(infercnv_obj,
                             num_groups=2,
                             hclust_method='complete') {

        
    flog.info(paste("::split_references:Start", sep=""))
    
    ref_expr_matrix = infercnv_obj@expr.data[ , get_reference_grouped_cell_indices(infercnv_obj) ]
    
    hc <- hclust(dist(t(ref_expr_matrix)), method=hclust_method)
    
    split_groups <- cutree(hc, k=num_groups)
        
    ref_groups <- list()
        
    grp_counter = 0
    for (cut_group in unique(split_groups)) {
        grp_counter = grp_counter + 1
        grp_name = sprintf("refgrp-%d", grp_counter)
        cell_names <- names(split_groups[split_groups == cut_group])
        ref_groups[[grp_name]] <- which(colnames(ref_expr_matrix) %in% cell_names)
    }
    
    infercnv_obj@reference_grouped_cell_indices <- ref_groups
    
    return(infercnv_obj)
}


#' @title remove_outliers_norm()
#' 
#' @description Set outliers to some upper or lower bound.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param out_method method for computing the outlier bounds (default: "average_bound", involving
#'                   determining the range of values for each cell, and then taking the mean of those bounds.)
#'
#' @param lower_bound setting the lower bound for the data (default: NA, uses out_method above)
#'
#' @param upper_bound setting the upper bound for the data (default: NA, uses out_method above)
#'
#' @return infercnv_obj with data bounds set accordingly.
#' 
#' @export
#'

remove_outliers_norm <- function(infercnv_obj,
                                 out_method="average_bound",
                                 lower_bound=NA,
                                 upper_bound=NA) {

    flog.info(paste("::remove_outlier_norm:Start",
                    "out_method:", out_method,
                    "lower_bound:" , lower_bound,
                    "upper_bound:", upper_bound))


    infercnv_obj@expr.data <- .remove_outliers_norm(data=infercnv_obj@expr.data,
                                                    out_method=out_method,
                                                    lower_bound=lower_bound,
                                                    upper_bound=upper_bound)
    
    return(infercnv_obj)
    
}


.remove_outliers_norm <- function(data,
                                  out_method="average_bound",
                                  lower_bound=NA,
                                  upper_bound=NA) {
    
    flog.info(paste("::remove_outlier_norm:Start",
                    "out_method:", out_method,
                    "lower_bound:" , lower_bound,
                    "upper_bound:", upper_bound))

    

    if(is.null(data) || nrow(data) < 1 || ncol(data) < 1){
        flog.error("::remove_outlier_norm: Error, something is wrong with the data, either null or no rows or columns")
        stop("Error, something is wrong with the data, either null or no rows or columns")
    }
    
    
    
    if ( (! is.na(lower_bound)) & (! is.na(upper_bound)) ) {
        
        ## using the values as specified.
        flog.info(paste("::remove_outlier_norm: using hard thresholds: ",
                        "lower_bound:" , lower_bound,
                        "upper_bound:", upper_bound) )
        
    } else if (! is.na(out_method)) {
        
        # using out_method instead of specified bounds.
        flog.info(paste("::remove_outlier_norm using method:", out_method, "for defining outliers."))
                
        if (out_method == "average_bound"){

            bounds = .get_average_bounds(data)
            lower_bound = bounds[1]
            upper_bound = bounds[2]

            flog.info(sprintf("outlier bounds defined between: %g - %g", lower_bound, upper_bound))
            
        } else {
            flog.error(paste("::remove_outlier_norm:Error, please",
                                    "provide an approved method for outlier",
                                    "removal for visualization."))
            stop(991)
        }
    } else {
        flog.error("::remove_outlier_norm:Error, must specify outmethod or define exact bounds")
        stop(992)
    }
    
    # apply bounds
    data[data < lower_bound] <- lower_bound
    data[data > upper_bound] <- upper_bound


    return(data)
}




#' @title center_cell_expr_across_chromosome()
#'
#' @description Centers expression data across all genes for each cell, using the cell mean or median expression
#' value as the center.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param method method to select the center of the cell expression value. (default: 'mean', options: 'mean,'median')
#'
#' @return infercnv_object
#' 
#' @export
#'

center_cell_expr_across_chromosome <- function(infercnv_obj, method="mean") { # or median

    flog.info(paste("::center_smooth across chromosomes per cell"))
    
    infercnv_obj@expr.data <- .center_columns(infercnv_obj@expr.data, method)
    
    return(infercnv_obj)
}

.center_columns <- function(expr_data, method) {

    # Center within columns (cells)
    if (method == "median") {
        row_median <- apply(expr_data, 2, function(x) { median(x, na.rm=T) } )
        
        expr_data <- t(apply(expr_data, 1, "-", row_median))
    }
    else {
        # by mean
        row_means <- apply(expr_data, 2, function(x) { mean(x, na.rm=T) } )
        
        expr_data <- t(apply(expr_data, 1, "-", row_means))
    }
    return(expr_data)
}




#' @title require_above_min_mean_expr_cutoff ()
#'
#' @description Filters out genes that have fewer than the corresponding mean value across all cell values.
#'
#' @param infercnv_obj  infercnv_object
#'
#' @param min_mean_expr_cutoff  the minimum mean value allowed for a gene to be retained in the expression matrix.
#'
#' @return infercnv_obj  the infercnv_object with lowly or unexpressed genes removed.
#' 
#' @export
#'

require_above_min_mean_expr_cutoff <- function(infercnv_obj, min_mean_expr_cutoff) {

    flog.info(paste("::above_min_mean_expr_cutoff:Start", sep=""))


    indices <-.below_min_mean_expr_cutoff(infercnv_obj@expr.data, min_mean_expr_cutoff)
    if (length(indices) > 0) {
        flog.info(sprintf("Removing %d genes from matrix as below mean expr threshold: %g",
                          length(indices), min_mean_expr_cutoff))
        
        infercnv_obj <- remove_genes(infercnv_obj, indices)

        expr_dim = dim(infercnv_obj@expr.data)
        flog.info(sprintf("There are %d genes and %d cells remaining in the expr matrix.",
                          expr_dim[1], expr_dim[2]))
        
    }
    
    return(infercnv_obj)
        
}


.below_min_mean_expr_cutoff <- function(expr_data, min_mean_expr) {

    average_gene <- rowMeans(expr_data)

    # Find averages above a certain threshold
    indices <- which(average_gene < min_mean_expr)

    return(indices)

}
    



#' @title require_above_min_cells_ref()
#'
#' @description Filters out genes that have fewer than specified number of cells expressing them.
#'
#' @param infercnv_obj infercnv_object
#' 
#' @param min_cells_per_gene int indicating number of cells required per gene for both obs and ref data
#'
#' @return infercnv_obj infercnv_object with corresponding genes removed.
#' 
#' @export
#'

require_above_min_cells_ref <- function(infercnv_obj, min_cells_per_gene) {
    
    genes_passed = which(apply(infercnv_obj@expr.data, 1, function(x) { sum(x>0 & ! is.na(x)) >= min_cells_per_gene}))

    num_genes_total = dim(infercnv_obj@expr.data)[1]
    num_removed = num_genes_total - length(genes_passed)
    if (num_removed > 0) {

        flog.info(sprintf("Removed %d genes having fewer than %d min cells per gene = %g %% genes removed here",
                          num_removed, min_cells_per_gene, num_removed / num_genes_total * 100))
        
        
        if (num_removed == num_genes_total) {

            flog.warn(paste("::All genes removed! Must revisit your data..., cannot continue here."))
            stop(998)
        }
        
        
        infercnv_obj <- remove_genes(infercnv_obj, -1 * genes_passed)
        
                
    }
    else {

        flog.info("no genes removed due to min cells/gene filter")
        
    }
    
    return(infercnv_obj)
    
}


#' @title clear_noise()
#'
#' @description Remove values that are too close to the reference cell expr average and are considered noise.
#' 
#' @param infercnv_obj infercnv_object
#'
#' @param threshold values within reference mean +- threshold are set to zero. 
#'
#' @return infercnv_obj 
#'
#' @export
#'

clear_noise <- function(infercnv_obj, threshold) {
    
    flog.info(paste("********* ::clear_noise:Start. threshold: ", threshold, sep=""))

    if (threshold == 0) {
        return(infercnv_obj); # nothing to do
    }

    ref_idx = get_reference_grouped_cell_indices(infercnv_obj)
    vals = infercnv_obj@expr.data[,ref_idx]
    
    mean_ref_vals = mean(vals)

    infercnv_obj@expr.data <- .clear_noise(infercnv_obj@expr.data, threshold, center_pos=mean_ref_vals)

    return(infercnv_obj)
}


.clear_noise <- function(expr_data, threshold, center_pos=0) {
    
    upper_bound = center_pos + threshold
    lower_bound = center_pos - threshold
    
    expr_data[expr_data > lower_bound & expr_data < upper_bound] = center_pos

    return(expr_data)
}



#' @title clear_noise_via_ref_mean_sd()
#'
#' @description Define noise based on the standard deviation of the reference cell expression data.
#' The range to remove noise would be mean +- sdev * sd_amplifier
#' where sd_amplifier expands the range around the mean to be removed as noise.
#' Data points defined as noise are set to zero.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param sd_amplifier multiplicative factor applied to the standard deviation to alter the noise
#'                     range (default: 1.5)
#' 
#' @export
#'

clear_noise_via_ref_mean_sd <- function(infercnv_obj, sd_amplifier=1.5) {

    ref_idx = get_reference_grouped_cell_indices(infercnv_obj)
    vals = infercnv_obj@expr.data[,ref_idx]
    
    mean_ref_vals = mean(vals)
    
    mean_ref_sd <- mean(apply(vals, 2, function(x) sd(x, na.rm=T))) * sd_amplifier

    upper_bound = mean_ref_vals + mean_ref_sd
    lower_bound = mean_ref_vals - mean_ref_sd
    
    flog.info(paste(":: **** clear_noise_via_ref_quantiles **** : removing noise between bounds: ",
                           lower_bound, "-", upper_bound, sep=" "))


    smooth_matrix <- infercnv_obj@expr.data
        
    smooth_matrix[smooth_matrix > lower_bound & smooth_matrix < upper_bound] = mean_ref_vals
    
    infercnv_obj@expr.data <- smooth_matrix
    
    return(infercnv_obj)
}




# Remove the tails of values of a specific chromosome.
# The smooth_matrix values are expected to be in genomic order.
# If the tail is too large and no contig will be left 1/3 of the
# contig is left.
#
# Args:
# smooth_matrix Smoothed values in genomic order.
#                          Row = Genes, Col = Cells.
# chr Indices of the chr in which the tails are to be removed.
# tail_length Length of the tail to remove on both ends of the
#                        chr indices.
# Returns:
# Indices to remove.
.remove_tails <- function(smooth_matrix, chr, tail_length){

    #flog.info(paste("::remove_tails:Start.", sep=""))
    chr_length <- length(chr)
    if ((tail_length < 3) || (chr_length < 3)){
        return(c())
    }
    if (chr_length < (tail_length * 2)){
         tail_length <- floor(chr_length / 3)
    }
    remove_indices <- chr[1:tail_length]
    remove_indices <- c(remove_indices,
                        chr[ ( (chr_length + 1) - tail_length):
                             chr_length])
    
    return(remove_indices)
}


#' @title smooth_by_chromosome()
#' 
#' @description Smooth expression values for each cell across each chromosome by using a
#' moving average with a window of specified length.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param window_length length of window (number of genes) for the moving average
#'
#' @param smooth_ends perform smoothing at the ends of the chromosomes (default: TRUE)
#'
#' @return infercnv_obj
#' 
#' @export
#'

smooth_by_chromosome <- function(infercnv_obj, window_length, smooth_ends=TRUE) {
    
    gene_chr_listing = infercnv_obj@gene_order[[C_CHR]]
    chrs = unlist(unique(gene_chr_listing))

    for (chr in chrs) {
        chr_genes_indices = which(gene_chr_listing == chr)
        flog.debug(paste0("smooth_by_chromosome: chr: ",chr)) 
        
        input_data = data=infercnv_obj@expr.data[chr_genes_indices, , drop=F]
        flog.debug(paste0("dim subset:", paste(dim(input_data), collapse=",")))
        
        chr_data=infercnv_obj@expr.data[chr_genes_indices, , drop=F]

        if (nrow(chr_data) > 1) {
            smoothed_chr_data = .smooth_window(data=chr_data,
                                               window_length=window_length)
            
            flog.debug(paste0("smoothed data: ", paste(dim(smoothed_chr_data), collapse=",")))
            
            infercnv_obj@expr.data[chr_genes_indices, ] <- smoothed_chr_data
        }
    }
    
    return(infercnv_obj)
}


.smooth_window <- function(data, window_length) {

    flog.info(paste("::smooth_window:Start.", sep=""))
    
    if (window_length < 2){
        flog.warn("window length < 2, returning original unmodified data")
        return(data)
    }

    num_genes <- nrow(data)
    flog.debug(paste("::smooth_window: dim data_sm: ", dim(data_sm), sep=" "))
    
                                        # Fix ends that couldn't be smoothed since not spanned by win/2 at ends.
        # data_sm <- apply(data_sm,
    data_sm <- apply(data,
                     2,
                     .smooth_helper,
                     window_length=window_length)
#                    tail_length=tail_length)
    
    # Set back row and column names
    row.names(data_sm) <- row.names(data)
    colnames(data_sm) <- colnames(data)
    
    return(data_sm)
}


# Helper function for smoothing with a moving average.
#
# Args:
# obs_data: Data to smooth
# window_length:  Length of the window to smooth.
#
# Returns:
# Data smoothed.
##.smooth_helper <- function(obs_data, tail_length) {

.smooth_helper <- function(obs_data, window_length) {
    # strip NAs out and replace after smoothing
    orig_obs_data = obs_data

    nas = is.na(obs_data)

    obs_data = obs_data[!nas]

    obs_length <- length(obs_data)
    end_data <- obs_data

    tail_length = (window_length - 1)/2
    if (obs_length >= window_length) {
        end_data <- .smooth_center_helper(obs_data, window_length)
    }
    
    # end_data will have the end positions replaced with mean values, smoothing just at the ends.

    obs_count <- length(obs_data)

    numerator_counts_vector = c(c(1:tail_length), tail_length + 1, c(tail_length:1))
    
    # defining the iteration range in cases where the window size is larger than the number of genes. In that case we only iterate to the half since the process is applied from both ends.
    iteration_range = ifelse(obs_count > window_length, tail_length, ceiling(obs_count/2))

    for (tail_end in 1:iteration_range) {
        end_tail = obs_count - tail_end + 1

        d_left = tail_end - 1
        d_right = obs_count - tail_end
        d_right = ifelse(d_right > tail_length, tail_length, d_right)

        r_left = tail_length - d_left
        r_right = tail_length - d_right

        denominator = (((window_length - 1)/2)^2 + window_length) - ((r_left * (r_left + 1))/2) - ((r_right * (r_right + 1))/2)

        left_input_vector_chunk = obs_data[1:(tail_end + d_right)]
        right_input_vector_chunk = obs_data[(end_tail - d_right):obs_length]

        numerator_range = numerator_counts_vector[(tail_length + 1 - d_left):(tail_length + 1 + d_right)]

        end_data[tail_end] = sum(left_input_vector_chunk * numerator_range)/denominator
        end_data[end_tail] = sum(right_input_vector_chunk * rev(numerator_range))/denominator
    }

    orig_obs_data[! nas] = end_data  # replace original data with end-smoothed data

    return(orig_obs_data)
}

# Smooth vector of values over the given window length.
#
# Args:
# obs_data Vector of data to smooth with a moving average.
# window_length Length of the window for smoothing.
#        Must be and odd, positive, integer.
#
# Returns:
# Vector of values smoothed with a moving average.
.smooth_center_helper <- function(obs_data, window_length){

    nas = is.na(obs_data)
    vals = obs_data[! nas]

    custom_filter_denominator = ((window_length-1)/2)^2 + window_length
    custom_filter_numerator = c(c(1:((window_length-1)/2)), ((window_length-1)/2)+1, c(((window_length-1)/2):1))

    custom_filter = custom_filter_numerator/rep(custom_filter_denominator, window_length)

#    flog.info(paste("custom filter = ", custom_filter, "\n and window_length =", window_length, "\n and nrow(data) = ", nrow(obs_data), sep=""))

    #smoothed = filter(vals, rep(1 / window_length, window_length), sides=2)
    smoothed = filter(vals, custom_filter, sides=2)

    ind = which(! is.na(smoothed))
    vals[ind] = smoothed[ind]

    obs_data[! nas] = vals

    return(obs_data)
}


#' @title get_average_bounds()
#'
#' @description Computes the mean of the upper and lower bound for the data across all cells.
#'
#' @param infercnv_obj infercnv_object
#'
#' @return (lower_bound, upper_bound)
#' 
#' @export
#'

get_average_bounds <- function (infercnv_obj) {
    
    return(.get_average_bounds(infercnv_obj@expr.data))
    
}

.get_average_bounds <- function(expr_matrix) {
    
    lower_bound <- mean(apply(expr_matrix, 2,
                              function(x) quantile(x, na.rm=TRUE)[[1]]))
    upper_bound <- mean(apply(expr_matrix, 2,
                              function(x) quantile(x, na.rm=TRUE)[[5]]))
    
    return(c(lower_bound, upper_bound))
}

#' @title log2xplus1()
#'
#' @description Computes log(x+1), updates infercnv_obj@expr.data
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#' 
#' @export
#'

log2xplus1 <- function(infercnv_obj) {

    flog.info("transforming log2xplus1()")
    
    infercnv_obj@expr.data <- log2(infercnv_obj@expr.data + 1)

    return(infercnv_obj)
        
}



#' @title invert_log2xplus1()
#'
#' @description Computes 2^x - 1
#' Updates infercnv_obj@expr.data
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#' 
#' @export
#'

invert_log2xplus1 <- function(infercnv_obj) {

    flog.info("inverting log2xplus1()")
    
    infercnv_obj@expr.data <- 2^infercnv_obj@expr.data - 1

    return(infercnv_obj)
}


#' @title invert_log2()
#'
#' @description Computes 2^x
#' Updates infercnv_obj@expr.data
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#' 
#' @export
#'

invert_log2 <- function(infercnv_obj) {

    flog.info("invert_log2(), computing 2^x")
    
    infercnv_obj@expr.data <- 2^infercnv_obj@expr.data

    return(infercnv_obj)
}


#' @title make_zero_NA()
#'
#' @description Converts zero to NA
#' Updates infercnv_obj@expr.data
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#' 
#' @export
#'

make_zero_NA <- function(infercnv_obj) {

    flog.info("make_zero_NA()")
    
    infercnv_obj@expr.data <- infercnv_obj@expr.data[infercnv_obj@expr.data == 0] <- NA
    
    return(infercnv_obj)

}


#' @title transform_to_reference_based_Zscores()
#'
#' @description Computes mean and standard deviation for the reference cells, then uses these values
#' to compute gene Z-scores for both the reference and observation cells.
#'
#' Note, reference cell gene expression values will be centered at zero after this operation.
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#' 
#' @export
#'

transform_to_reference_based_Zscores <- function(infercnv_obj) {

    # center and convert to z-scores
    flog.info(paste("::center_and_Zscore_conversion", sep=""))
    
    # remember, genes are rows, cells are cols
    
    # centering and z-scores based on the reference (normal) cells:
    
    # ref data represent the null distribution
    ref_idx = get_reference_grouped_cell_indices(infercnv_obj)
    
    ref_data = infercnv_obj@expr.data[,ref_idx]
    
    gene_ref_mean = apply(ref_data, 1, function(x) {mean(x, na.rm=T)})
    gene_ref_sd = apply(ref_data, 1, function(x) {sd(x, na.rm=T)})

    # assume at least Poisson level variation
    gene_ref_sd = pmax(gene_ref_sd, gene_ref_mean)
    
    # center all genes at the ref (normal) center:
    infercnv_obj@expr.data = sweep(infercnv_obj@expr.data, 1, gene_ref_mean, FUN="-")
    
    # convert to z-scores based on the ref (normal) distribution
    infercnv_obj@expr.data = sweep(infercnv_obj@expr.data, 1, gene_ref_sd, FUN="/") # make all data z-scores based on the ref data distribution.
    
    
    return(infercnv_obj)
    
}


#' @title mean_center_gene_expr()
#'
#' @description mean-center all gene expression values across all cells
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#' 
#' @export
#'

mean_center_gene_expr <- function(infercnv_obj) {

    flog.info(paste("::centering", sep=""))

    infercnv_obj@expr.data <- sweep(infercnv_obj@expr.data, 1, rowMeans(infercnv_obj@expr.data, na.rm=T), FUN="-")
        
    return(infercnv_obj)
}


#' @title get_reference_grouped_cell_indices()
#'
#' @description Retrieves the matrix indices for the columns correspoinding to the reference cells.
#'
#' @param infercnv_obj infercnv_object
#'
#' @return vector of column indices
#' 
#' @export
#'

get_reference_grouped_cell_indices <- function(infercnv_obj) {

    return( unlist(infercnv_obj@reference_grouped_cell_indices) )

}


#' @title apply_max_threshold_bounds()
#'
#' @description Assumes centered at zero and sets bounds to +- threshold value.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param threshold value to threshold the data
#'
#' @export
#'

apply_max_threshold_bounds <- function(infercnv_obj, threshold) {

    flog.info(paste("::process_data:setting max centered expr, threshold set to: +/-: ", threshold))

    infercnv_obj@expr.data[infercnv_obj@expr.data > threshold] <- threshold
    infercnv_obj@expr.data[infercnv_obj@expr.data < (-1 * threshold)] <- -1 * threshold

    return(infercnv_obj)
}


#' @title remove_genes_at_ends_of_chromosomes()
#'
#' @description Removes genes that are within window_length/2 of the ends of each chromosome.
#' 
#' @param infercnv_obj infercnv_object
#'
#' @param window_length length of the window to use.
#' 
#' @export
#'

remove_genes_at_ends_of_chromosomes <- function(infercnv_obj, window_length) {

    contig_tail = (window_length - 1) / 2
    
    remove_indices <- c()
    gene_chr_listing = infercnv_obj@gene_order[[C_CHR]]
    chrs = unlist(unique(gene_chr_listing))
    for (chr in chrs){
        #flog.info(paste("::process_data:Remove tail contig ",chr, ".", sep=""))
        remove_chr <- .remove_tails(infercnv_obj@expr.data,
                                    which(gene_chr_listing == chr),
                                    contig_tail)
        
        #flog.debug(paste("::process_data:Remove tail - removing indices for chr: ", chr, ", count: ", length(remove_chr), sep=""))

        remove_indices <- c(remove_indices, remove_chr)

    }
    if (length(remove_indices) > 0){

        infercnv_obj = remove_genes(infercnv_obj, remove_indices)

        flog.info(paste("::process_data:Remove genes at chr ends, ",
                        "new dimensions (r,c) = ",
                        paste(dim(infercnv_obj@expr.data), collapse=","),
                        " Total=", sum(infercnv_obj@expr.data, na.rm=TRUE),
                        " Min=", min(infercnv_obj@expr.data, na.rm=TRUE),
                        " Max=", max(infercnv_obj@expr.data, na.rm=TRUE),
                        ".", sep=""))
        
    }
    else {
        flog.error("No genes removed at chr ends.... something wrong here")
        stop(1234)
    }

        
    return(infercnv_obj)

}


#' @title normalize_counts_by_seq_depth()
#'
#' @description Normalizes count data by total sum scaling
#'
#' For single cell data, a typical normalization factor is 1e5, providing counts per 100k total counts.
#' If a normalization factor is not provided, one is estimated based on:
#'     10^round(log10(mean(column_sums))) 
#' 
#' @param infercnv_obj infercnv_object
#'
#' @param normalize_factor  total counts to scale the normalization to (default: NA, computed as described above)
#' 
#' @export
#'

normalize_counts_by_seq_depth <- function(infercnv_obj, normalize_factor=NA) {

    data <- infercnv_obj@expr.data
    
    cs = colSums(data)

    # make fraction of total counts:
    data <- sweep(data, STATS=cs, MARGIN=2, FUN="/")
    
    if (is.na(normalize_factor)) {
        
        normalize_factor = .compute_normalization_factor_from_column_sums(cs)
        
        flog.info(sprintf("Computed total sum normalization factor as: %f", normalize_factor))
                
    } else {
        flog.info(sprintf("Using specified normalization factor: %f", normalize_factor))
    }

    data <- data * normalize_factor

    infercnv_obj@expr.data <- data

    return(infercnv_obj)
        
}

#' @title compute_normalization_factor()
#'
#' @description computes norm factor as:
#'    normalize_factor = 10^round(log10(mean(cs)))
#'
#' @param infercnv_obj infercnv_object
#'
#' @return normalization_factor
#'
#' @export
#'

compute_normalization_factor <- function(infercnv_obj) {

    data <- infercnv_obj@expr.data
    
    cs = colSums(data)
    
    normalize_factor = .compute_normalization_factor_from_column_sums(cs)

    return(normalize_factor)

}
    
#' @keywords internal
#' @noRd
#'
.compute_normalization_factor_from_column_sums <- function(cs) {

    normalize_factor = 10^round(log10(mean(cs)))

    return(normalize_factor)
}



#' @title anscombe_transform()
#'
#' @description Performs Anscombe's transformation:
#'    y = 2 * sqrt(x + 3/8)
#' as per
#' https://en.wikipedia.org/wiki/Anscombe_transform
#'
#' @param infercnv_obj infercnv_object
#'
#' @export
#'

anscombe_transform <- function(infercnv_obj) {

    infercnv_obj@expr.data <- 2 * sqrt(infercnv_obj@expr.data + 3/8)
    
    return(infercnv_obj)
    
}

#' @keywords internal
#' @noRd
#'
add_pseudocount <- function(infercnv_obj, pseudocount) {

    flog.info(sprintf("Adding pseudocount: %g", pseudocount))
        
    infercnv_obj@expr.data = infercnv_obj@expr.data + pseudocount

    return(infercnv_obj)
}
