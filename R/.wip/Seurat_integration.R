#' @include inferCNV.R
NULL

#' @title mask_non_DE_genes_via_Seurat
#'
#' @description Masks expression data values for those genes defined as *not* differentially expressed.
#' Seurat FindMarkers() is used to define differentially expressed genes.
#' All tumor types are each compared pairwise to the normal cell type groups.
#' Any genes not found differentially expressed according to the defined parameters in any
#' of the tumor/normal comparisons are set to the center_val, effectively masking them.
#'
#' @param infercnv_obj infercnv object
#'
#' @param thresh.use  from Seurat::FindMarkers : "Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25 Increasing thresh.use speeds up the function, but can miss weaker signals."
#'
#' @param min.pct   from Seurat::FindMarkers : "only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed."
#'
#' @param p_val_thresh maximum allowed p-value
#'
#' @param test.use from Seurat::FindMarkers : "Denotes which test to use. Seurat currently implements "bimod" (likelihood-ratio test for single cell gene expression, McDavid et al., Bioinformatics, 2013, default), "roc" (standard AUC classifier), "t" (Students t-test), and "tobit" (Tobit-test for differential gene expression, as in Trapnell et al., Nature Biotech, 2014), 'poisson', and 'negbinom'. The latter two options should only be used on UMI datasets, and assume an underlying poisson or negative-binomial distribution"
#'
#' @param center_val value to set the non-DE gene expression values to in the matrix.  By default, sets to the mean of all input expression values.
#'
#' @param normalize_factor normalization factor used by Seurat in total sum count scaling. Typically set to 1e4 or 1e5. Default: estimates from the input data.
#' 
#' @export
#'

mask_non_DE_genes_via_Seurat <- function(infercnv_obj,
                                         thresh.use=0.25,
                                         min.pct=0.50,
                                         p_val_thresh = 0.05,
                                         test.use='bimod',
                                         center_val = mean(infercnv_obj@expr.data),
                                         normalize_factor=compute_normalization_factor(infercnv_obj) ) {
    
    
  

    all_DE_results = get_DE_genes_via_Seurat(infercnv_obj,
                                             thresh.use=thresh.use,
                                             min.pct=min.pct,
                                             p_val_thresh=p_val_thresh,
                                             test.use=test.use,
                                             normalize_factor=normalize_factor)


    infercnv_obj <- .mask_DE_genes(infercnv_obj, all_DE_results, mask_val=center_val)

    return(infercnv_obj)

}



#' @title seurat_from_infercnv_obj()
#'
#' @description Generate a Seurat object from an infercnv_obj
#'
#' @param infercnv_obj infercnv_object
#'
#' @param normalize_factor scale factor for Seurat
#'
#' @export
#'


seurat_from_infercnv_obj <- function(infercnv_obj, normalize_factor) {

    seurat_obj = CreateSeuratObject(raw.data = infercnv_obj@count.data, project="infcnv")
    
    seurat_obj = NormalizeData(object=seurat_obj, normalization.method = "LogNormalize",
                               scale.factor = normalize_factor)

    ## store cell type info as metadata

    all_cell_classes = c(infercnv_obj@reference_grouped_cell_indices, infercnv_obj@observation_grouped_cell_indices)

    cell_classes = names(all_cell_classes)

    cell_column_annotations = rep(NA, ncol(infercnv_obj@expr.data))

    for (cell_class in cell_classes) {
        cell_column_annotations[ all_cell_classes[[cell_class]] ] <- cell_class
    }
    cell_column_annotations = as.factor(cell_column_annotations)
    names(cell_column_annotations) = colnames(seurat_obj@data)


    patient.ident = seurat_obj@ident

    seurat_obj = StashIdent(seurat_obj, save.name='patient.ident')

    seurat_obj@ident <- as.factor(cell_column_annotations)

    return(seurat_obj)


}



#' @title get_DE_genes_via_Seurat()
#'
#' 
#' @description Seurat FindMarkers() is used to define differentially expressed genes.
#' All tumor types are each compared pairwise to the normal cell type groups.
#' Any genes not found differentially expressed according to the defined parameters in any
#' of the tumor/normal comparisons are set to the center_val, effectively masking them.
#'
#' @param infercnv_obj infercnv object
#'
#' @param thresh.use  from Seurat::FindMarkers : "Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25 Increasing thresh.use speeds up the function, but can miss weaker signals."
#'
#' @param min.pct   from Seurat::FindMarkers : "only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed."
#'
#' @param p_val_thresh maximum allowed p-value
#'
#' @param test.use from Seurat::FindMarkers : "Denotes which test to use. Seurat currently implements "bimod" (likelihood-ratio test for single cell gene expression, McDavid et al., Bioinformatics, 2013, default), "roc" (standard AUC classifier), "t" (Students t-test), and "tobit" (Tobit-test for differential gene expression, as in Trapnell et al., Nature Biotech, 2014), 'poisson', and 'negbinom'. The latter two options should only be used on UMI datasets, and assume an underlying poisson or negative-binomial distribution"
#'
#' @param normalize_factor normalization factor used by Seurat in total sum count scaling. Typically set to 1e4 or 1e5. Default: estimates from the input data.
#'
#' @return DE_results_list, structured like so:
#'
#'             all_DE_results[[condition_pair]] = list(tumor=tumor_type,
#'                                                     normal=normal_type,
#'                                                     markers=markers,
#'                                                     de_genes=genes)
#'
#' @export
#'

get_DE_genes_via_Seurat <- function(infercnv_obj,
                                    thresh.use=0.25,
                                    min.pct=0.50,
                                    p_val_thresh = 0.05,
                                    test.use='bimod',
                                    normalize_factor=1e5 ) {
    

    seurat_obj = seurat_from_infercnv_obj(infercnv_obj, normalize_factor)
    
    ## Find DE genes by comparing the mutant types to normal types
    normal_types = names(infercnv_obj@reference_grouped_cell_indices)
    tumor_types = names(infercnv_obj@observation_grouped_cell_indices)

    all_DE_results = list()
    
    ## turn on only DE genes in tumors
    for (tumor_type in tumor_types) {
        for (normal_type in normal_types) {
            flog.info(sprintf("Finding DE genes between %s and %s", tumor_type, normal_type))
            markers = FindMarkers(seurat_obj,
                                  ident.1=tumor_type,
                                  ident.2=normal_type,
                                  thresh.use=thresh.use,
                                  min.pct=min.pct,
                                  only.pos=FALSE,
                                  test.use=test.use)
            

            
            genes = rownames(markers[markers$p_val<0.05,])
            flog.info(sprintf("Found %d genes as DE", length(genes)))

            condition_pair = paste(tumor_type, normal_type, sep=",")
            
            all_DE_results[[condition_pair]] = list(tumor=tumor_type,
                                                    normal=normal_type,
                                                    markers=markers,
                                                    de_genes=genes)
            
            
        }
    }
    
    return(all_DE_results)

}


