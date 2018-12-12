#' @include inferCNV.R
NULL



#' @title .mask_DE_genes()
#'
#' @description private function that does the actual masking of expression values in the matrix
#' according to the specified mask value
#'
#' @param infercnv_obj infercnv_object
#'
#' @param all_DE_results  DE results list with structure:
#'                           all_DE_results[[ tumor_type ]] = list(tumor=tumor_type, de_genes=c(geneA, geneB, ...))
#'
#' @param mask_val float value to assign to genes that are in the mask set (complement of the DE gene set)
#'
#' @return infercnv_obj
#' 
#' @keywords internal
#' @noRd
#'

.mask_DE_genes <- function(infercnv_obj, all_DE_results, mask_val) {

    
    all_DE_genes_matrix = matrix(data=FALSE, nrow=nrow(infercnv_obj@expr.data),
                                 ncol=ncol(infercnv_obj@expr.data),
                                 dimnames = list(rownames(infercnv_obj@expr.data),
                                                 colnames(infercnv_obj@expr.data) ) )



    all_DE_genes_matrix[,] = FALSE
    
    ## turn on all normal genes
    all_DE_genes_matrix[, unlist(infercnv_obj@reference_grouped_cell_indices)] = TRUE
    
    for (DE_results in all_DE_results) {
        tumor_type = DE_results$tumor
        genes = DE_results$de_genes
        
    
        all_DE_genes_matrix[rownames(all_DE_genes_matrix) %in% genes,
                            infercnv_obj@observation_grouped_cell_indices[[ tumor_type ]] ] = TRUE
    
    }
        
    infercnv_obj@expr.data[ ! all_DE_genes_matrix ] = mask_val
    
    return(infercnv_obj)
    

}

#' @title mask_non_DE_genes_basic()
#'
#' @description Mask gene expression in infercnv_obj based on expression values that are found
#' statistically significantly different between tumor/normal based on specified
#' statistical test.
#'
#' Uses the expression data stored in infercnv_obj@expr.data
#'
#' @param infercnv_obj infercnv object
#'
#' @param p_val_thresh  p-value threshold for assigning DE statistically significant 
#'
#' @param test.use statistical test to use.  (default: "wilcoxon") alternatives include 'perm' or 't'.
#'
#' @param center_val value to assign to those genes that are not found to be statistically DE.
#'
#' @return infercnv_obj
#' 
#' @export
#'

mask_non_DE_genes_basic <- function(infercnv_obj,
                                    p_val_thresh = 0.05,
                                    test.use="wilcoxon",
                                    center_val=mean(infercnv_obj@expr.data) ) {
    

    all_DE_results = get_DE_genes_basic(infercnv_obj,
                                        p_val_thresh=p_val_thresh,
                                        test.use=test.use)

    #save('all_DE_results', file='all_DE_results.obj')
    
    infercnv_obj <- .mask_DE_genes(infercnv_obj,
                                   all_DE_results,
                                   mask_val=center_val)

    return(infercnv_obj)
}



#' @title get_DE_genes_basic
#'
#' @description Retrieves genes identified as significantly DE based on the corresponding statistical test as applied
#' to the infercnv_obj@expr.data with tumor/normal comparison
#'
#' @param infercnv_obj infercnv object
#'
#' @param p_val_thresh maximum p-value for defining statistical signficance
#'
#' @param test.use  statistical test to use. (default: "t" for t-test)
#' 
#' @return all_DE_results[[condition_pair]] = list(tumor=tumor_type,
#'                                                 normal=normal_type,
#'                                                 pvals=pvals,
#'                                                 de_genes=genes)
#'
#' @export
#'

get_DE_genes_basic <- function(infercnv_obj,
                               p_val_thresh = 0.05,
                               test.use="wilcoxon" # other options:  (wilcoxon, t, perm)
                               ) {
    
    ## Find DE genes by comparing the mutant types to normal types
    normal_types = names(infercnv_obj@reference_grouped_cell_indices)
    tumor_types = names(infercnv_obj@observation_grouped_cell_indices)

    all_DE_results = list()

    statfxns = list()
    statfxns[[ "t" ]] <- function(x, idx1, idx2) {
        vals1 = x[idx1]
        vals2 = x[idx2]

        # useful way of handling tests that may fail:
        # https://stat.ethz.ch/pipermail/r-help/2008-February/154167.html
        res = try(t.test(vals1, vals2), silent=TRUE)
        
        if (is(res, "try-error")) return(NA) else return(res$p.value)
        
    }

    statfxns[[ "perm" ]] <- function(x, idx1, idx2) {

        vals1 = x[idx1]
        vals2 = x[idx2]

        allvals = c(vals1, vals2)
        facts = factor(rep(c("A","B"), c(length(vals1), length(vals2))))

        perm = coin::oneway_test(allvals ~ facts)
        pval = coin::pvalue(perm)
        return(pval)
    }

    statfxns[[ "wilcoxon" ]] <- function(x, idx1, idx2) {

        vals1 = x[idx1]
        vals2 = x[idx2]
                
        w = wilcox.test(vals1, vals2)

        return(w$p.value)
    }
    
    
    statfxn = statfxns[[ test.use ]]
    
    ## turn on only DE genes in tumors
    for (tumor_type in tumor_types) {

        tumor_indices = infercnv_obj@observation_grouped_cell_indices[[ tumor_type ]] 
        
        for (normal_type in normal_types) {
            flog.info(sprintf("Finding DE genes between %s and %s", tumor_type, normal_type))

            normal_indices = infercnv_obj@reference_grouped_cell_indices[[ normal_type ]]


            pvals = apply(infercnv_obj@expr.data, 1, statfxn, idx1=normal_indices, idx2=tumor_indices)
            pvals = unlist(pvals)
            pvals = p.adjust(pvals, method="BH")
            
            names(pvals) = rownames(infercnv_obj@expr.data)

            genes = names(pvals)[pvals<p_val_thresh]
            
            flog.info(sprintf("Found %d genes / %d total as DE", length(genes), length(pvals)))
            
            condition_pair = paste(tumor_type, normal_type, sep=",")
            
            all_DE_results[[condition_pair]] = list(tumor=tumor_type,
                                                    normal=normal_type,
                                                    pvals=pvals,
                                                    de_genes=genes)
            
            
        }
    }
    
    return(all_DE_results)

}

