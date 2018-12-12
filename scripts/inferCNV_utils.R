
library(tidyverse)
library(futile.logger)

# plot expression density by chromosome for each observation group, reference groups are shown as single 'normal' group.
plot_density_by_chr <- function(infercnv_obj, pdf_filename=NULL, exclude_range=NULL, include_range = NULL, chrs=NULL) {

    ref_group_cell_indices = infercnv:::get_reference_grouped_cell_indices(infercnv_obj)
    
    
    if (is.null(chrs)) {
        chrs = unique(infercnv_obj@gene_order$chr) 
    }
    
    if (! is.null(pdf_filename)) {
        pdf(pdf_filename)
    }


    chr_expr_vals = list()
    
    for (chr in chrs) {

        
        gene_idx = which(infercnv_obj@gene_order$chr == chr)
        
        ref_data_pts = as.numeric(infercnv_obj@expr.data[gene_idx,ref_group_cell_indices])
        
        df = data.frame(class='normal', vals=ref_data_pts)
        
        for (tumor in names(infercnv_obj@observation_grouped_cell_indices) ) {
            
            tumor_cell_idx = infercnv_obj@observation_grouped_cell_indices[[ tumor ]]
            tumor_data_pts = as.numeric(infercnv_obj@expr.data[gene_idx, tumor_cell_idx])
            
            df = rbind(df, data.frame(class=tumor, vals=tumor_data_pts))
        }

        flog.info(sprintf("Plotting data for chr: %s", chr))

        if (! is.null(exclude_range)) {
            excl_range_left = exclude_range[1]
            excl_range_right = exclude_range[2]

            df = df %>% filter(vals < excl_range_left | vals > excl_range_right)
        } else if (! is.null(include_range)) {
            include_range_left = include_range[1]
            include_range_right = include_range[2]

            df = df %>% filter(vals >= include_range_left & vals <= include_range_right)
        }
                
        p = df %>% ggplot(aes(vals, fill=class)) + geom_density(alpha=0.3) + scale_y_continuous(trans='log10', limits=c(1,NA)) + ggtitle(chr)
        plot(p)

        chr_expr_vals[[ chr ]] = df
        
    }

    if (! is.null(pdf_filename)) {
        dev.off()
    }

    return(chr_expr_vals)
    
}



# plot the spike distribution for each specified chromosome in a single density plot
plot_spike_dist <- function(infercnv_obj, chrs) {


    spike_cell_idx = infercnv_obj@observation_grouped_cell_indices[[ 'SPIKE' ]]

    spike_expr = infercnv_obj@expr.data[ , spike_cell_idx ]

    df = data.frame(class='rest', vals=as.numeric(spike_expr[ -1 * which(infercnv_obj@gene_order$chr %in% chrs), ]))

    for (chr in chrs) {

        df = rbind(df, data.frame(class=chr, vals=as.numeric(spike_expr[ which(infercnv_obj@gene_order$chr == chr), ])))
    }

    p = df %>% ggplot(aes(vals, fill=class)) + geom_density(alpha=0.3) + scale_y_continuous(trans='log10', limits=c(1,NA)) + ggtitle('spike')
    plot(p)
    
}

## examine dist of counts of non-zero valued genes per cell per grouping 
plot_dist_counts_expr_genes_by_chr <- function(infercnv_obj, pdf_filename=NULL, chrs=NULL) {

    group_indices = c(infercnv_obj@observation_grouped_cell_indices, infercnv_obj@reference_grouped_cell_indices)

    if (is.null(chrs)) {
        chrs = unique(infercnv_obj@gene_order$chr) 
    }
    
    if (! is.null(pdf_filename)) {
        pdf(pdf_filename)
    }

    gene_counts_dfs = list()
    
    for (chr in chrs) {
        gene_idx = which(infercnv_obj@gene_order$chr == chr)

        df = NULL
        for (group in names(group_indices)) {
            cell_idx = group_indices[[group]]
            expr.data = infercnv_obj@expr.data[gene_idx, cell_idx]
            gene_counts = apply(expr.data, 2, function(x) { sum(x != 0) } )
            if (is.null(df)) {
                df = data.frame(class=group, gene_counts=gene_counts)
            } else {
                df = rbind(df, data.frame(class=group, gene_counts=gene_counts))
            }
        }
        p = df %>% ggplot(aes(gene_counts, fill=class)) + geom_density(alpha=0.3) + ggtitle(chr)
        plot(p)

        gene_counts_dfs[[ chr ]] = df
    }

    if (! is.null(pdf_filename)) {
        dev.off()
    }

    return(gene_counts_dfs)
}



#' takes the mean expr per gene per group
#' returns dataframe with mean_gene_grpA, mean_gene_grpB 
compare_gene_expr_means_by_group_pair <- function(infercnv_obj, groupA, groupB, chr=NULL) {

    group_indices = c(infercnv_obj@observation_grouped_cell_indices, infercnv_obj@reference_grouped_cell_indices)

    group_indices[[ "normal" ]] = infercnv:::get_reference_grouped_cell_indices(infercnv_obj)
    
    expr.data = infercnv_obj@expr.data
    
    if (! is.null(chr)) {
        gene_idx = which(infercnv_obj@gene_order$chr == chr)
        expr.data = expr.data[gene_idx,]
    }
    groupA.expr.data = expr.data[, group_indices[[ groupA ]] ]
    groupB.expr.data = expr.data[, group_indices[[ groupB ]] ]

    groupA.gene_mean = rowMeans(groupA.expr.data)
    groupB.gene_mean = rowMeans(groupB.expr.data)

    #plot(groupA.gene_mean, groupB.gene_mean)
    smoothScatter(groupA.gene_mean, groupB.gene_mean)
    abline(a=0, b=1, col='magenta')
    
    df=data.frame(groupA=groupA.gene_mean, groupB=groupB.gene_mean)

    return(df)
    
}

#' compare spike vs cancer, both to normal

compare_gene_expr_means_spike_vs_cancer_to_normal <- function(infercnv_obj, tumor_type, chr, xlim=NULL, ylim=NULL) {

    df_normal_vs_spike = compare_gene_expr_means_by_group_pair(infercnv_obj, 'normal', 'SPIKE', chr)
    df_tumor_vs_spike = compare_gene_expr_means_by_group_pair(infercnv_obj, 'normal', tumor_type, chr)

    plot(df_tumor_vs_spike[,1], df_tumor_vs_spike[,2], xlab='normal', ylab=tumor_type, xlim=xlim, ylim=ylim)
    points(df_normal_vs_spike[,1], df_normal_vs_spike[,2], col='red')
    abline(a=0,b=1, col='blue')
    
}
                                                              

#' model the mean-to-variance relationship

get_mean_var <- function(infercnv_obj) {

    group_indices = c(infercnv_obj@observation_grouped_cell_indices, infercnv_obj@reference_grouped_cell_indices)

    mean_var_table = NULL
    
    for (group_name in names(group_indices)) {
        flog.info(sprintf("processing group: %s", group_name))
        expr.data = infercnv_obj@expr.data[, group_indices[[ group_name ]] ]
        m = rowMeans(expr.data)
        v = apply(expr.data, 1, var)
        if (is.null(mean_var_table)) {
            mean_var_table = data.frame(g=group, m=m, v=v)
        } else {
            mean_var_table = rbind(mean_var_table, data.frame(g=group, m=m, v=v))
        }
    }

    
    
    return(mean_var_table)
}

plot_mean_var_table <- function(mvtable) {
    s = smooth.spline(log2(mvtable$m+1), log2(mvtable$v+1))
    p = predict(s, log2(mvtable$m+1))
    smoothScatter(log2(mvtable$m+1), log2(mvtable$v+1))
    points(p, col='green', pch='.')
}
