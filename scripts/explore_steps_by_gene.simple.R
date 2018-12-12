#!/usr/bin/env Rscript

options(error = function() traceback(2))

main = function() {

    args<-commandArgs(TRUE)
    
    if (length(args) == 0) {
        stop("\n\n\tUsage:\texplore_steps_by_gene.simple.R gene_name\n\n");
    }
    
    gene_name = args[1]
        
    data_bundles = load_all_data()

    pdf_filename = paste("by_step.", gene_name, ".pdf", sep="")
    pdf(pdf_filename)
    
    make_plots(gene_name, data_bundles)

    quit(save = "no", status = 0, runLast = FALSE)

}


plot_gene_dist = function(gene_name, data_bundle, drop_zeros=F) {

    ref_gene_expr = na.omit(as.numeric(data_bundle$ref_matrix[gene_name,]))
    obs_gene_expr = na.omit(as.numeric(data_bundle$obs_matrix[gene_name,,drop=T]))

    if (drop_zeros) {
        ref_gene_dens = density(ref_gene_expr[ref_gene_expr!=0])
        obs_gene_dens = density(obs_gene_expr[obs_gene_expr!=0])
    }
    else {
        ref_gene_dens = density(ref_gene_expr)
        obs_gene_dens = density(obs_gene_expr)
    }
    xrange = range(ref_gene_dens$x, obs_gene_dens$x)
    yrange = range(ref_gene_dens$y, obs_gene_dens$y)

    plot(ref_gene_dens, xlim=xrange, ylim=yrange, t='l', col='blue',
         main=data_bundle$filename)
    points(obs_gene_dens, t='l', col='green')
}



make_plots = function(gene_name, data_bundles) {

    #dev.new()
    par(mfrow=c(3,3))

    count = 0
    for (data_bundle in data_bundles) {
        count = count + 1
        remove_zeros = F
        if (count < 3) { remove_zeros = T }
        plot_gene_dist(gene_name, data_bundle, remove_zeros)
    }

}

load_all_data = function() {
    
    obj_files = list.files(".", "*.infercnv_obj")
    for (obj_file in obj_files) {
        message(paste("loading file:", obj_file))
        load(obj_file)
    }
    # put variables into the global environment.
    # https://stackoverflow.com/questions/41193543/r-set-all-variables-to-global-environment
    vars <- ls(all = TRUE)
    for (i in 1:length(vars)){
        assign(vars[i], get(vars[i]), envir = .GlobalEnv)
    }
}


plot_ref_obs_comparison = function(data_bundle, gene_name) {
    S0x = data_bundle
    obs = na.omit(as.numeric(S0x$obs_matrix[gene_name,]))
    ref = na.omit(as.numeric(S0x$ref_matrix[gene_name,]))

    den_ref = density(ref)
    den_obs = density(obs)

    xrange = range(den_ref$x, den_obs$x)
    yrange = range(den_ref$y, den_obs$y)
    
    plot(density(ref), t='l', xlim=xrange, ylim=yrange);
    points(density(obs), col='blue', t='l');
    
}


plot_mean_chr_expr_lineplot = function(infercnv_obj,
                                       num_random_cells=0,
                                       ylim=NA, xlim=NA,
                                       plot_separate=F,
                                       sep_obs_types=F,
                                       incl_sd=F,
                                       incl_combined=F) {
    
    
    data_bundle <- make_data_bundle(infercnv_obj)
    
    ref_data= data_bundle$ref_matrix
    obs_data = data_bundle$obs_matrix
    
    if (num_random_cells > 0) {
        ref_rand_cells = sample(x=colnames(ref_data), size=num_random_cells)
        ref_data = ref_data[,ref_rand_cells]

        obs_rand_cells = sample(x=colnames(obs_data), size=num_random_cells)
        obs_data = obs_data[,obs_rand_cells]
    }
    
    mean_expr_ref = rowMeans(ref_data, na.rm=T)
    mean_expr_obs = rowMeans(obs_data, na.rm=T)

    sd_expr_ref = NULL
    sd_expr_obs = NULL
    
    if (incl_sd) {
        sd_expr_ref = apply(ref_data, 1, sd)
        sd_expr_obs = apply(obs_data, 1, sd)

    }

    
    
    num_features = length(mean_expr_ref)
    idx = seq(1, num_features)

    
    if (length(xlim)==1 && is.na(xlim)) {
        xlim = c(1,num_features)
    }
    
    
    if (length(ylim)==1 && is.na(ylim)) {
        ylim = range(mean_expr_ref, mean_expr_obs)

        if (incl_sd) {
            ylim = range(ylim,
                         mean_expr_ref + sd_expr_ref,
                         mean_expr_ref - sd_expr_ref,
                         sd_expr_obs + sd_expr_obs,
                         sd_expr_obs - sd_expr_obs)            
        }

    }
    
    
    plot_chr_expr_line_by_type = function(data) {

        cells = colnames(data)
        types = sapply(cells, function(x) { vals = strsplit(x, "_"); vals[[1]][1]})
        df = data.frame(types=types, name=names(types))
        by_types = split(df, df$types)
        colors = rainbow(length(by_types))
        col_idx = 1
        for (type in names(by_types)) {
            data_for_type = data[,colnames(data) %in% by_types[[type]]$name]
            mean_expr = rowMeans(data_for_type, na.rm=T)
            points(idx, mean_expr, col=colors[col_idx], t='l')
            col_idx = col_idx + 1
        }
    }
    
    
    if (plot_separate) {
        prevpar = par(mfrow=c(2,1))

        
        plot(idx, mean_expr_ref, col='gray', t='l', ylim=ylim, xlim=xlim, main=data_bundle$filename)

        plot(idx, mean_expr_obs, col='salmon', t='l', ylim=ylim, xlim=xlim, main=data_bundle$filename)

        if (incl_sd) {
            points(idx, mean_expr_ref + sd_expr_ref, col='gray', t='l', lty=3)
            points(idx, mean_expr_ref + -1*sd_expr_ref, col='gray', t='l', lty=3)
            
            points(idx, mean_expr_obs + sd_expr_obs, col='salmon', t='l', lty=3)
            points(idx, mean_expr_obs + -1*sd_expr_obs, col='salmon', t='l', lty=3)

        }
        

        if (sep_obs_types) {
            plot_chr_expr_line_by_type(obs_data)
        }

        par(prevpar)

    }
    else {
    
        plot(idx, mean_expr_ref, col='gray', t='l', ylim=ylim, xlim=xlim, main=data_bundle$filename)

        points(idx, mean_expr_obs, col='salmon', t='l')

        if (incl_combined) {
            points(idx, rowMeans(infercnv_obj@expr.data), col='magenta', t='l')
        }
        
        if (incl_sd) {
            points(idx, mean_expr_ref + sd_expr_ref, col='gray', t='l', lty=3)
            points(idx, mean_expr_ref + -1*sd_expr_ref, col='gray', t='l', lty=3)
            
            points(idx, mean_expr_obs + sd_expr_obs, col='salmon', t='l', lty=3)
            points(idx, mean_expr_obs + -1*sd_expr_obs, col='salmon', t='l', lty=3)
            
        }
        


        if (sep_obs_types) {
            plot_chr_expr_line_by_type(obs_data)
        }
        
    }
}




# plot number of gene-expressing cells by chr
plot_chr_num_cells_lineplot = function(infercnv_obj) {

    data_bundle <- make_data_bundle(infercnv_obj)
    
    ref_num_cells_expr = apply(data_bundle$ref_matrix, 1, function(x) { x[is.na(x)] = 0; sum(x!=0)} )
    obs_num_cells_expr = apply(data_bundle$obs_matrix, 1, function(x) { x[is.na(x)] = 0; sum(x!=0)} )

    idx = seq(1, length(ref_num_cells_expr))
    
    ylim = range(ref_num_cells_expr, obs_num_cells_expr, na.rm=T)

    oldpar = par(mfrow=c(2,1))
    plot(idx, ref_num_cells_expr, col='gray', t='l', ylim=ylim)
    plot(idx, obs_num_cells_expr, col='blue', t='l', ylim=ylim)
    
    par(oldpar)

}

make_data_bundle <- function(infercnv_obj) {

    ref_indices <- unlist(infercnv_obj@reference_grouped_cell_indices)
    obs_indices <- unlist(infercnv_obj@observation_grouped_cell_indices)

    data_bundle <- list()
    data_bundle$ref_matrix <- infercnv_obj@expr.data[,ref_indices]
    data_bundle$obs_matrix <- infercnv_obj@expr.data[,obs_indices]

    return(data_bundle)
    
}


boxplot_mean_expr_distr <- function(infercnv_obj, by="gene" ) { # or cell
    db = make_data_bundle(infercnv_obj)

    if (by == "gene") {
        boxplot(rowMeans(db$ref_matrix), rowMeans(db$obs_matrix), outline=F, names=c('ref', 'obs'))
    }
    else {
        # by 'cell'
        boxplot(colMeans(db$ref_matrix), colMeans(db$obs_matrix), outline=F, names=c('ref', 'obs'))
    }
}


#if (!interactive()) {
#    main()
#}

