#!/usr/bin/env Rscript

library("infercnv")

# create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix="../oligodendroglioma_expression_downsampled.counts.matrix",
                                    annotations_file="../oligodendroglioma_annotations_downsampled.txt",
                                    delim="\t",
                                    gene_order_file="../gencode_downsampled.txt",
                                    ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)"))

out_dir="output_dir.no_spike"
# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir, 
                             cluster_by_groups=T, 
                             plot_steps=F,
                             include.spike=F  # used for final scaling to fit range (0,2) centered at 1.
                             )

