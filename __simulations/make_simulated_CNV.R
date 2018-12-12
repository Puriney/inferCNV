#!/usr/bin/env Rscript


filename="data/normal.counts.matrix"

data = read.table(filename, header=T, row.names=1)


ordered_genes = read.table("../example/full_dataset/gencode_v19_gene_pos.txt")
ordered_gene_names = ordered_genes[,1]
data = data[rownames(data) %in% ordered_gene_names,]
ordered_gene_names = ordered_gene_names[ordered_gene_names %in% rownames(data)]
ordering = match(ordered_gene_names, rownames(data))
data = data[ordering,]



mutate_gene_expr = function(data, genes_to_mutate, multiplier, colname_prefix) {

    data_m = data
    data_m[genes_to_mutate,] =  data_m[genes_to_mutate,] * multiplier
    colnames(data_m) = paste(colname_prefix, colnames(data_m), sep="")
    
    return(data_m)

}

sample_annots = data.frame(patient=colnames(data), type="normal")

# lower ranges

genes_to_mutate = 1000:3000
data_pt25 = mutate_gene_expr(data, genes_to_mutate, 0.25, "A_")
sample_annots = rbind(sample_annots, data.frame(patient=colnames(data_pt25), type="Apt25"))


genes_to_mutate = 5000:7000
data_pt5 = mutate_gene_expr(data, genes_to_mutate, 0.5, "B_")
sample_annots = rbind(sample_annots, data.frame(patient=colnames(data_pt5), type="Bpt5"))


# upper ranges
genes_to_mutate = 9000:11000
data_1pt5 = mutate_gene_expr(data, genes_to_mutate, 1.5, "C_")
sample_annots = rbind(sample_annots, data.frame(patient=colnames(data_1pt5), type="C1pt5"))

genes_to_mutate = 13000:15000
data_1pt75 = mutate_gene_expr(data, genes_to_mutate, 1.75, "D_")
sample_annots = rbind(sample_annots, data.frame(patient=colnames(data_1pt75), type="D1pt75"))



## make new data
new_data = cbind(data, data_pt25, data_pt5, data_1pt5, data_1pt75)

sample_annots = cbind(sample_annots, t2=sample_annots$type)

# write outputs

write.table(new_data, "sim.data", quote=F, sep="\t")

write.table(sample_annots, "sim.sample.annots.txt", quote=F, sep="\t", row.names=FALSE, col.names=FALSE)

