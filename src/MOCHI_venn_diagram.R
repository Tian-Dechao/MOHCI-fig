rm(list=ls())
library(venn)
comp_cell_venn <- function(gene_table, pdf_file) {
#gene_table <- read.table(file=filename, header = TRUE, sep = '\t')
    cat(colnames(gene_table))
    colnames(gene_table) <- c("gene", "GM12878", "HeLa", "HUVEC", "NHEK", "K562")
    gene_table = gene_table[order(gene_table$gene),]
    gene_table = gene_table[!duplicated(gene_table$gene),]
    venn_table <- gene_table[, 2:6]
    rownames(venn_table) <- gene_table[, 1]
    venn_matrix <- as.matrix(venn_table)
    venn_matrix[venn_matrix>1] <- 0
    venn_table <- as.data.frame(venn_matrix)
    pdf(file=pdf_file, width=5, height=5)
    venn(venn_table, ilabel = TRUE, zcolor = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"))
    dev.off()
    gene_table$label <- paste(venn_table$GM12878, venn_table$HeLa, venn_table$HUVEC, venn_table$NHEK, venn_table$K562, sep='')
    gene_table$label_ori <- paste(gene_table$GM12878, gene_table$HeLa, gene_table$HUVEC, gene_table$NHEK, gene_table$K562, sep='')
    return(gene_table)
}

# compare cell types
data = read.table('data/venn.txt', header=T, sep='\t')
data2 = data; data2[data2 == 2]  = 1
venn_table <- comp_cell_venn(data2, "sup_fig/5_cell_venn_grn.pdf")
venn_table <- comp_cell_venn(data, "sup_fig/5_cell_venn_him.pdf")
#
venn_table <- comp_cell_venn(paste(anno_result, "venn.txt", sep='/'), paste(prefix, "5_cell_venn.pdf", sep='/'))
venn_table <- comp_cell_venn(paste(anno_result, "venn.txt", sep='/'), paste(prefix, "5_cell_venn.pdf", sep='/'))
go_select_gene(venn_table, c("11111", "11110", "11101", "11011", "10111", "01111"), paste(prefix, "5_way_shared.txt", sep='/'))
go_select_gene(venn_table, c("11111"), paste(prefix, "5_way_shared_strict.txt", sep='/'))
go_select_gene(venn_table, c("10000"), paste(prefix, "5_way_GM12878_uniq.txt", sep='/'))
go_select_gene(venn_table, c("01000"), paste(prefix, "5_way_HeLa_uniq.txt", sep='/'))
go_select_gene(venn_table, c("00100"), paste(prefix, "5_way_HUVEC_uniq.txt", sep='/'))
go_select_gene(venn_table, c("00010"), paste(prefix, "5_way_NHEK_uniq.txt", sep='/'))
go_select_gene(venn_table, c("00001"), paste(prefix, "5_way_K562_uniq.txt", sep='/'))
write.table(unlist(lapply(as.vector(venn_table[, 'gene']), function(x) unlist(strsplit(x, ';')))), file=paste(prefix, "5_way_background.txt", sep='/'), col.names = FALSE, row.names = FALSE, quote = FALSE)
