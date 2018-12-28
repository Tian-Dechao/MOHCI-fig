rm(list=ls())
source('src/boxdot.R')
df = read.table('data/venn.txt', header=T, row.names = 1, stringsAsFactors = F, sep='\t')
cells = colnames(df)
# pick up genes that are expressed in at least n cell types
df.sub = df[rowSums(df != 0) >= 2, ]
genes = rownames(df.sub)
genes = unlist(sapply(genes, function(z) strsplit(z, split=';')))
write.table(genes, file='inter_results/genes_expressed_in_all_5_cell_types.txt', col.names = F, row.names = F, quote=F)
ind = df.sub[, 'nhek'] == 1 & (rowSums(df.sub == 1) == 1 ) # genes only assigned to HIMs in NHEK 
genes.nhek = genes[ind]
genes.nhek = unlist(sapply(genes.nhek, function(z) strsplit(z, split=';')))
write.table(genes.nhek, file='inter_results/genes_only_assigned_HIMs_in_NHEK.txt', col.names = F, row.names = F, quote=F)
ind = df.sub[, 'gm12878'] == 1 & (rowSums(df.sub == 1) == 1 ) # genes only assigned to HIMs in NHEK 
genes.gm = genes[ind]
genes.gm = unlist(sapply(genes.gm, function(z) strsplit(z, split=';')))
write.table(genes.gm, file='inter_results/genes_only_assigned_HIMs_in_gm12878.txt', col.names = F, row.names = F, quote=F)
# GO terms in the main text are gone.

# compute how many genes that care in HIMs in one cell type are acutally only expressed in the cell type
cseg2cshg = function(X, cell){
    # NO. of genes that are only in the GRN in the cell type
    ind = (rowSums(X!=0) == 1) & (X[, cell] != 0)
    n1 = sum(ind)
    # NO. of genes that are only assigned to HIMs in the cell type
    ind2 = (rowSums(X==1) == 1) & (X[, cell] == 1)
    n2 = sum(ind2)
    # overlap between the two types of genes
    ind3 = ind & ind2
    n3 = sum(ind3)
    
    res = c(n1, n2, n3)
    names(res) = c('cseg', 'cshg', 'csehg')
    return(res)
}
res = sapply(cells, function(z) cseg2cshg(X=df, cell=z))
# proportion of cell type specifically expressed genes in the genes uniquely assigned to HIMs in the cell type
p = res[3, ] / res[2, ] * 100
p = round(p, 2)
p; range(p)
#Note that less than half (27.12\%-48.36\%) of the genes uniquely assigned to HIMs  are cell type-specific expressed genes.
library(ggplot2)
prop = data.frame(cell=rename_features_v2(names(p), nl=nl), prop=p, stringsAsFactors = F)
pdf('sup_fig/prop_cseg_in_cshg.pdf', width=3, height=3)
ggplot(prop, aes(x=cell, y=p)) + geom_bar(stat='identity', width=0.3, fill='blue', alpha=0.3) + 
    xlab('') + ylab('% of cell type-specific expressed genes\nin the genes uniquely assigned to HINs in a cell type') + 
    theme_minimal() +
    theme(axis.title.y = element_text(size=8),
          axis.title.x = element_blank(),
          axis.text = element_text(size=7)) 
dev.off()
