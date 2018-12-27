rm(list=ls())
df = read.table('data/venn.txt', header=T, row.names = 1, stringsAsFactors = F, sep='\t')
cells = colnames(df)
# pick up genes that are expressed in at least n cell types
df.sub = df[rowSums(df != 0) >= 2, ]
genes = rownames(df.sub)
genes = unlist(sapply(genes, function(z) strsplit(z, split=';')))
write.table(genes, file='data/genes_expressed_in_all_5_cell_types.txt', col.names = F, row.names = F, quote=F)
ind = df.sub[, 'nhek'] == 1 & (rowSums(df.sub == 1) == 1 ) # genes only assigned to HIMs in NHEK 
genes.nhek = genes[ind]
genes.nhek = unlist(sapply(genes.nhek, function(z) strsplit(z, split=';')))
write.table(genes.nhek, file='data/genes_only_assigned_HIMs_in_NHEK.txt', col.names = F, row.names = F, quote=F)
ind = df.sub[, 'gm12878'] == 1 & (rowSums(df.sub == 1) == 1 ) # genes only assigned to HIMs in NHEK 
genes.gm = genes[ind]
genes.gm = unlist(sapply(genes.gm, function(z) strsplit(z, split=';')))
write.table(genes.gm, file='data/genes_only_assigned_HIMs_in_gm12878.txt', col.names = F, row.names = F, quote=F)
# GO terms in the main text are gone.

# compute how many genes that care in HIMs in one cell type are acutally only expressed in the cell type
ind1 = df[, 'gm12878'] == 1 & (rowSums(df == 1) == 1)
ind2 = ind1 & ( rowSums(df) == 1)
table(ind1)
table(ind2)
table(ind1, ind2)
