rm(list=ls())
df = read.table('data/venn.txt', header=T, row.names = 1, stringsAsFactors = F, sep='\t')
cells = colnames(df)
geneinfor = read.table('data/gene_chrom_bin_num_hg19_combined.bed', header=F, sep='\t', row.names = 4, stringsAsFactors = F)
# remove genes on chr 9 
# df rownames has "gene1;gene2" etc. while geneinfor rownames only has "gene".
df.gene = sapply(rownames(df), function(z) gsub(';.*$', '', z))
ind0 = sapply(rownames(df), function(z) grepl(';', z))
ind = geneinfor[df.gene, 1] != 'chr9'
df = df[ind, ]
sum(!ind); # number of genes on chr9
sum(ind0)
sum(ind0 & !ind); # number of genes on chr9 and has ;
## compare two cell types; compute the number of genes that are assigned to HIMs in both cell types 
nogene_in_hims_2cell = function(cell1, cell2, df, type=1){
    ind = df[, cell1] == type & df[, cell2] == type
    return(sum(ind))
}
cellcomb = combn(length(cells), 2)
cellcombname = paste(cells[cellcomb[1, ]], cells[cellcomb[2, ]], sep='_vs_')
n1 = apply(df, 2, function(z) sum(z == 1))

res = matrix(0, nrow=length(cells), ncol=length(cells))
dimnames(res) = list(cells, cells)
for( i in 1:(length(cells) -1 )){
    for(j in (i+1):length(cells)){
        ng = nogene_in_hims_2cell(cell1 = cells[i], cell2=cells[j], df=df, type=1)
        res[i, j] = ng
        res[j, i] = ng
    }
}
res
library(corrplot)
corrplot(res, is.corr = F)
range(res[row(res) != col(res)])

