##### significance of overlapped genes between genes that are assigned to HIMs in two cell types. Hypergeometric test is used.
rm(list=ls())
source('src/boxdot.R') # rename cell types
df = read.table('data/venn.txt', header=T, row.names = 1, stringsAsFactors = F, sep='\t')
cells = colnames(df)
geneinfor = read.table('data/gene_chrom_bin_num_hg19_combined.bed', header=F, sep='\t', row.names = 4, stringsAsFactors = F)
# remove genes on chr 9 
# df rownames has "gene1;gene2" etc. while geneinfor rownames only has "gene".
df.gene = sapply(rownames(df), function(z) gsub(';.*$', '', z))
ind0 = sapply(rownames(df), function(z) grepl(';', z))
ind = geneinfor[df.gene, 1] != 'chr9'
df = df[ind, ]
## compare two cell types; compute the number of genes that are assigned to HIMs in both cell types 
nogene_in_hims_2cell = function(cell1, cell2, df, type=1){
    # hypergeometric test which ajust the differnece in the number of genes among the different cell types
    ind = df[, cell1] != 0 & df[, cell2] != 0
    ind1 = df[, cell1] == type
    ind2 = df[, cell2] == type
    N = sum(ind) #  # white + black balls in the urn 
    N1 = sum(ind1 & ind) # white balls; m
    N2 = sum(ind2 & ind) # balls draw; k
    N3 = sum(ind1 & ind2 & ind) # whilte balls draw; x
    pval = phyper(q=N3-1, m=N1, n=N-N1, k=N2, lower.tail = F)
    return(pval)
}

#cellcomb = combn(length(cells), 2)
#cellcombname = paste(cells[cellcomb[1, ]], cells[cellcomb[2, ]], sep='_vs_')
#n1 = apply(df, 2, function(z) sum(z == 1))

res = matrix(0, nrow=length(cells), ncol=length(cells))
for( i in 1:(length(cells) -1 )){
    for(j in (i+1):length(cells)){
        ng = nogene_in_hims_2cell(cell1 = cells[i], cell2=cells[j], df=df, type=1)
        res[i, j] = ng
        res[j, i] = ng
    }
}
res[row(res) == col(res)] = min(res[row(res) != col(res)]) # random number
# transfer res
res = -1 * log10(res)
library(corrplot)
cells2 = rename_features_v2(x=cells, nl=nl)
dimnames(res) = list(cells2, cells2)
pdf('sup_fig/gene_assigment_2_hims_agreements_2cells.pdf', width=3, height=3)
corrplot(res, is.corr = F, order='hclust', diag=F, cl.lim=c(0, max(res)), type='upper', mar=rep(0, 4), tl.cex=0.8,
         cl.length=5, cl.ratio=0.4, cl.align.text = 'c', cl.cex=0.6, title=NULL)
dev.off()
# Pairwise comparison shows that GM12878 and K562 (both are blood cell lines) 
# share the highest number of genes that are assigned to HIMs, which reconfirms biological relavance of HIMs.
