rm(list=ls())
library(splitstackshape)
library(ggplot2)
source('src/boxdot.R')
df = read.table('data/venn.txt', header=T, row.names = 1, stringsAsFactors = F, sep='\t')
cells = colnames(df)
genes = rownames(df); genes = unlist(sapply(genes, function(z) strsplit(z, ';')))
#### part 1. Essential genes and genes assigned to HIMs across cell types 
####  Housekeeping genes also;
gess = read.table('data/essential_combined.txt', header=F, stringsAsFactors = F)[, 1]
gess = gess[ gess %in% genes]
gene_freq_in_HIM = function(X){
    # be careful with genes in the 10kb bins; gene1;gene2;gene3
    genes = rownames(X)
    freq = rowSums(X==1)
    res = data.frame(genes, nhim=freq, stringsAsFactors = F)
    res = cSplit(res, 'genes', sep=';', direction='long')
    res = as.data.frame(res)
    return(res)
}
himfreq = gene_freq_in_HIM(X=df)
him_freq_ess = himfreq[himfreq$genes %in% gess, ]
# prepare numbers for main text
prop.table(table(him_freq_ess$nhim))
ggplot(data=him_freq_ess, aes(nhim)) + geom_bar(aes(y=(..count..)/sum(..count..)), width=0.2, fill='blue', alpha=0.5) + 
    xlab('# cell types that a essential gene is in a HIM') +  
    theme_classic() + 
    scale_y_continuous(labels=percent) + 
    theme(axis.title.y=element_blank(), axis.title.x=element_text(size=8), axis.text=element_text(size=7)) 
##### part2  cell type spceific genes vs genes assigned to HIMs only in one cell type
### NHEK him 107
genes107 = c('DSC1', 'DSC3', 'DSG1')
df[genes107, ] # DSC1 and DSG1 are hela specifically expressed genes; DSC3 is expressed only at nhek and hela
### Five-way Venn diagram 

#### Validate GO results
## orginally, GO is done on genes uniquely assigned to HIMs in one cell type (cshg)
## now GO analysis is done by controling cell type-specifically expressed genes (cseg)
# how to do it?
#cshg_c_cseg = function(X, cell){
#    # What's the logic???
#    # observations
#        # 1. cshg - cseg cannot reproduce GO terms on cshg
#    # cshg - cseg. Why
#    # definition of cseg
#}
#df.sub = df[rowSums(df != 0) >= 2, ]
#genes = rownames(df.sub)
#genes = unlist(sapply(genes, function(z) strsplit(z, split=';')))
#write.table(genes, file='inter_results/genes_expressed_in_all_5_cell_types.txt', col.names = F, row.names = F, quote=F)
#ind = df.sub[, 'nhek'] == 1 & (rowSums(df.sub == 1) == 1 ) # genes only assigned to HIMs in NHEK 
#genes.nhek = genes[ind]
#genes.nhek = unlist(sapply(genes.nhek, function(z) strsplit(z, split=';')))
#write.table(genes.nhek, file='inter_results/genes_only_assigned_HIMs_in_NHEK.txt', col.names = F, row.names = F, quote=F)
#ind = df.sub[, 'gm12878'] == 1 & (rowSums(df.sub == 1) == 1 ) # genes only assigned to HIMs in NHEK 
#genes.gm = genes[ind]
#genes.gm = unlist(sapply(genes.gm, function(z) strsplit(z, split=';')))
#write.table(genes.gm, file='inter_results/genes_only_assigned_HIMs_in_gm12878.txt', col.names = F, row.names = F, quote=F)
# GO terms in the main text are gone.

#### compute how many genes that are in HIMs in one cell type are acutally only expressed in the cell type
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
    # prepare the results
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
pdf('sup_fig/prop_cseg_in_cshg.pdf', width=2.7, height=2.3)
ggplot(prop, aes(x=cell, y=p)) + geom_bar(stat='identity', width=0.3, fill='blue', alpha=0.3) + 
    xlab('') + ylab('% of cell type-specific expressed genes\nin the genes uniquely assigned to\nHINs in a cell type') + 
    theme_minimal() +
    theme(axis.title.y = element_text(size=8),
          axis.title.x = element_blank(),
          axis.text = element_text(size=7)) 
dev.off()
