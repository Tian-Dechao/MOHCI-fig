rm(list=ls())
library(splitstackshape)
library(ggplot2)
source('src/boxdot.R')
source('src/cell_type_specific_genes.R')
## extract cell type specific genes 
ifile='data/expression/PH_gene_expression.txt'
cells = c('gm12878', 'hela', 'huvec', 'k562', 'nhek')
CELLs = c('GM12878', 'HeLa', 'HUVEC', 'K562', 'NHEK'); names(CELLs) = cells
gcs = load_rna_cs(file=ifile, cs=cells)
df = read.table('data/venn.txt', header=T, row.names = 1, stringsAsFactors = F, sep='\t')
genes = rownames(df); genes = unlist(sapply(genes, function(z) strsplit(z, ';')))
#### part 1. Essential genes and genes assigned to HIMs across cell types 
####  Housekeeping genes from Human housekeeping genes, revisited; 
gess = read.table('data/essential_combined.txt', header=F, stringsAsFactors = F)[, 1]
gess = gess[ gess %in% genes]
gess_k = read.table('data/essential_genes_k562.txt', header=F, stringsAsFactors = F)[, 1]
gess_k = gess_k[ gess_k %in% genes]
ghk = read.table('data/HK_genes.txt', header=F, stringsAsFactors = F)[, 1]
ghk = ghk[ghk %in% genes]
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
him_freq_ess_k = himfreq[himfreq$genes %in% gess_k, ]
him_freq_hk = himfreq[himfreq$genes %in% ghk, ]
# prepare numbers for main text
round(prop.table(table(him_freq_ess$nhim)) * 100, 2)
round(prop.table(table(him_freq_hk$nhim)) * 100, 2)
#30.45 + 36.77 
#31.84 +  36.22 
# data frame for essential genes and housekeeping genes
res = rbind(data.frame(him_freq_ess, type='Essential genes'), data.frame(him_freq_hk, type='HK genes'))
# data frame for cell type specific genes 
res_cs = c()
for(cell in cells){
    tmp = himfreq[himfreq$genes %in% gcs[[cell]], ]
    tmp = data.frame(tmp, type=cell)
    res_cs = rbind(res_cs, tmp)
}
compute_prop_per_category = function(df){
    cats = unique(df$type)
    res = c()
    for(c in cats){
        tmp = as.data.frame(prop.table(table(df[df$type == c, 'nhim'])))
        tmp = data.frame(type=c, tmp)
        res = rbind(res, tmp)
    }
    colnames(res) = c('type', 'nhim', 'prop')
    res$nhim = factor(res$nhim, levels=5:0)
    return(res)
}
res_cs2 = compute_prop_per_category(res_cs)
res2 = compute_prop_per_category(res)
res3 = rbind(res2, res_cs2)
pdf('main_fig/assignment_hims_gene_types.pdf', height=3, width=3)
ggplot(data=res3, aes(x=type, y=prop, fill=nhim)) + geom_bar(stat='identity') + 
    theme_classic() + 
    theme(legend.text = element_text(size=7), legend.title = element_text(size=8)) + 
    ylab('Proportion of genes') + 
    xlab('Cell type specific genes') +
    scale_fill_discrete(name='# cell\ntypes') + 
    scale_x_discrete(labels=CELLs) + 
    scale_y_continuous(labels=percent) + 
    theme(axis.title=element_text(size=8), axis.text=element_text(size=7), axis.text.x= element_text(angle = 45, hjust=1)) 
dev.off()
#ggplot(data=him_freq_ess_k, aes(nhim)) + geom_bar(aes(y=(..count..)/sum(..count..)), width=0.2, fill='blue', alpha=0.5) + 
pdf('main_fig/assiment_hims_ess_hk_genes.pdf', height=3, width=3)
ggplot(data=res, aes(nhim, fill=type)) + geom_bar(position='dodge', aes(y=..prop..)) + 
    xlab('# cell types that a gene is in a HIM') +  
    scale_fill_manual(values=c('#d95f02', '#7570b3')) + 
    theme_classic() + 
    theme(legend.position = c(0.25, 0.85), legend.text = element_text(size=7), legend.title = element_text(size=8)) + 
    scale_y_continuous(labels=percent) + 
    theme(axis.title.y=element_blank(), axis.title.x=element_text(size=8), axis.text=element_text(size=7)) 
dev.off()
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
