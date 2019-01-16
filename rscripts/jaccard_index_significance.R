rm(list=ls())
library(jaccard)
source('src/jaccard_index_sig.R')
# higher than expcted (conserved); close to expected; lower than expected (cell type specific) 
# what are the background? All the genes in the chromosome? All the genes assigned to HIMs per chrom?

#df = read.table('data/venn.txt', header=T, row.names = 1, stringsAsFactors = F, sep='\t')
#genes = rownames(df); genes = unlist(sapply(genes, function(z) strsplit(z, ';')))
# step 0. compaute JI indices to select pairs of hims 
# step 1. load the genes, TFs, and chromosome per HIMs
himinfo = extract_gene_tf_chr_him()
# create the background genes; genes in the heterogeneous networks per chromosome 
gene_bg_chr = background_genes()
# create the background TFs which is universal
TFs_bg = unique(unlist(himinfo[['tfs']]))
# load the hims pairs with non-zero JI on genes 
hims_ji = load_him_pairs_JI()
# compute the Jaccard index and pvals only if two overlaps 
hims_ji_2 = compute_ji_pval()
write.table(hims_ji_2, file='inter_results/ji_pval.txt', col.names =T, row.names=F, sep='\t', quote=F)
hims_dyna = cbind(hims_ji, hims_ji_2)
# summrize 
library(ggplot2)
ggplot(hims_dyna, aes(x=log2(jiGene/ expectation_gene), y = -1*log10(pvalue_gene))) + geom_point() + 
    facet_grid(cell1 ~ cell2)
# defnie consered and cell type specific genes 
# conserved
fivenum(hims_dyna$pvalue_gene)
ind1 = hims_dyna$jiGene > 2 * hims_dyna$expectation_gene 
ind2 = hims_dyna$pvalue_gene <=0.001
table(ind1, ind2)
fivenum(hims_dyna[ind1&ind2, 'jiGene'])
# cell type-specific
ind1 = hims_dyna$jiGene < hims_dyna$expectation_gene 
ind2 = hims_dyna$pvalue_gene <=0.05
table(ind1, ind2)
fivenum(hims_dyna[ind1&ind2, 'jiGene'])
# tfs
ind1 = hims_dyna$jiTF > 2* hims_dyna$expectation_tf 
ind2 = hims_dyna$pvalue_tf <=0.05
table(ind1, ind2)
fivenum(hims_dyna[ind1&ind2, 'jiTF'])
ind1 = hims_dyna$jiTF < hims_dyna$expectation_tf 
ind2 = hims_dyna$pvalue_tf <=0.05
table(ind1, ind2)
fivenum(hims_dyna[ind1&ind2, 'jiTF'])
# cell type-specific
# redefine conserved and cell type specific HIMs 
# multiple test adjustment?
# summarize the numbers 
