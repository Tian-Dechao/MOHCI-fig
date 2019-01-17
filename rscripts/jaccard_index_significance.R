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
himids = himinfo[['himids']]; himids2 = gsub('_.*', '', himids); table(himids2)
# create the background genes; genes in the heterogeneous networks per chromosome 
gene_bg_chr = background_genes()
# create the background TFs which is universal
TFs_bg = unique(unlist(himinfo[['tfs']]))
# load the hims pairs with non-zero JI on genes 
hims_ji = load_him_pairs_JI()
# compute the Jaccard index and pvals only if two overlaps 
#hims_ji_2 = compute_ji_pval()
#write.table(hims_ji_2, file='inter_results/ji_pval_fix_x.txt', col.names =T, row.names=F, sep='\t', quote=F)
#hims_ji_2 = read.table('inter_results/ji_pval_fix_x.txt', header=T, sep='\t')
#hims_dyna = cbind(hims_ji, hims_ji_2)
# hypergeometric test 
hims_hyper = compute_ji_pval_hypergeo()
hims_hyper_adj = apply(hims_hyper, 2, function(z) p.adjust(z, method='bonferroni'))
hims_dyna = cbind(hims_ji, hims_hyper_adj)
head(hims_dyna)
