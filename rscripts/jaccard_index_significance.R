rm(list=ls())
library(jaccard)
source('src/jaccard_index_sig.R')
#df = read.table('data/venn.txt', header=T, row.names = 1, stringsAsFactors = F, sep='\t')
#genes = rownames(df); genes = unlist(sapply(genes, function(z) strsplit(z, ';')))
# step 0. compaute JI indices to select pairs of hims 
# step 1. load the genes, TFs, and chromosome per HIMs
himinfo = extract_gene_tf_chr_him()
himids = extraact_himid_per_cell()
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
# extract conserved HIMs per cell line 
stable_hims = hims_stable()
# create a new variable to denote stable or cs hims
dyna_cat = ifelse(himids$himid %in% stable_hims, 'stable', 'cs')
himids  = data.frame(himids, conserve_staVScs=dyna_cat, stringsAsFactors = F)
res = aggregate(conserve_staVScs~cell, himids, FUN=function(z) round(prop.table(table(z)) * 100, 2))
res
range(res$conserve_staVScs[, 2])
range(res$conserve_staVScs[, 1])
subset(himids, himid=='k562_712')
write.table(himids, 'inter_results/him_staVScs_ji_pval.txt', col.names = T, row.names = F, sep='\t', quote=F)