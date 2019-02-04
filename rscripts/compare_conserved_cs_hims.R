rm(list=ls())
#library(jaccard)
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
# compute the expected jiTF and jiGene between two HIMs
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

### draw supplementary figures
source('src/boxdot.R')
source('src/load_him_summary_allinone.R')
# replace the conserve_staVScs with the new one
rownames(himids) = himids$himid
X$conserve_staVScs = himids[rownames(X), 'conserve_staVScs']
feat = c('gene_expression_mean', 'csgp')
nrow=length(feat); ncol=5
pdf(file='sup_fig/sup_stable_vs_cs_functional.pdf', width=1.2 * ncol, height=1.8 * nrow)
fi = 'conserve_staVScs'
Xf = X[, c(feat, 'cell', fi)]
print(boxdot2(df=Xf, fs='cell', fn=fi, g1='stable', g2='cs', ylabsize=1.1, pvalalter = F, psize=1))
dev.off()

# Apercent for constitutive and cell type-specific HIMs
cutoffs = c(0, 50, 80,   99.99999999999,101)
labels = c('[0, 50)', '[50, 80)', '[80, 100)', '100' )
Acon = extract_freq(df=subset(X, conserve_staVScs == 'stable'), cutoffs=cutoffs, labels=labels)
Acon[, 'conserve_staVScs'] = 'stable'
Acs = extract_freq(df=subset(X, conserve_staVScs == 'cs'), cutoffs=cutoffs, labels=labels)
Acs[, 'conserve_staVScs'] = 'cs'
Adynamic = rbind(Acon, Acs)
Adynamic$conserve_staVScs = factor(Adynamic$conserve_staVScs, levels=c('stable', 'cs'), ordered = T)
# report some numbers
Acon[Acon$category == '100', 'proportion'] - Acs[Acs$category == '100', 'proportion']
# compute pvalue
cells = unique(Adynamic$cell)
pval = c()
for( i in 1:length(cells)){
    tab = cbind(Acon[Acon$cell== cells[i], 'proportion'] * sum(X$conserve_staVScs == 'stable' & X$cell== cells[i], na.rm=T), 
                Acs[Acs$cell == cells[i], 'proportion'] * sum(X$conserve_staVScs == 'cs' & X$cell== cells[i], na.rm=T))
    tab = tab / 100
    pval[i] = chisq.test(tab)$p.value
}
pval = format(pval, scientific=T, digits=3); pval = paste('P=', pval, sep='')
AdynamicPval = data.frame(cell=cells, pval=pval, x=2.5, y=60)
col_list = c('stable' = "#2b83ba", 'cs' = "#d7191c")
pdf('sup_fig/Apercent_stable_cs.pdf', width=5.85, height=2.3)
p2 = ggplot(data=Adynamic, aes(category, proportion, fill=conserve_staVScs)) + 
    geom_col(width=0.5,  alpha=0.5, position=position_dodge()) + 
    ylab('% of HIMs') + xlab('% of genes in A compartments') + 
    scale_y_continuous(expand=c(0, 0.5)) + 
    scale_fill_manual(name="",values = col_list, labels=rename_features('A', nl)) +
    theme(axis.title = element_text(size=10), axis.text=element_text(size=8), 
          axis.text.x=element_text(angle=25, hjust=1), axis.ticks.x=element_blank()) + 
    facet_grid(~cell, labeller = as_labeller(rename_features('A', nl))) + 
    theme( strip.background = element_blank(), strip.text.x = element_blank() ) + 
    geom_text(data=AdynamicPval, aes(x=x, y=y, label=pval), size=3, inherit.aes = F)  + 
    theme(legend.position='bottom', legend.box.margin=margin(t=-20, r=0, b=-10, l=0),
                                      legend.text=element_text(size=10))
print(p2)
dev.off()
# Apercent for constitutive and cell type-specific HIMs stops here 
feat = c('edge_density', 'tfgene_density', 'motif_density', 'repli_mean', 'repli_cv')
nrow=length(feat); ncol=5
#pdf(file='sup_fig/sup_stable_vs_cs_spatial.pdf', width=1.2 * ncol, height=1.8 * nrow)
png(file='sup_fig/sup_stable_vs_cs_spatial.png', width=1.2 * ncol, height=1.8 * nrow, units='in', res=300)
fi = 'conserve_staVScs'
Xf = X[, c(feat, 'cell', fi)]
p1 = boxdot2(df=Xf, fs='cell', fn=fi, g1='stable', g2='cs', ylabsize=1.1, pvalalter = F, psize=1, amaxy=1.05)
print(p1)
dev.off()
