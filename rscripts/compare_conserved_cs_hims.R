rm(list=ls())
source('src/boxdot.R')
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
#### cell identify 
# comment from Zhijun Duan
# super-enhancer related 
#feat = c('SE_num_100kb', 'SE_num_100kb_per_Mb', 'SE_density_100kb','SE_num_500kb', 'SE_num_500kb_per_Mb', 'SE_density_500kb', 'SE_num_1Mb', 'SE_num_1Mb_per_Mb','SE_density_1Mb')
feat = c('SE_num_20kb_per_Mb', 'SE_num_50kb_per_Mb', 'SE_num_100kb_per_Mb','SE_num_500kb_per_Mb','SE_num_1Mb_per_Mb')
nrow=length(feat); ncol=5
pdf(file='sup_fig/sup_stable_vs_cs_se.pdf', width=1.2 * ncol, height=1.8 * nrow)
#png(file='sup_fig/sup_stable_vs_cs_spatial.png', width=1.2 * ncol, height=1.8 * nrow, units='in', res=300)
fi = 'conserve_staVScs'
Xf = X[, c(feat, 'cell', fi)]
p1 = boxdot2(df=Xf, fs='cell', fn=fi, g1='stable', g2='cs', ylabsize=1.1, pvalalter = F, psize=1, amaxy=1.05)
print(p1)
dev.off()
# enseential genes 
feat = c('ess_all_gene_p', 'ess_k562_gene_p')
nrow=length(feat); ncol=5
pdf(file='sup_fig/sup_stable_vs_cs_ess.pdf', width=1.2 * ncol, height=1.8 * nrow)
fi = 'conserve_staVScs'
Xf = X[, c(feat, 'cell', fi)]
p1 = boxdot2(df=Xf, fs='cell', fn=fi, g1='stable', g2='cs', ylabsize=1.1, pvalalter = F, psize=1, amaxy=1.05)
print(p1)
dev.off()
# master TFs 
feat = c('master_in', 'master_in_p')
nrow=length(feat); ncol=5
pdf(file='sup_fig/sup_stable_vs_cs_master.pdf', width=1.2 * ncol, height=1.8 * nrow)
fi = 'conserve_staVScs'
Xf = X[, c(feat, 'cell', fi)]
p1 = boxdot2(df=Xf, fs='cell', fn=fi, g1='stable', g2='cs', ylabsize=1.1, pvalalter = F, psize=1, amaxy=1.05)
print(p1)
dev.off()
# gene centered comparison 
# load gene info in G
source('src/load_gene_info.R')
# create gene list in conserved and cell type specific HIMs 
source('src/him_short2long.R')
CN = c('cell', 'conserve_staVScs', 'genes')
hims_long = short2long(X[X$source == 'him', ], CN=CN, cn='genes')
# merge G and hims_long 
rownames(G) = paste(as.character(G$cell), as.character(G$gene_name), sep='_')
rownames(hims_long) = paste(hims_long$cell, hims_long$genes, sep='_')
####### double check
hims_long[c('gm12878_NEDD4', 'gm12878_NEDD4L', 'www', 'gm12878_NEDD4'), ]
# try hims_long
rn1 = rownames(hims_long); rn2= rownames(G); rn1 = rn1[rn1 %in% rn2]
hims_long = hims_long[rn1, ]
tmp = matrix(, nrow=nrow(G), ncol=ncol(hims_long)); tmp = as.data.frame(tmp); 
dimnames(tmp) = list(rownames(G), colnames(hims_long))
tmp[rownames(hims_long), ] = hims_long
G = cbind(G, tmp)
# double check gene assignment
table(G$HIM, G$conserve_staVScs)
ind = G$HIM == 'No' & G$conserve_staVScs == 'cs'
G[is.na(G$conserve_staVScs), 'conserve_staVScs'] = 'NoHIM'
G$conserve_staVScs = factor(G$conserve_staVScs, levels=c('stable', 'cs', 'NoHIM'))
#### visualize expression 
#col_list = c('stable' = "#2b83ba", 'cs'='red', 'No' = "#d7191c")
source('src/color.R')
col_list = gg_color_hue(n=3); names(col_list) = c('stable', 'cs', 'NoHIM')
feat = c('expression')
ncol=5; nrow=length(feat)
pdf(file='sup_fig/sup_stable_vs_cs_vs_not_assigned_expression.pdf', width=6.65, height=2 * nrow)
fi = 'conserve_staVScs'
Gf = G[, c(feat, fi, 'cell')]
ylim = boxplot.stats(Gf[, feat])$stats[c(1, 5)]; ylim[2] = 1.05 * ylim[2]
ggplot(Gf, aes(x=conserve_staVScs, y=expression, fill=conserve_staVScs)) + geom_boxplot(outlier.shape = NA) + 
#    geom_violin() + 
    ylab('Expression') +  xlab('') + 
    scale_y_continuous(limits = ylim) + 
    scale_fill_manual(values=col_list, labels=rename_features(names(col_list),nl), name='Gene assignment to HIMs') + 
    facet_grid(.~cell, labeller=as_labeller(rename_features('', nl))) + 
    theme(legend.position = 'bottom', axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank())
dev.off()
# compute the pvals 
cells = unique(G$cell)
pval = c()
for(ic in cells){
    Gsub = G[G$cell == ic, ]
    x = Gsub[Gsub$conserve_staVScs == 'stable', 'expression']
    y = Gsub[Gsub$conserve_staVScs == 'cs', 'expression']
    z = Gsub[Gsub$conserve_staVScs == 'NoHIM', 'expression']
    res1 = compare2vect(x=x, y=z)
    res2 = compare2vect(x=y, y=z)
    res = rbind(res1, res2)
    res = cbind(ic, res)
    pval = rbind(pval, res)
}
pval
max(pval[, 'pval'])
### stop here 
### need to redraw barplot
feat = c('SE_1000K')
ncol=5; nrow=length(feat)
pdf(file='sup_fig/sup_stable_vs_cs_vs_not_assigned_SE.pdf', width=6.65, height=2 * nrow)
fi = 'conserve_staVScs'
Gf = G[, c(feat, fi, 'cell')]
ylim = boxplot.stats(Gf[, feat])$stats[c(1, 5)]; ylim[2] = 1.05 * ylim[2]
ggplot(Gf, aes(x=SE_1000K, fill=conserve_staVScs)) + 
    geom_bar( position='dodge')  
    #geom_boxplot(outlier.shape = NA) + 
    ylab('# SE within 1Mb') +  xlab('') + 
    scale_y_continuous(limits = ylim) + 
    scale_fill_manual(values=col_list, labels=rename_features(names(col_list),nl), name='Gene assignment to HIMs') + 
    facet_grid(.~cell, labeller=as_labeller(rename_features('', nl))) + 
    theme(legend.position = 'bottom', axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x = element_blank())
dev.off()

