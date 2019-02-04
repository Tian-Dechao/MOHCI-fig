rm(list=ls())
library(reshape2)
library(plyr)
library(ggplot2)
source('src/jaccard_index_sig.R')
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
hims_ji_expected = compute_ji_expected_ev()
hims_ji_fc = hims_ji[, c('jiTF', 'jiGene')] / hims_ji_expected[, c('jiExp_tf', 'jiExp_gene')]
hims_ji_fc = log2(hims_ji_fc)
colnames(hims_ji_fc) = c('jiTF_fc', 'jiGene_fc')
hims_ji_comb = data.frame(hims_ji, hims_ji_expected, hims_ji_fc)
head(hims_ji_comb)
write.table(hims_ji_comb, file='data/ji_fc.txt', col.names = T, row.names = F, sep='\t', quote=F)
#### stop here 
##################################################
##################################################
##################################################
ggplot(hims_ji_comb, aes(x=jiTF_fc, y=jiGene_fc)) + geom_point(size=0.1) +  geom_smooth(method='lm') + 
    facet_grid(cell1 ~ cell2)
fji = c('jiGene_fc', 'jiTF_fc')
Y = melt(hims_ji_comb, id.vars = c('cell1', 'cell2'), measure.vars = fji, variable.name = "feature", value.name="value")
Y[, 'cell1'] = factor(Y[, 'cell1'])
Y[, 'cell2'] = factor(Y[, 'cell2'])
#Y.mu = ddply(Y, 'feature', summarise, grp.mean=mean(value))
Y.me = ddply(Y, 'feature', summarise, grp.mean=median(value))
#pdf('main_fig/jaccard_index_hist.pdf', width=2, height=2)
hist.plot = ggplot(Y, aes(x=value, color=feature)) + 
    geom_histogram(aes(y=..density..), fill='white', alpha=0.2, position='identity') +
#    geom_density(alpha=0.2)+
#    facet_grid(feature~.) +
#    xlim(0, 0.5) + 
#    ylim(0, 10) +
    theme_minimal() +  theme_classic() + 
    theme(axis.text.y = element_text(color='black', size=7)) +
    theme(axis.text.x = element_text(color='black', size=7)) +
    theme(legend.position='bottom' ) +
    geom_vline(data=Y.me, aes(xintercept=grp.mean, color=feature), linetype='dashed', size=1)
print(hist.plot)
#dev.off()
