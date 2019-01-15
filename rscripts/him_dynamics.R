rm(list=ls())
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(plyr)
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
compare_2vects = function(x1, x2, paired=F){
     x1[is.infinite(x1)] = max(x1[is.finite(x1)])
     m1 = median(x1, na.rm=T)
     m2 = median(x2, na.rm=T)
     if(length(x2) > 1){
         x2[is.infinite(x2)] = max(x2[is.finite(x2)])
         pl = wilcox.test(x1, x2, alternative='l', paired=paired)$p.value
         pg = wilcox.test(x1, x2, alternative='g', paired=paired)$p.value
     } else {
        pl = wilcox.test(x1, mu=x2, alternative='l')$p.value
         pg = wilcox.test(x1, mu=x2, alternative='g')$p.value
     }
     if(pl <= pg){
         pl = format(pl, scientific=T, digits=3)
         p = c(pl, 'l')
     } else {
         pg = format(pg, scientific=T, digits=3)
         p = c(pg, 'g')
     }
     if(as.numeric(p[1]) <= 0.05){
         p = c(p, '***')
     } else {
         p = c(p, '')
     }
     p2 = p[1]
     if(as.numeric(p2) <= 2.22e-16){
         p2 = 'P<\n2.22e-16'
     } else{
         p2 = paste('P=\n', p2, sep='')
     }
     result = c(m1, m2, p, p2)
     return(result)
 }
X = read.table('data/him_dynamic_allinone.txt', header=T, stringsAsFactors = F, sep='\t')
fji = c('jiTF', 'jiGene')
####################
## him dynamics stratified by housekeeping genes and essential genes 
####################
Xsub = subset(X, select=c('himid1', 'cell1', 'himid2', 'cell2', fji))
df = read.table('data/venn.txt', header=T, row.names = 1, stringsAsFactors = F, sep='\t')
genes = rownames(df); genes = unlist(sapply(genes, function(z) strsplit(z, ';')))
gess = read.table('data/essential_combined.txt', header=F, stringsAsFactors = F)[, 1]
gess = gess[ gess %in% genes]
gess_k = read.table('data/essential_genes_k562.txt', header=F, stringsAsFactors = F)[, 1]
gess_k = gess_k[ gess_k %in% genes]
ghk = read.table('data/HK_genes.txt', header=F, stringsAsFactors = F)[, 1]
ghk = ghk[ghk %in% genes]
# step 1. extract essential genes or housekeeping genes in each HIMs
himsum = read.table('data/him_summary_allinone.txt', header=T, stringsAsFactors = F, row.names = 1, sep='\t')
himsum = subset(himsum, source=='him', select=c('genes'))
himsum_g = lapply(himsum[, 'genes'], function(z) unlist(strsplit(z, ',|;')))
himsum_ess = sapply(himsum_g, function(z) z[z %in% gess])
himsum_hk = sapply(himsum_g, function(z) z[z %in% ghk])
names(himsum_ess) = rownames(himsum)
names(himsum_hk) = rownames(himsum)
# step 2. compute the number of shared essential genes 
shared_ess = mapply(function(x, y) {sum(x %in% y)}, x=himsum_ess[Xsub[, 'himid1']], y=himsum_ess[Xsub[, 'himid2']])
shared_hk = mapply(function(x, y) {sum(x %in% y)}, x=himsum_hk[Xsub[, 'himid1']], y=himsum_hk[Xsub[, 'himid2']])
shared_ess[shared_ess>=5] = 5
shared_hk[shared_hk>=5] = 5
shared_hk = as.character(shared_hk)
unique(shared_hk)
Xsub = data.frame(Xsub, ness=shared_ess, nhk=shared_hk, stringsAsFactors = F)
Xsub$nhk = factor(Xsub$nhk, levels=as.character(0:5), labels=c(as.character(0:4), ">=5") )
# step 3. compute p values between groups
nhk_cat = unique(Xsub$nhk)
pval = c()
for(i in 1:(length(nhk_cat)-1)){
    tmp = compare_2vects(x1=Xsub[Xsub$nhk == nhk_cat[i], 'jiTF'], x2=Xsub[Xsub$nhk == nhk_cat[i+1], 'jiTF'], paired=F)
    pval = rbind(pval, tmp)
}
pval
ylim_hk = boxplot.stats(Xsub$jiTF)$stats[c(1, 5)]
ylim_hk[2] = ylim_hk[2] + 0.5 * diff(ylim_hk)
pdf('main_fig/jiTF_hk.pdf', width=2, height=3)
ggplot(Xsub, aes(x=nhk, y=jiTF)) + 
    geom_boxplot(outlier.shape = NA, fill=gg_color_hue(n=2)[2], color='grey') + 
    coord_cartesian(ylim = ylim_hk) + 
    annotate('text', x=1:5 + 0.5, y=ylim_hk[2] * (6:10) / 10, label=pval[, 6], size= 4 * 5 / 14) +
    scale_x_discrete(labels=c(0:4, expression("">=5))) +
    theme_classic() + 
    xlab('# HK genes shared between\ntwo HIMs from two cell types') + 
    ylab(expression(JI[TF])) + 
    theme(axis.title=element_text(size=8), axis.text=element_text(size=7)) 
dev.off()
ness_cat = unique(Xsub$ness)
pval_ess = c()
for(i in 1:(length(ness_cat)-1)){
    tmp = compare_2vects(x1=Xsub[Xsub$ness == ness_cat[i], 'jiTF'], x2=Xsub[Xsub$ness == ness_cat[i+1], 'jiTF'], paired=F)
    pval_ess = rbind(pval_ess, tmp)
}
Xsub$ness = factor(Xsub$ness, levels=as.character(0:5), labels=c(as.character(0:4), ">=5") )
ylim_ness = boxplot.stats(Xsub$jiTF)$stats[c(1, 5)]
ylim_ness[2] = ylim_ness[2] + 0.5 * diff(ylim_ness)
pdf('sup_fig/jiTF_ess.pdf', width=2, height=2.2)
ggplot(Xsub, aes(x=ness, y=jiTF)) + 
    geom_boxplot(outlier.shape = NA, fill=gg_color_hue(n=2)[2], color='grey') + 
    coord_cartesian(ylim = ylim_ness) + 
    annotate('text', x=1:5 + 0.5, y=ylim_ness[2] * (5:9) / 10, label=pval_ess[, 6], size=5 * 5 / 14) +
    scale_x_discrete(labels=c(0:4, expression("">=5))) +
    theme_classic() + 
    xlab('# essential genes shared between\ntwo HIMs from two cell types') + 
    ylab(expression(JI[TF])) + 
    theme(axis.title=element_text(size=8), axis.text=element_text(size=7)) 
dev.off()
#
## nucleous  and speckle connections
hubtype = read.table('nucleoli/sprite_hub_annotation_gm12878.txt', header=T, sep='\t', stringsAsFactors = F)
hubtype$Group.1 = gsub('HIM', 'gm12878', hubtype$Group.1)
tsa.gm12878 = read.table('nucleoli/gm12878_tsa.txt', header=T, sep='\t', stringsAsFactors = F)
tsa.gm12878 = tsa.gm12878[tsa.gm12878$source == 'him', ]
row.names(tsa.gm12878) = paste('gm12878', tsa.gm12878$index, sep='_')
hubtype$tsa = tsa.gm12878[hubtype[, 'Group.1'], 'SON_Sucrose_gene_mean']
ggplot(hubtype, aes(x = x.nucleous, y=tsa)) + geom_point() + geom_smooth(method='loess')
# test stop here
## jiGene and jiTF connections  
ggplot(X, aes(x=jiGene, y=jiTF)) + geom_point() + geom_smooth(method='loess') + 
    ylim(0, 0.4) + 
    facet_grid(cell1~cell2)
## test stop here
fhic = c('hic_density_ratio', 'A_percent_ratio', 'tad_num_ratio', 'loop_num_ratio', 'repli_mean_ratio')
fgrn = c('tfgene_density_ratio', 'motif_density_ratio')
fexp = c('gene_expression_mean_ratio', 'tf_expression_mean_ratio', 'all_expression_mean_ratio') 
fedge = c('hic_gain', 'hic_lost', 'grn_gain', 'grn_lost')
fall = c(fji, fhic, fgrn, fexp)
for(i in 1:ncol(X)){
    X[is.infinite(X[, i]), i] = max(X[is.finite(X[, i]), i], na.rm = T)
}
X = X[!apply(X[, fall], 1, function(z) any(is.na(z))), ]
# GO plot
himtype = c(rep('Constitutive genes', 5), 
            rep('GM12878 specific genes', 2), 
            rep('NEHK specific genes',3))
himtype = factor(himtype, levels=c('Constitutive genes','GM12878 specific genes', 'NEHK specific genes'))
goterm = c('chromosome organization', 'nuclesome organization',  'DNA packaging', 
           'regulation of gene expression, epigenetic', 'RNA splicing',
           'regulation of lymphocyte activation', 'regulation of T cell activation',
           'keratinocyte differentiation', 'keratinization', 'epidermis development')
goterm = factor(goterm, levels=rev(goterm))
pvals = c(5.7e-19, 4.2e-12, 5.2e-10, 3.4e-12, 3.4e-9,
          6.5e-8, 7.4e-8,
          2.5e-16, 2.0e-15, 2.6e-14)
pvals = -1* log10(pvals)
godf = data.frame(HIMtype=himtype, GOterm=goterm, Pvals=pvals, stringsAsFactors = F)
pdf('main_fig/go_pval_barplot.pdf', width=3.5, height=2.5)
ggplot(godf, aes(GOterm, Pvals, colour=HIMtype, fill=HIMtype)) + geom_col(width=0.6) + 
    xlab('') +  ylab('-log10(P value)') + 
    scale_y_continuous(expand=c(0, 0)) + 
    theme_minimal() +  theme_classic()  + 
    theme(legend.position = 'bottom') + 
    theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size=8)) + 
    coord_flip() 
dev.off() 
# create the table in latex for all the GO terms 
gotab = read.csv('GO_term_result_V6.csv', stringsAsFactors = F)
gotab = gotab[1:(nrow(gotab)-1),]
str(gotab)
gotab$P.value = format(gotab$P.value, scientific=T, digits=3)
library(xtable)
print(xtable(gotab), include.rownames = F)


# density plot
fji = c('jiGene', 'jiTF')
Y = melt(X, id.vars = c('cell1', 'cell2'), measure.vars = fji, variable.name = "feature", value.name="value")
Y[, 'cell1'] = factor(Y[, 'cell1'])
Y[, 'cell2'] = factor(Y[, 'cell2'])
#Y.mu = ddply(Y, 'feature', summarise, grp.mean=mean(value))
Y.me = ddply(Y, 'feature', summarise, grp.mean=median(value))
head(Y)
pdf('main_fig/jaccard_index_hist.pdf', width=2, height=2)
hist.plot = ggplot(Y, aes(x=value, color=feature)) + 
    geom_histogram(aes(y=..density..), fill='white', alpha=0.2, position='identity') +
#    geom_density(alpha=0.2)+
#    facet_grid(feature~.) +
    xlim(0, 0.5) + 
    ylim(0, 10) +
    theme_minimal() +  theme_classic() + 
    theme(axis.text.y = element_text(color='black', size=7)) +
    theme(axis.text.x = element_text(color='black', size=7)) +
    theme(legend.position='bottom' ) +
    geom_vline(data=Y.me, aes(xintercept=grp.mean, color=feature), linetype='dashed', size=1)
print(hist.plot)
dev.off()

    geom_vline(aes(xintercept=mean(value)), color='blue', linetype='dashed' , size=1.5) +
pdf('main_fig/jaccard_index.pdf', width=3, height=1)
pcomp = ggplot(Y, aes(factor(feature), value)) + 
    geom_boxplot(outlier.shape=NA, fill='white', colour='blue') + 
    coord_flip() + 
    ylim(0, 0.4) + 
    theme(axis.text.y = element_text(color='black', size=7)) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pcomp
dev.off()
pdf('main_fig/jaccard_index_vertical.pdf', width=1.5, height=3)
pcomp = ggplot(Y, aes(factor(feature), value)) + 
    geom_boxplot(outlier.shape=NA, fill='white', colour='blue') + 
    ylim(0, 0.4) + 
    theme(axis.text.y = element_text(color='black', size=7)) +
    theme(axis.text.x = element_text(color='black', size=7)) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pcomp
dev.off()
# histogram plot
X.

# changes 
#Z = melt(X, id.vars = c('cell1', 'cell2'), measure.vars = c('gene_expression_mean_ratio', 'hic_density_ratio', 'tfgene_density_ratio', 'motif_density_ratio'), variable.name = 'feature', value.name='value')
#Z = melt(X, id.vars = c('cell1', 'cell2'), measure.vars = c('hic_density_ratio', 'tfgene_density_ratio'), variable.name = 'feature', value.name='value')
Z = melt(X, id.vars = c('cell1', 'cell2'), measure.vars = c('hic_density_ratio', 'motif_density_ratio'), variable.name = 'feature', value.name='value')
Z[, 'value'] = log2(Z[, 'value'])
Z[, 'cell1'] = factor(Z[, 'cell1'])
Z[, 'cell2'] = factor(Z[, 'cell2'])
Z.me = ddply(Z, 'feature', summarise, grp.mean=median(value))
range(X$motif_density_ratio, finite=T)
fivenum(X$motif_density_ratio)
m1=median(X$motif_density_ratio)
log2(m1)
m2 = median(X$hic_density_ratio)
log2(m2)
head(Z)
pdf('main_fig/density_fc_hist.pdf', width=2, height=2)
hist.plot = ggplot(Z, aes(x=value, color=feature)) + 
    geom_histogram(aes(y=..density..), fill='white', alpha=0.2, position='identity') +
    #xlim(-1, 2.5) + 
    xlim(-1, 5.5) + 
    theme_minimal() +  theme_classic() + 
    theme(axis.text.y = element_text(color='black', size=7)) +
    theme(axis.text.x = element_text(color='black', size=7)) +
    theme(legend.position='bottom' ) +
    geom_vline(data=Z.me, aes(xintercept=grp.mean, color=feature), linetype='dashed', size=1)
print(hist.plot)
dev.off()
pdf('main_fig/density_fc.pdf', width=3, height=1)
pcomp = ggplot(Z, aes(factor(feature), value)) + 
    geom_boxplot(outlier.shape=NA, fill='white', colour='blue') + 
    coord_flip() + 
    ylim(-1, 3) + 
    theme(axis.text.y = element_text(color='black', size=7)) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pcomp
dev.off()
pdf('main_fig/density_fc_vertical.pdf', width=1.5, height=3)
pcomp = ggplot(Z, aes(factor(feature), value)) + 
    geom_boxplot(outlier.shape=NA, fill='white', colour='blue') + 
    ylim(-1, 3) + 
    theme(axis.text.y = element_text(color='black', size=7)) +
    theme(axis.text.x = element_text(color='black', size=7)) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
pcomp
dev.off()
# compute P-value
compare_2vects(X$hic_density_ratio, X$tfgene_density_ratio)
# create color scale 
z = runif(98, -2.3, 2.3)
z = c(z, -2.3, 2.3)
df = data.frame(x = runif(100), y=runif(100), z=z)
pdf('main_fig/color_key.pdf')
ggplot(df, aes(x, y)) + 
    geom_point(aes(colour=z)) + 
    scale_colour_gradient2(low='#3F4498', mid='white', high = 'red') + 
    theme(legend.position = 'bottom')
dev.off()
png('main_fig/color_key.png')
ggplot(df, aes(x, y)) + 
    geom_point(aes(colour=z)) + 
    scale_colour_gradient2(low='#3F4498', mid='white', high = 'red') + 
    theme(legend.position = 'bottom')
dev.off()

pcomp = ggplot(subset(Z, feature!='motif_density_ratio' & feature!='gene_expression_mean_ratio'), aes(feature, value)) + 
    geom_boxplot(outlier.shape = NA) + 
    ylim(-1, 3)
pcomp
pcomp = ggplot(subset(Z, feature!='motif_density_ratio' & feature!='gene_expression_mean_ratio'), aes(feature, value)) + 
    geom_boxplot(outlier.shape = NA) + 
    ylim(-1, 3) +
    facet_grid(cell1 ~ cell2)
pcomp
# rewire: gain or lost
W = melt(X, id.vars = c('cell1', 'cell2'), measure.vars = fedge, variable.name = 'feature', value.name='value')
W[, 'cell1'] = factor(W[, 'cell1'])
W[, 'cell2'] = factor(W[, 'cell2'])
pcomp = ggplot(W, aes(feature, value)) +
    geom_boxplot(outlier.shape = NA) + 
    facet_grid(cell1~cell2, switch='y') + 
    ylim(0, 1)
pcomp
# stop here

pcomp = ggplot(Y, aes(x=value, colour=feature, fill=feature)) + 
    geom_density(alpha=0.1) + 
    xlim(0, 0.3) + 
    facet_grid(~cell1+cell2)
pcomp
pcomp = ggplot(Y, aes(x=value, colour=feature, fill=feature)) + 
    geom_density(alpha=0.1) + 
    xlim(-0.3, 0.3)  
pcomp
pcomp = ggplot(subset(X, jiGene>=0.05), aes(x=jiTF, y=jiGene)) + 
    geom_hex(bins=100) + 
    scale_color_gradient2(midpoint=0.15, low="blue", mid="white", high="red") + 
    geom_smooth()
pcomp
# plot changes for each feature 
pdf('him_dynamic.pdf', width=6, height=7)
col = brewer.pal(n=5, name='Dark2')
par(mfrow=c(1,1))
layout(matrix(c(1, 2, 3, 3), nrow=2, byrow=T))
# jigene and jiTF
plot(-1, -1, xlim=c(0, 2), ylim=c(0, 0.5), xlab='', ylab='Jaccard index', las=1,  xaxt = 'n', tck=-0.02, mgp=c(1.5, 0.5,0) )
boxplot(X[, 'jiGene'], col=col[4], at=0.5, add=T, box.wex=0.75, outline = F, yaxt='n')
boxplot(X[, 'jiTF'], col=col[5], at=1.5, add=T, box.wex=0.75, outline=F, yaxt='n')
#vioplot(X[, 'jiGene'], col=col[4], at=0.5, add=T, wex=0.75, border=NA)
#vioplot(X[, 'jiTF'], col=col[5], at=1.5, add=T, wex=0.75, border=NA)
axis(1, at=c(0.5, 1.5), labels = c('Gene', 'TF'), tck=-0.02, mgp=c(1.5, 0.25,0))
compare_2vects(X[, 'jiGene'], X[, 'jiTF'], paired = T)
# hic_lost and grn_lost
plot(-1, -1, xlim=c(0, 2), ylim=c(0, 100), xlab='', ylab='% edges lost', las=1,  xaxt = 'n', tck=-0.02, mgp=c(1.5, 0.5,0) )
boxplot(X[, 'hic_lost'] * 100, col=col[4], at=0.5, add=T, box.wex=0.75, outline = F, yaxt='n')
boxplot(X[, 'grn_lost'] * 100, col=col[5], at=1.5, add=T, box.wex=0.75, outline = F, yaxt='n')
#vioplot(X[, 'hic_lost'], col=col[4], at=0.5, add=T, wex=0.75, border=NA, pchMed = 2)
#vioplot(X[, 'grn_lost'], col=col[5], at=1.5, add=T, wex=0.75, border=NA, pchMed = 2)
axis(1, at=c(0.5, 1.5), labels = c('Hi-C network', 'GRN'), tck=-0.02, mgp=c(1.5, 0.25,0))
compare_2vects(X[, 'hic_lost'], X[, 'grn_lost'], paired=T)
## hic_lost and grn_lost
#plot(-1, -1, xlim=c(0, 2), ylim=c(0, 1), xlab='', ylab='% edges gained', las=1,  xaxt = 'n', tck=-0.02, mgp=c(1.5, 0.25,0) )
#vioplot(X[, 'hic_gain'], col=col[4], at=0.5, add=T, wex=0.75, border=NA)
#vioplot(X[, 'grn_gain'], col=col[5], at=1.5, add=T, wex=0.75, border=NA)
#axis(1, at=c(0.5, 1.5), labels = c('Hi-C network', 'GRN'))
#compare_2vects(X[, 'hic_gain'], X[, 'grn_gain'], paired=T)
# gene expression

plot(-1, -1, xlim=c(0, 4), ylim=c(-3, 15), xlab='', ylab='log2(ratio)', las=1,  xaxt = 'n', tck=-0.02, mgp=c(1.5, 0.5,0) )
boxplot(log2(X[, 'gene_expression_mean_ratio']), col=col[3], at=0.5, add=T, box.wex=0.75, outline=F, yaxt='n')
boxplot(log2(X[, 'hic_density_ratio']), col=col[4], at=1.5, add=T, box.wex=0.75, outline=F, yaxt='n')
boxplot(log2(X[, 'tfgene_density_ratio']), col=col[5], at=2.5, add=T, box.wex=0.75, outline=F, yaxt='n')
boxplot(log2(X[, 'motif_density_ratio']), col=col[2], at=3.5, add=T, box.wex=0.75, outline=F, yaxt='n')
#vioplot(log2(X[, 'gene_expression_mean_ratio']), col=col[3], at=0.5, add=T, wex=0.75, border=NA)
#vioplot(log2(X[, 'hic_density_ratio']), col=col[4], at=1.5, add=T, wex=0.75, border=NA)
#vioplot(log2(X[, 'tfgene_density_ratio']), col=col[5], at=2.5, add=T, wex=0.75, border=NA)
#vioplot(log2(X[, 'motif_density_ratio']), col=col[2], at=3.5, add=T, wex=0.75, border=NA)
abline(h=1, lwd=2, lty='dashed')
axis(1, at=c(0.5, 1.5, 2.5, 3.5), labels = c('Gene expression','Hi-C', 'GRN','Motif'), tck=-0.02, mgp=c(1.5, 0.25,0))
dev.off()

compare_2vects((X[, 'gene_expression_mean_ratio']), (X[, 'hic_density_ratio']), paired=T)
compare_2vects((X[, 'hic_density_ratio']), (X[, 'tfgene_density_ratio']), paired=T)
compare_2vects((X[, 'tfgene_density_ratio']), (X[, 'motif_density_ratio']), paired=T)


plot(X[, 'jiTF'], log(X[, 'motif_density_ratio']))
plot(X[, 'jiGene'], log(X[, 'hic_density_ratio']))
y = cut(X$jiTF, breaks=c(0, 0.05 , 0.1, 0.3, 1), include.lowest = T)
boxplot(log(X$motif_density_ratio)~y)
cor(1:4, by(X$motif_density_ratio, y, median))
 # correlations
fivenum(X$jiGene)
Xsub = subset(X, jiGene>0.05)
c = c()
#for(i in c('jiTF', 'jiGene')){
 for(i in c('jiTF', 'jiGene', 'gene_expression_mean_ratio', 'tf_expression_mean_ratio', 'gene_expression_cv_ratio', 'tf_expression_cv_ratio')){
     for(j in c(fhic, fgrn, fexp, fedge)){
         rho = cor(Xsub[, i], Xsub[, j], use='pairwise.complete.obs', method='spearman')
         rho = c(i, j, round(rho, 3))
         p = cor.test(Xsub[, i], Xsub[, j], method='spearman', alternative='t')$p.value
         p = format(p, scientific=T, digits=3)
         c = rbind(c, c(rho, p))
     }
 }
 print(c)
 plot(Xsub$jiGene, Xsub$hic_lost)
 df = data.frame(motif_density_ratio=log2(Xsub$motif_density_ratio),hic_edge_gain=Xsub$hic_gain, hic_edge_lost=(Xsub$hic_lost),hic_edge_density_ratio=log2(Xsub$hic_density_ratio), gene_expression_mean_ratio=log2(Xsub$gene_expression_mean_ratio)) 
 plot(log2(Xsub$gene_expression_mean_ratio), log2(Xsub$hic_density_ratio), pch=16, xlim=c(-5, 5))
 d = ggplot(df, aes(x=hic_edge_density_ratio,y=gene_expression_mean_ratio)) + geom_hex(bins=100)
 d = ggplot(df, aes(x=hic_edge_lost,y=gene_expression_mean_ratio)) + geom_hex(bins=100) + geom_smooth(method='lm')
 d = ggplot(df, aes(x=hic_edge_gain,y=gene_expression_mean_ratio)) + geom_hex(bins=100) + geom_smooth(method='lm')
 d = ggplot(df, aes(x=motif_density_ratio,y=gene_expression_mean_ratio)) + geom_hex(bins=100) + geom_smooth(method='loess', colour='red', size=3)
 d
 cn_dep = 'gene_expression_mean_ratio'
 cn_ind = c('motif_density_ratio', 'hic_edge_gain', 'hic_edge_lost', 'hic_edge_density_ratio')
 # too slow
formu = paste(cn_dep, '~', paste(cn_ind, collapse='+'), collapse='')
formu = as.formula(formu)
bw=npregbw(formu, data=df)
model = npreg(bw)
plot(model)
library(randomForest)
library(doParallel)

cn_dep = 'gene_expression_mean_ratio'
cn_dep = 'jiGene'
cn_ind= c(fedge, fgrn, 'hic_density_ratio')
formu = paste(cn_dep, '~', paste(cn_ind, collapse='+'), collapse='')
formu = as.formula(formu)
dfrf = X[, c(cn_dep, cn_ind)]
rf = randomForest(formu, data=dfrf, ntree=200,  importance=T, proximity=F, keep.forest=F, do.trace=F )
varImpPlot(rf)
rf$importance[, '%IncMSE']
fivenum(rf$rsq)
fivenum(mean(rf$mse))
# do not run
rf = foreach(ntree_p=rep(50,4), .combine=combine, .multicombine=TRUE, .packages='randomForest') %dopar% {
    randomForest(formu, data=dfrf, ntree=ntree_p, importance=T, keep.forest=F, proximity=F, do.trace=F)
}
varImpPlot(rf)


 # five number summary on the gene expresion
apply(X[, fexp], 2, fivenum)
apply(X[, fexp], 2, quantile, probs=c(0.75, 0.9, 0.95))
boxplot(gene_expression_mean_ratio ~ conserve_All, X, ylim=c(0, 5))
abline(h=2)
# what does this mean? df hims have higher proportion of HIMs that have >2-fold expression change
Xuniq = X[!duplicated(X[, 'himid1']), ]
aggregate(cbind(gene_expression_mean_ratio, tf_expression_mean_ratio)~conserve_All, data=Xuniq, function(z){
    c(sum(z>=2),length(z), sum(z>=2)/length(z))
    } )

boxplot(gene_expression_mean ~ conserve_All + cell, Xuniq, ylim=c(0, 5))

 
# analysis interaction dynamics
apply(X[, fedge], 2, fivenum)
# gain
compare_2vects(X[, fedge[1]], X[, fedge[3]], paired=T)
# lost
compare_2vects(X[, fedge[2]], X[, fedge[4]], paired=T)
# gain vs lost
compare_2vects(X[, fedge[1]], X[, fedge[2]], paired=T)
compare_2vects(X[, fedge[3]], X[, fedge[4]], paired=T)
# lost exclude GM
indgm = X$cell1!='gm12878' & X$cell2!='gm12878'
compare_2vects(X[indgm, fedge[2]], X[indgm, fedge[4]], paired=T)
par(mfrow=c(1,1))
layout(matrix(c(1:3), nrow=1))
boxplot(X[, fedge[c(1, 3, 2, 4)]])
boxplot(X[, fedge[1]] ~ X[, 'cstypepair'], at = 1:4 - 0.25, col='red', boxwex=0.25)
boxplot(X[, fedge[2]] ~ X[, 'cstypepair'], at = 1:4 + 0.25, add=TRUE, col='blue', axes=F, boxwex=0.25)
legend('topright', fedge[1:2], fill=c('red', 'blue'))
boxplot(X[, fedge[3]] ~ X[, 'cstypepair'], at = 1:4 - 0.25, col='red', boxwex=0.25)
boxplot(X[, fedge[4]] ~ X[, 'cstypepair'], at = 1:4 + 0.25, add=TRUE, col='blue', axes=F, boxwex=0.25)
legend('topright', fedge[3:4], fill=c('red', 'blue'))


# misc
x = c(38,  349,  292,  17)
x1 = c(38,  349,  246,  17)
y = c(19,     327,      15,     288)
sum(x1) - sum(y)
hela1 = c(59,  457,  274,  16)
hela2 = c(39,     443,      15,     309)
sum(hela1) - sum(hela2)
# test on expression data
x = read.table('~/Downloads/ENCFF051IJZ.tsv ', header = T, stringsAsFactors = F, sep='\t')
