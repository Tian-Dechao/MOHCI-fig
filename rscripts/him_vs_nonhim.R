rm(list=ls())
source('src/boxdot.R')
source('src/load_him_summary_allinone.R')
source('src/load_gene_info.R')
##### check the cell type specific HIMs with cstyp 0, i.e., HIMs with majority of their genes only assigned to one cell type
## This part is optional
himcs0 = read.table('data/him_cell_type_specific_cstype_0.txt', header=T, sep='\t', stringsAsFactors = F)
himcs0 = X[himcs0$himid, ]
himcs0[, fchara]
himcs0[, c('gene_num', 'TF_number', fessg)]
# cell type specific genes are mostly in NHEK cstype 0 HIMs
himcs0[, c('gene_num', 'TF_number', fcsg)]
himcs0[, c('gene_num', 'TF_number', fexpg)]
# K562_128 has one master TF, and the genes are at most 203 away
himcs0[, c('gene_num', 'TF_number', fmas)]
# gm12878_63  loop_num_ratio = 3, tad_num_ratio=1; but only has 4 genes 
# nhek_107 loop_num_ratio=2, tad_num_ratio=0.5; gene_expression_mean_ratio=4.6, 
himcs0[, c('gene_num', 'TF_number', featratio)]
himcs0['nhek_107', c('gene_num', 'TF_number', featratio)]
## stop here

##### TSA 
Xk = X[X$cell == 'k562' & X$source == 'him', ]
Xk2 = X[X$cell == 'k562', ]
colnames(tsa)
##### for the main text
# proportion of HIMs that are close to speckle (mean tsa score >= 0.284)
n1 = sum(Xk$SON_Sucrose_gene_mean >= 0.284); n2 = nrow(Xk)
n1 / n2 * 100
# A_percent = 100
X2 = X[X$source == 'him', ]
aggregate(x=X2$A_percent, by=list(X2$cell), function(z) sum(z ==100)/length(z))
aggregate(x=X2$A_percent, by=list(X2$cell), function(z) length(z))
##### main text stop here 

##### find some hims with high lamina in K562; picked one to show in the main figure 
fivenum(Xk$LaminB_gene_mean)
# there is no him with positive lamina B and positive SON TSA-seq
indlamina = Xk$LaminB_gene_mean > 0.5
himlam =  Xk[indlamina, ]
himlam$LaminB_gene_mean
# 
#himlam[, c('gene_num', 'TF_number', fchara)]
#himlam[, c('gene_num', 'TF_number', fessg)]
#himlam[, c('gene_num', 'TF_number', fexpg)]
# multiple HIMs has master TFs
#himlam[, c('gene_num', 'TF_number', fmas)]
#himlam[, c('gene_num', 'TF_number', featratio)]
#himlam[, c('gene_num', 'TF_number', fcsg)]
# interesting featuers: csgp, A_percent, repli_mean, master_in, ess_k562_gene_p
fintersted = c('gene_num', 'TF_number', 'edge_density', 'tfgene_density', 'motif_density',
               'A_percent', 'repli_mean', 'master_in', 
               'ess_k562_gene_p', 'csgp', 'gene_expression_mean','gene_distance',  
               'LaminB_gene_mean' , 'SON_Sucrose_gene_mean')
himsub = himlam[, fintersted]; himsub = himsub[order(himsub$LaminB_gene_mean), ]
himsub = himsub[himsub$motif_density > 0.1, ]
himsub = himsub[himsub$repli_mean < 60, ]
coor_541 = himsub['k562_541', c('SON_Sucrose_gene_mean', 'LaminB_gene_mean')]
himsub['k562_541', c('gene_distance')]
(17798291-24581956) / 10^6
himsub['k562_541', c('gene_num','TF_number')]
himsub['k562_541', c('edge_density', 'tfgene_density', 'motif_density')]
# check K562_123, 541, which has essential genes and are highly expressed
# gene RPL15 in him_541 is a K562 essential genes
# gene SATB1 related to T-cell by genecards
# TF CBX5 in heterchromatin domain 
# TF CDX1 adult acute lymphocytic leukemia
# TF HOXA10 regulation of hematopoietic lineage commitment.
# TF Diseases associated with HOXA9 include Myeloid Leukemia

#### TSA related main figure and supplementary figures
source("src/tsa_related.R")
##### HIMs ~ A compartments 
source("src/compartment_related.R")
##### HIMSs ~ repli-seq on replication timing. 
source('src/replication_timing_related.R')
##### HIMSs ~ essential genes . 
source('src/essential_gene_related.R')
##### HIMSs ~ gene expression . 
source('src/gene_expression_related.R')



  

# 
# Main figures start here
#feat = c('A_percent', 'repli_mean', 'repli_cv', 'SON_Sucrose_gene_mean')
#feat = c('A_percent', 'repli_mean', 'SON_Sucrose_gene_mean')
#ncol=length(feat); nrow=ceiling(length(feat)/ncol)
#pdf(file='main_fig/him_vs_nonhim_spatial.pdf', width=1.5 * ncol, height=2.5 * nrow)
#fi = 'source'
#Xf = X[, c(feat, fi)]
#i=4
#Xf_c = Xf[X[, 'cell'] == cells[i], ]
#print(boxdot(df=Xf_c, fn=fi, g1='him', g2='merged non-him', cell=cells[i], fnrow=NULL, fncol=ncol, ylabsize=1.1, pvalalter = F))
#dev.off()
# mhim vs mnon-him in functional properties
#feat = c('ess_k562_gene_p', 'gene_expression_mean', 'SE_num_50kb_per_Mb', 'cosmic_snp_frequency')
#ncol=length(feat); nrow=ceiling(length(feat)/ncol)
#pdf(file='main_fig/him_vs_nonhim_functional.pdf', width=1.5 * ncol, height=3 * nrow)
#fi = 'source'
#Xf = X[, c(feat, fi)]
#Xf = Xf[Xf$SE_num_50kb_per_Mb <=10, ]
## k562 for essential genes and higher number of super-enhancers
#i=4
#Xf_c = Xf[X[, 'cell'] == cells[i], ]
##print(boxdot(df=Xf_c, fn=fi, g1='merged him', g2='merged non-him', cell=cells[i], fnrow=NULL, fncol=ncol, ylabsize=1.1, pvalalter = F, psize=2))
#print(boxdot(df=Xf_c, fn=fi, g1='him', g2='merged non-him', cell=cells[i], fnrow=NULL, fncol=ncol, ylabsize=1.1, pvalalter = F, psize=2))
#dev.off()
# try density plot
#ggplot(subset(Xf_c, source%in% c('him', 'merged non-him')), aes(SE_num_50kb_per_Mb, color=source)) + 
#     geom_density()  + xlim(0, 2)
pdf(file='main_fig/mhim_vs_nonhim_functional.pdf', width=1.5 * ncol, height=3 * nrow)
fi = 'source'
Xf = X[, c(feat, fi)]
# k562 for essential genes and higher number of super-enhancers
i=4
Xf_c = Xf[X[, 'cell'] == cells[i], ]
print(boxdot(df=Xf_c, fn=fi, g1='merged him', g2='merged non-him', cell=cells[i], fnrow=NULL, fncol=ncol, ylabsize=1.1, pvalalter = F, psize=2))
dev.off()
feat = c('SE_num_50kb_per_Mb', 'gene_expression_mean', 'gene_expression_cv')
ncol=length(feat); nrow=ceiling(length(feat)/ncol)
pdf(file='main_fig/mhim_vs_nonhim_functional_slides.pdf', width=1.5 * ncol, height=3 * nrow)
fi = 'source'
Xf = X[, c(feat, fi)]
# k562 for essential genes and higher number of super-enhancers
i=4
Xf_c = Xf[X[, 'cell'] == cells[i], ]
print(boxdot(df=Xf_c, fn=fi, g1='merged him', g2='merged non-him', cell=cells[i], fnrow=NULL, fncol=ncol, ylabsize=1.1, pvalalter = F, psize=2))
dev.off()
# ppi
#source('~/Dropbox/r_codes/boxdot.R')
boxdot_ppi(df=X, fm='main_fig/him_ppi.pdf', fs='sup_fig/sup_him_ppi.pdf')
# compute P values 
Pval = matrix(0,nrow=5, ncol=3)
tfcut = c(0, 5, 10)
cells = unique(X$cell)
for( i in 1:length(cells)){
    Xsub = X[X$cell == cells[i] & X$source == 'him', c('TF_number', 'ppi_density')]
    for(j in 1:3){
        Xsub2 = Xsub[Xsub$TF_number >= tfcut[j], ]
        pval = wilcox.test(Xsub2$ppi_density, mu=0.158, alternative = 'g')$p.value
        Pval[i, j]  = pval
    }
}
Pval
max(Pval)
# stable vs cs hims
feat = c('motif_density',"SON_Sucrose_gene_mean", 'gene_expression_mean', 'csgp')
#feat = c("LaminAC_bin_mean", "LaminAC_gene_mean", "LaminB_bin_mean",  "LaminB_gene_mean","SON_Sucrose_bin_mean","SON_Sucrose_gene_mean")
source('~/Dropbox/r_codes/boxdot.R')
ncol=4; nrow=ceiling(length(feat)/ncol)
pdf(file='main_fig/stable_vs_cs_functional.pdf', width=1.5 * ncol, height=1.8 * nrow)
fi = 'conserve_staVScs'
Xf = X[, c(feat, fi)]
# k562 for essential genes and higher number of super-enhancers
i=4
Xf_c = Xf[X[, 'cell'] == cells[i], ]
print(boxdot(df=Xf_c, fn=fi, g1='stable', g2='cs', cell=cells[i], fnrow=NULL, fncol=ncol, ylabsize=1.1, pvalalter = F, psize=1))
dev.off()
# Main figures stop here

# Sup figures start here
feat = c('motif_density',"SON_Sucrose_gene_mean", 'gene_expression_mean', 'csgp')
source('~/Dropbox/r_codes/boxdot.R')
ncol=4; nrow=length(feat)
pdf(file='sup_fig/stable_vs_cs_functional.pdf', width=1.5 * ncol, height=1.8 * nrow)
fi = 'conserve_staVScs'
Xf = X[, c(feat, fi, 'cell')]
print(boxviolin2(df=Xf, fs='cell', fn=fi, g1='stable', g2='cs', facet.switch='y', ylabsize=1.1, pvalalter = F, psize=1, amaxy=1.05))
dev.off()
# Apercent for constitutive and cell type-specific HIMs
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
# Apercent for constitutive and cell type-specific HIMs stops here 
col_list = c('stable' = "#2b83ba", 'cs' = "#d7191c")
pdf('sup_fig/Apercent_stable_cs.pdf', width=5.85, height=2.3)
ggplot(data=Adynamic, aes(category, proportion, fill=conserve_staVScs)) + 
    geom_col(width=0.5,  alpha=0.5, position=position_dodge()) + 
    ylab('% of HIMs') + xlab('% of genes in A compartments') + 
    scale_y_continuous(expand=c(0, 0.5)) + 
    scale_fill_manual(name="",values = col_list, labels=rename_features('A', nl)) +
    theme(axis.title = element_text(size=10), axis.text=element_text(size=8), 
          axis.text.x=element_text(angle=25, hjust=1), axis.ticks.x=element_blank()) + 
    facet_grid(~cell, labeller = as_labeller(rename_features('A', nl))) + 
    geom_text(data=AdynamicPval, aes(x=x, y=y, label=pval), size=3, inherit.aes = F)  + 
    theme(legend.position='bottom', legend.box.margin=margin(t=-20, r=0, b=-10, l=0),
                                      legend.text=element_text(size=10))
dev.off()
# ChIA-PET
feat = c('ChIA_num_p', 'ChIA_num')
ncol=length(feat); nrow=ceiling(length(feat)/ncol)
pdf(file='sup_fig/him_vs_nonhim_spatial.pdf', width=1.5 * ncol, height=3 * nrow)
fi = 'source'
Xf = X[, c(feat, fi)]
i=1
Xf_c = Xf[X[, 'cell'] == cells[i], ]
print(boxdot(df=Xf_c, fn=fi, g1='him', g2='merged non-him', cell=cells[i], fnrow=NULL, fncol=ncol, ylabsize=1.1, pvalalter = F))
print(boxdot(df=Xf_c, fn=fi, g1='merged him', g2='merged non-him', cell=cells[i], fnrow=NULL, fncol=ncol, ylabsize=1.1, pvalalter = F))
dev.off()

source('~/Dropbox/r_codes/boxdot.R')
feat = c('A_percent', 'repli_mean', 'repli_cv')
nrow=length(feat); ncol=5
pdf(file='sup_fig/sup_him_vs_nonhim_spatial.pdf', width=1.5 * ncol, height=2.2 * nrow)
fi = 'source'
Xf = X[, c(feat, 'cell', fi)]
#Xf = Xf[Xf$cell != 'gm12878', ]
print(boxdot2(df=Xf, fs='cell', fn=fi, g1='him', g2='merged non-him', ylabsize=1.1, pvalalter = F, psize=1, amaxy=1.05))
dev.off()
aggregate(A_percent~cell, data=X, function(z) sum(z>=100)/length(z)*100)

feat = c('SON_Sucrose_gene_mean', 'A_percent', 'repli_mean', 'repli_cv')
nrow=4; ncol=5
pdf(file='sup_fig/sup_mhim_vs_nonhim_spatial.pdf', width=1.5 * ncol, height=2.2 * nrow)
fi = 'source'
Xf = X[, c(feat, 'cell', fi)]
print(boxdot2(df=Xf, fs='cell', fn=fi, g1='merged him', g2='merged non-him', ylabsize=1.1, pvalalter = F, psize=2))
dev.off()
feat = c('ess_all_gene_p', 'ess_k562_specific_gene_p', 'gene_expression_mean', 'gene_expression_cv','SE_num_50kb_per_Mb', 'cosmic_snp_frequency')
nrow=length(feat); ncol=5
source('~/Dropbox/r_codes/boxdot.R')
pdf(file='sup_fig/sup_mhim_vs_nonhim_functional.pdf', width=1.5 * ncol, height=1.8 * nrow)
fi = 'source'
Xf = X[, c(feat, 'cell', fi)]
print(boxdot2(df=Xf, fs='cell', fn=fi, g1='merged him', g2='merged non-him', ylabsize=1.1, pvalalter = F, psize=2))
dev.off()
pdf(file='sup_fig/sup_him_vs_nonhim_functional.pdf', width=1.5 * ncol, height=2.2 * nrow)
fi = 'source'
# remove hims with few genes
ind = X$gene_num <= 6
sum(ind)
X2 = X[!ind, ]
Xf = X2[, c(feat, 'cell', fi)]
print(boxdot2(df=Xf, fs='cell', fn=fi, g1='him', g2='merged non-him', ylabsize=1.1, pvalalter = F, psize=1))
dev.off()
# sup_fig for gene_expression
feat = c('gene_expression_mean', 'gene_expression_cv')
nrow=length(feat); ncol=5
pdf(file='sup_fig/sup_mhim_vs_nonhim_gene_expression.pdf', width=6.5, height=1.4 * nrow)
fi = 'source'
Xf = X[, c(feat, 'cell', fi)]
print(boxdot2(df=Xf, fs='cell', fn=fi, g1='merged him', g2='merged non-him', ylabsize=1.1, pvalalter = F, psize=2,amaxy=1.05 ))
dev.off()
# sup_fig for gene vs gene on gene expression 
feat = c('expression', 'expression')
ncol=5; nrow=length(feat)
pdf(file='sup_fig/gene_gene_gene_expression.pdf', width=6.4, height=1.4 * nrow)
fi = 'HIM'
Gf = G[, c(feat, fi, 'cell')]
Gf = Gf[Gf$expression <= 10, ]
print(boxviolin2(df=Gf, fs='cell', fn=fi, g1='Yes', g2='No', facet.switch='y' ,ylabsize=1.1, pvalalter = F, psize=1, amaxy=1.05))
dev.off()
# sup_fig for gene_expression stops here
# sup_fig for super-enhancers
nl['SE_num_50kb_per_Mb'] = '# super-enhancers per Mb\n(window size=50kb)'
feat = c('SE_num_20kb_per_Mb', 'SE_num_50kb_per_Mb', 'SE_num_100kb_per_Mb','SE_num_500kb_per_Mb','SE_num_1Mb_per_Mb')
# compute fold-change
for(fe in feat){
    print(fe)
    for(cell in unique(X$cell)[1:4]){
        x = X[X$cell == cell & X$source == 'merged him', fe]
        m1 = median(x, na.rm=T)
        y = X[X$cell == cell & X$source == 'merged non-him', fe]
        m2 = median(y, na.rm=T)
        print(c(m1, m2, m1/m2))
    }
}
nrow=length(feat); ncol=4
pdf(file='sup_fig/sup_mhim_vs_nonhim_SE.pdf', width=1.5 * ncol, height=1.8 * nrow)
fi = 'source'
Xf = X[, c(feat, 'cell', fi)]
# remove nhek due to lack of super-enhancer data
Xf = Xf[Xf$cell != 'nhek', ]
print(boxdot2(df=Xf, fs='cell', fn=fi, g1='merged him', g2='merged non-him', ylabsize=1.1, pvalalter = F, psize=2))
dev.off()
# sup_fig for super-enhancers stop here
# sup_fig for SNVs 
feat = c('cosmic_snp_frequency','snp_single_frequency', 'snp_ins_frequency', 'snp_del_frequency')
nrow=length(feat); ncol=5
# reorder cells, HeLa + K562 vs Normal cell types 
pdf(file='sup_fig/sup_mhim_vs_nonhim_snp.pdf', width=1.5 * ncol, height=1.8 * nrow)
fi = 'source'
Xf = X[, c(feat, 'cell', fi)]
Xf$cell = factor(Xf$cell, levels=c('hela', 'k562', 'gm12878', 'huvec', 'nhek'))
print(boxdot2(df=Xf, fs='cell', fn=fi, g1='merged him', g2='merged non-him', ylabsize=1.1, pvalalter = F, psize=2))
dev.off()
# sup_fig for SNVs stops here
feat = c('ess_k562_specific_gene_p')
nrow=length(feat); ncol=1
pdf(file='sup_fig/sup_mhim_vs_nonhim_ess_k562_specific.pdf', width=2 * ncol, height=2.5 * nrow)
fi = 'source'
Xf = X[X$cell=='k562', c(feat, fi)]
print(boxdot_1panel(df=Xf, fn=fi, g1='merged him', g2='merged non-him', ylabsize=1.1, pvalalter = F, psize=2))
dev.off()
feat = c('SE_num_20kb_per_Mb', 'SE_num_100kb_per_Mb','SE_num_500kb_per_Mb','SE_num_1Mb_per_Mb')
nrow=length(feat); ncol=4
pdf(file='sup_fig/sup_mhim_vs_nonhim_SE.pdf', width=1.5 * ncol, height=1.8 * nrow)
fi = 'source'
Xf = X[, c(feat, 'cell', fi)]
# remove nhek due to lack of super-enhancer data
Xf = Xf[Xf$cell != 'nhek', ]
print(boxdot2(df=Xf, fs='cell', fn=fi, g1='merged him', g2='merged non-him', ylabsize=1.1, pvalalter = F, psize=2))
dev.off()
feat = c('snp_single_frequency', 'snp_ins_frequency', 'snp_del_frequency')
nrow=length(feat); ncol=5
pdf(file='sup_fig/sup_mhim_vs_nonhim_snp.pdf', width=1.5 * ncol, height=1.8 * nrow)
fi = 'source'
Xf = X[, c(feat, 'cell', fi)]
print(boxdot2(df=Xf, fs='cell', fn=fi, g1='merged him', g2='merged non-him', ylabsize=1.1, pvalalter = F, psize=2))
dev.off()
# stable vs cs
#feat = c('edge_density', 'tfgene_density', 'motif_density', 'A_percent', 'repli_mean', 'repli_cv')
feat = c('edge_density', 'tfgene_density', 'motif_density', 'repli_mean', 'repli_cv')
nrow=length(feat); ncol=5
#pdf(file='sup_fig/sup_stable_vs_cs_spatial.pdf', width=1.2 * ncol, height=1.8 * nrow)
png(file='sup_fig/sup_stable_vs_cs_spatial.png', width=1.2 * ncol, height=1.8 * nrow, units='in', res=300)
fi = 'conserve_staVScs'
Xf = X[, c(feat, 'cell', fi)]
print(boxdot2(df=Xf, fs='cell', fn=fi, g1='stable', g2='cs', ylabsize=1.1, pvalalter = F, psize=1, amaxy=1.05))
dev.off()
#feat = c('gene_expression_mean', 'csgp', 'ess_all_gene_p', 'SE_num_50kb_per_Mb')
#feat = c('gene_expression_mean', 'csgp', 'SON_Sucrose_gene_mean')
feat = c('gene_expression_mean', 'csgp')
nrow=length(feat); ncol=5
pdf(file='sup_fig/sup_stable_vs_cs_functional.pdf', width=1.2 * ncol, height=1.8 * nrow)
fi = 'conserve_staVScs'
Xf = X[, c(feat, 'cell', fi)]
print(boxdot2(df=Xf, fs='cell', fn=fi, g1='stable', g2='cs', ylabsize=1.1, pvalalter = F, psize=1))
dev.off()
# Sup figures stop here

feat = c(fchara, fexpg, fsnp, fessg, fcsg, fmas, fchia, fse)
ncol=5; nrow=7
pdf(file='fig/him_vs_nonhim.pdf', width=2.5 * ncol, height=4 * nrow)
fi = 'source'
Xf = X[, c(feat, fi)]
for(i in 1:length(cells)){
    Xf_c = Xf[X[, 'cell'] == cells[i], ]
    print(boxdot(df=Xf_c, fn=fi, g1='him', g2='merged non-him', cell=cells[i], fnrow=NULL, fncol=ncol))
}

for(i in 1:length(cells)){
    Xf_c = Xf[X[, 'cell'] == cells[i], ]
    print(boxdot(df=Xf_c, fn=fi, g1='merged him', g2='merged non-him', cell=cells[i], fnrow=NULL, fncol=ncol))
}
dev.off()

fivenum(X[X$source=='him', 'ChIA_num_p'])

pdf(file='fig/stable_vs_cs.pdf', width=2.5 * ncol, height=4 * nrow)
fi = 'conserve_staVScs'
Xf = X[, c(feat, fi)]
for(i in 1:length(cells)){
    Xf_c = Xf[X[, 'cell'] == cells[i], ]
    print(boxdot(df=Xf_c, fn=fi, g1='stable', g2='cs', cell=cells[i], fnrow=NULL, fncol=ncol))
}
dev.off()

# prepare figures for slides and main figures
ncol=3; nrow=2
pdf(file='fig/him_vs_nonhim_spatial.pdf', width=2.5 * ncol, height=4 * nrow)
fi = 'source'
feat = c('edge_density', 'motif_density', 'ChIA_num_p', 'A_percent', 'repli_mean', 'repli_cv')
Xf = X[, c(feat, fi)]
# gm12878
i=1
Xf_c = Xf[X[, 'cell'] == cells[i], ]
print(boxdot(df=Xf_c, fn=fi, g1='him', g2='merged non-him', cell=cells[i], fnrow=NULL, fncol=ncol))
dev.off()
#
ncol=4; nrow=1
pdf(file='fig/him_vs_nonhim_snp.pdf', width=2.5 * ncol, height=4 * nrow)
fi = 'source'
feat = c('snp_single_frequency', 'snp_ins_frequency', 'snp_del_frequency','cosmic_snp_frequency')
Xf = X[, c(feat, fi)]
i=1
Xf_c = Xf[X[, 'cell'] == cells[i], ]
print(boxdot(df=Xf_c, fn=fi, g1='merged him', g2='merged non-him', cell=cells[i], fnrow=NULL, fncol=ncol))
dev.off()
#
#feat = c('gene_expression_mean', 'gene_expression_cv', 'ess_all_gene_p','ess_con_gene_p','csgp', 'SE_num_500kb')
feat = c('gene_expression_mean', 'ess_all_gene_p', 'ess_k562_gene_p', 'csgp', 'gene_expression_cv','ess_con_gene_p', 'ess_k562_specific_gene_p', 'SE_num_500kb')
ncol=4; nrow=ceiling(length(feat) / ncol)
fi = 'source'
pdf(file='fig/him_vs_nonhim_functional.pdf', width=2.5 * ncol, height=4 * nrow)
Xf = X[, c(feat, fi)]
i=4
Xf_c = Xf[X[, 'cell'] == cells[i], ]
print(boxdot(df=Xf_c, fn=fi, g1='merged him', g2='merged non-him', cell=cells[i], fnrow=NULL, fncol=ncol))
dev.off()
# stable vs cs
feat = fchara
ncol=3; nrow=ceiling(length(feat) / ncol)
fi = 'conserve_staVScs'
pdf(file='fig/stable_vs_cs_spatial.pdf', width=2.5 * ncol, height=4 * nrow)
Xf = X[, c(feat, fi)]
i=4
Xf_c = Xf[X[, 'cell'] == cells[i], ]
print(boxdot(df=Xf_c, fn=fi, g1='stable', g2='cs', cell=cells[i], fnrow=NULL, fncol=ncol))
dev.off()
#
feat = fsnp
ncol=4; nrow=ceiling(length(feat) / ncol)
fi = 'conserve_staVScs'
pdf(file='fig/stable_vs_cs_snp.pdf', width=2.5 * ncol, height=4 * nrow)
Xf = X[, c(feat, fi)]
i=4
Xf_c = Xf[X[, 'cell'] == cells[i], ]
print(boxdot(df=Xf_c, fn=fi, g1='stable', g2='cs', cell=cells[i], fnrow=NULL, fncol=ncol))
dev.off()
#
feat = fexpg
ncol=2; nrow=ceiling(length(feat) / ncol)
fi = 'conserve_staVScs'
pdf(file='fig/stable_vs_cs_exp.pdf', width=2.5 * ncol, height=4 * nrow)
Xf = X[, c(feat, fi)]
i=4
Xf_c = Xf[X[, 'cell'] == cells[i], ]
print(boxdot(df=Xf_c, fn=fi, g1='stable', g2='cs', cell=cells[i], fnrow=NULL, fncol=ncol))
dev.off()
#
feat = c(fessg, 'csgp', 'SE_num_500kb')
ncol=4; nrow=ceiling(length(feat) / ncol)
fi = 'conserve_staVScs'
pdf(file='fig/stable_vs_cs_functional.pdf', width=2.5 * ncol, height=4 * nrow)
Xf = X[, c(feat, fi)]
i=4
Xf_c = Xf[X[, 'cell'] == cells[i], ]
print(boxdot(df=Xf_c, fn=fi, g1='stable', g2='cs', cell=cells[i], fnrow=NULL, fncol=ncol))
dev.off()
# figures for slides stop here
# Main figures
feat = c('snp_single_frequency', 'gene_expression_mean', 'ess_all_gene_p', 'ess_k562_specific_gene_p', 'SE_num_500kb', 'csgp')
ncol=3; nrow=ceiling(length(feat) / ncol)
fi = 'source'
pdf(file='fig/main_him_vs_nonhim_functional.pdf', width=2.5 * ncol, height=4 * nrow)
Xf = X[, c(feat, fi)]
i=4
Xf_c = Xf[X[, 'cell'] == cells[i], ]
print(boxdot(df=Xf_c, fn=fi, g1='merged him', g2='merged non-him', cell=cells[i], fnrow=NULL, fncol=ncol, psize=3, ylabsize=2))
dev.off()

ncol=3; nrow=2
pdf(file='fig/main_him_vs_nonhim_spatial.pdf', width=2.5 * ncol, height=4 * nrow)
fi = 'source'
feat = c('edge_density', 'motif_density', 'ChIA_num_p', 'A_percent', 'repli_mean', 'repli_cv')
Xf = X[, c(feat, fi)]
# gm12878
i=1
Xf_c = Xf[X[, 'cell'] == cells[i], ]
print(boxdot(df=Xf_c, fn=fi, g1='him', g2='merged non-him', cell=cells[i], fnrow=NULL, fncol=ncol, psize=1, ylabsize = 2))
dev.off()

ncol=3; nrow=2
pdf(file='fig/main_mhim_vs_mnonhim_spatial.pdf', width=2.5 * ncol, height=4 * nrow)
fi = 'source'
feat = c('edge_density', 'motif_density', 'ChIA_num_p', 'A_percent', 'repli_mean', 'repli_cv')
Xf = X[, c(feat, fi)]
# gm12878
i=1
Xf_c = Xf[X[, 'cell'] == cells[i], ]
print(boxdot(df=Xf_c, fn=fi, g1='merged him', g2='merged non-him', cell=cells[i], fnrow=NULL, fncol=ncol, psize=3, ylabsize = 2))
dev.off()

feat = c('snp_single_frequency', 'gene_expression_mean', 'ess_all_gene_p', 'ess_k562_specific_gene_p', 'SE_num_500kb', 'csgp')
ncol=3; nrow=ceiling(length(feat) / ncol)
fi = 'conserve_staVScs'
pdf(file='fig/main_stable_vs_cs_functional.pdf', width=2.5 * ncol, height=4 * nrow)
Xf = X[, c(feat, fi)]
i=4
Xf_c = Xf[X[, 'cell'] == cells[i], ]
print(boxdot(df=Xf_c, fn=fi, g1='stable', g2='cs', cell=cells[i], fnrow=NULL, fncol=ncol, psize=1, ylabsize = 2))
dev.off()

feat = fchara
ncol=3; nrow=ceiling(length(feat) / ncol)
fi = 'conserve_staVScs'
pdf(file='fig/stable_vs_cs_spatial.pdf', width=2.5 * ncol, height=4 * nrow)
Xf = X[, c(feat, fi)]
i=4
Xf_c = Xf[X[, 'cell'] == cells[i], ]
print(boxdot(df=Xf_c, fn=fi, g1='stable', g2='cs', cell=cells[i], fnrow=NULL, fncol=ncol, psize=1, ylabsize = 2))
dev.off()

# Test super-enhancers
#feat= fsef[!grepl('above', fsef)]
#feat = fsef[grepl('per_gene', fsef)]
#feat= fsef[grepl('above_10', fsef)]
feat= fsef[grepl('density', fsef) & !grepl('per', fsef)]
ncol=10; nrow=ceiling(length(feat) / ncol)
pdf(file='fig/him_vs_nonhim_SE.pdf', width=2.5 * ncol, height=4 * nrow)
fi = 'source'
Xf = X[, c(feat, fi)]
for(i in 1:length(cells)){
    Xf_c = Xf[X[, 'cell'] == cells[i], ]
    print(boxdot(df=Xf_c, fn=fi, g1='him', g2='merged non-him', cell=cells[i], fnrow=NULL, fncol=ncol))
}
for(i in 1:length(cells)){
    Xf_c = Xf[X[, 'cell'] == cells[i], ]
    print(boxdot(df=Xf_c, fn=fi, g1='merged him', g2='merged non-him', cell=cells[i], fnrow=NULL, fncol=ncol))
}
dev.off()
pdf(file='fig/stable_vs_cs_SE.pdf', width=2.5 * ncol, height=4 * nrow)
fi = 'conserve_staVScs'
Xf = X[, c(feat, fi)]
for(i in 1:length(cells)){
    Xf_c = Xf[X[, 'cell'] == cells[i], ]
    print(boxdot(df=Xf_c, fn=fi, g1='stable', g2='cs', cell=cells[i], fnrow=NULL, fncol=ncol))
}
dev.off()
# histogram
Xf_c = Xf[X[, 'cell'] == cells[4], ]
ft = 'SE_above_10_500kb_per_Mb'
p = ggplot(subset(Xf_c, SE_above_10_500kb_per_Mb>0), aes_string(x=ft, color=fi)) + 
    geom_histogram(fill='white', alpha=0.5, position='dodge')
p
# barplot
Xf_c = Xf[X[, 'cell'] == cells[4], ]
ft = 'SE_above_10_1Mb'
tmp = as.data.frame(table(Xf_c[, 'conserve_staVScs'], Xf_c[, ft]))
colnames(tmp) = c('conserve_staVScs', ft, 'Freq')
p = ggplot(tmp, aes_string(x=ft, y='Freq', fill='conserve_staVScs')) + 
    geom_bar(alpha=0.5, position=position_dodge(), stat='identity')  
p

# test super-enhancers while adjusting the genne num and gene distance
cn_ind = 'conserve_staVScs'
X[, cn_ind] = factor(X[, cn_ind], levels=c('stable', 'cs'))
cn_control = c('gene_num', 'gene_distance')
cn_ind = c(cn_control, 'conserve_staVScs')
cn_dep = 'SE_num_500kb'
np_adjust_test_staVScs(df=subset(X, cell!='nhek'), cn_dep='SE_num_100kb', cn_ind=cn_ind, cn_control = cn_control)
np_adjust_test_staVScs(df=subset(X, cell!='nhek'), cn_dep='SE_num_500kb', cn_ind=cn_ind, cn_control = cn_control)
np_adjust_test_staVScs(df=subset(X, cell!='nhek'), cn_dep='SE_num_1Mb', cn_ind=cn_ind, cn_control = cn_control)

##### GO on gene assignment
x = c(11627, 12036, 11927, 12391, 12161)
y = c(69.1, 77.2, 75.3, 76.5, 62.1) / 100
z = c(344, 441, 365, 335, 448)
w = round(z / (x * y) * 100, 2)
u = round(3025 / (x * y) * 100, 2)
range(u)
