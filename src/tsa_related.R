library(grid)
library(gtable)
Xk = X[X$cell == 'k562' & X$source == 'him', ]
Xk2 = X[X$cell == 'k562', ]
col_list = c('him' = "#2b83ba", 'merged non-him' = "#d7191c")
featn = c('him', 'merged non-him')
coor_541 = X['k562_541', c('SON_Sucrose_gene_mean', 'LaminB_gene_mean')]
dfarrow = data.frame(coor_541[1], coor_541[1] + 0.5, coor_541[2], coor_541[2]+0.5)
colnames(dfarrow) = c('x1', 'x2', 'y1', 'y2')
# scatter plot + CDF
tsa.scatter = ggplot(subset(Xk2, source=='him'), 
  aes(x=SON_Sucrose_gene_mean, y=LaminB_gene_mean, color=source)) +
  xlab('Mean SON TSA-seq') + 
  ylab('Mean Lamin B TSA-seq') + 
  geom_point(size=0.5, alpha=0.5) + 
    geom_vline(xintercept = 0.284, linetype='dotted', size=0.5, colour='red') + 
  scale_color_manual(name="",values = col_list, labels=rename_features(featn, nl)) +
  theme(legend.position = 'none') + 
    theme(axis.title = element_text(size=9), axis.text=element_text(size=7)) + 
  geom_segment(data=dfarrow, aes(xend=x1, x=x2, yend=y1, y=y2), size=1, alpha=0.5,
                 arrow=arrow(length=unit(0.05, 'npc')), inherit.aes=F )  
tsa.scatter
# CDF
tsa.son.density = ggplot(subset(Xk2, source=='him'),
 aes(SON_Sucrose_gene_mean, fill=source, colour=source)) + 
 stat_ecdf(pad=FALSE) + 
 geom_vline(xintercept = 0.284, linetype='dotted', size=0.5, colour='red') + 
 geom_hline(yintercept = 0.4, linetype='dotted', size=0.5, colour='red') + 
 scale_fill_manual(name="",values = col_list, labels=rename_features(featn, nl)) +
 scale_color_manual(name="",values = col_list, labels=rename_features(featn, nl)) +
 theme(legend.position = 'none') + 
 ylab('CDF')  + 
 theme(axis.ticks.x = element_blank(), axis.text.x=element_blank(), axis.title.x = element_blank(),
       axis.text.y=element_text(size=6), axis.title.y=element_text(size=6)) 
tsa.lam.density = ggplot(subset(Xk2, source=='him'),
 aes(LaminB_gene_mean, fill=source, colour=source)) + 
 stat_ecdf(pad=FALSE) + 
 scale_fill_manual(name="",values = col_list, labels=rename_features(featn, nl)) +
 scale_color_manual(name="",values = col_list, labels=rename_features(featn, nl)) +
 theme(legend.position = 'none') + 
 ylab('CDF')  + 
 coord_flip() + 
 theme(axis.ticks.y = element_blank(), axis.text.y=element_blank(), 
       axis.title.y = element_blank(), axis.text.x=element_text(size=6, angle=90, hjust=1),
       axis.title.x=element_text(size=6)) 
tsa.lam.density
# combine into one figure
g <- ggplotGrob(tsa.scatter)
panel_id <- g$layout[g$layout$name == "panel",c("t","l")]
g <- gtable_add_cols(g, unit(0.2,"npc"))    
g <- gtable_add_grob(g, ggplotGrob(tsa.lam.density), t = panel_id$t-4, b=panel_id$t+4, l=ncol(g))
g <- gtable_add_rows(g, unit(0.2,"npc"), 0)
g <- gtable_add_grob(g, ggplotGrob(tsa.son.density), t = 1, l=panel_id$l, r=panel_id$l + 2)
pdf('main_fig/tsa_scatter.pdf', width=3, height=3)
grid.newpage()
grid.draw(g)
dev.off()

feat = c('SON_Sucrose_gene_mean', 'LaminB_gene_mean')
ncol=length(feat); nrow=ceiling(length(feat)/ncol)
pdf(file='sup_fig/tsa_mhim_mnonhim.pdf', width=1.5 * ncol, height=1.8 * nrow)
fi = 'source'
Xf = X[, c(feat, fi)]
i=4; Xf_c = Xf[X[, 'cell'] == cells[i], ]
print(boxdot(df=Xf_c, fn=fi, g1='merged him', g2='merged non-him', cell=cells[i], fnrow=NULL, fncol=ncol, ylabsize=1.1, pvalalter = F))
dev.off()
feat = c('SON_Sucrose',"LaminB")
ncol=length(feat); nrow=ceiling(length(feat)/ncol)
pdf(file='sup_fig/tsa_gene_gene.pdf', width=1.5*ncol, height=1.8 * nrow)
fi = 'HIM'
Gf = G[, c(feat, fi)]
Gf_c = Gf[G$cell == 'k562', ]
print(boxviolin(df=Gf_c, fn=fi, g1='Yes', g2='No', cell=cells[i], fnrow=NULL, fncol=ncol, ylabsize=1.1, pvalalter = F))
dev.off()
median(Gf_c[Gf_c$HIM == 'Yes', 'SON_Sucrose'])
median(Gf_c[Gf_c$HIM == 'No', 'SON_Sucrose'])

# TSA on constitutive and cell type-specific HIMs
# there is no visible difference between the two types of HIMs
# commet out on 26 Dec 2018
#library(grid)
#library(gtable)
#col_list = c('stable' = "#2b83ba", 'cs' = "#d7191c")
#featn = c('stable', 'cs')
## scatter plot + CDF
#tsa.scatter = ggplot(subset(Xk2, source=='him'), 
#  aes(x=SON_Sucrose_gene_mean, y=LaminB_gene_mean, color=conserve_staVScs)) +
#  xlab('Mean SON TSA-seq') + 
#  ylab('Mean Lamin B TSA-seq') + 
#  geom_point(size=0.5, alpha=0.5) + 
#    geom_vline(xintercept = 0.284, linetype='dotted', size=0.5, colour='red') + 
#  scale_color_manual(name="",values = col_list, labels=rename_features(featn, nl)) +
#  theme(legend.position = 'none') + 
#    theme(axis.title = element_text(size=9), axis.text=element_text(size=7)) 
#tsa.scatter
#tsa.son.density = ggplot(subset(Xk2, source=='him'),
# aes(SON_Sucrose_gene_mean, fill=conserve_staVScs, colour=conserve_staVScs)) + 
# geom_density(alpha=0.5) + 
# geom_vline(xintercept = 0.284, linetype='dotted', size=0.5, colour='red') + 
# scale_fill_manual(name="",values = col_list, labels=rename_features(featn, nl)) +
# scale_color_manual(name="",values = col_list, labels=rename_features(featn, nl)) +
# theme(legend.position = 'none') + 
# ylab('Density')  + 
# theme(axis.ticks.x = element_blank(), axis.text.x=element_blank(), axis.title.x = element_blank(),
#       axis.text.y=element_text(size=6), axis.title.y=element_text(size=6)) 
#tsa.lam.density = ggplot(subset(Xk2, source=='him'),
# aes(LaminB_gene_mean, fill=conserve_staVScs, colour=conserve_staVScs)) + 
# geom_density(alpha=0.5) + 
# scale_fill_manual(name="",values = col_list, labels=rename_features(featn, nl)) +
# scale_color_manual(name="",values = col_list, labels=rename_features(featn, nl)) +
# theme(legend.position = 'none') + 
# ylab('Density')  + 
# coord_flip() + 
# theme(axis.ticks.y = element_blank(), axis.text.y=element_blank(), 
#       axis.title.y = element_blank(), axis.text.x=element_text(size=6, angle=90, hjust=1),
#       axis.title.x=element_text(size=6)) 
#tsa.lam.density
## combine into one figure
#g <- ggplotGrob(tsa.scatter)
#panel_id <- g$layout[g$layout$name == "panel",c("t","l")]
#g <- gtable_add_cols(g, unit(0.2,"npc"))    
#g <- gtable_add_grob(g, ggplotGrob(tsa.lam.density), t = panel_id$t-4, b=panel_id$t+4, l=ncol(g))
#g <- gtable_add_rows(g, unit(0.2,"npc"), 0)
#g <- gtable_add_grob(g, ggplotGrob(tsa.son.density), t = 1, l=panel_id$l, r=panel_id$l + 2)
#pdf('sup_fig/tsa_scatter_constitutive_cs.pdf', width=3, height=3)
#grid.newpage()
#grid.draw(g)
#dev.off()

#tsa.plot = ggplot(Xk, aes(x=motif_density, y=SON_Sucrose_gene_mean, color=conserve_staVScs)) +
#    geom_point() + geom_smooth(method=lm)
#print(tsa.plot)
#tsa.plot = ggplot(Xk, aes(x=edge_density, y=SON_Sucrose_gene_mean, color=conserve_staVScs)) +
#    geom_point()+  geom_smooth(method=lm)
#print(tsa.plot)
#tsa.plot = ggplot(Xk, aes(x=tfgene_density, y=SON_Sucrose_gene_mean, color=conserve_staVScs)) +
#    geom_point()+  geom_smooth(method=lm)
#print(tsa.plot)
#tsa.plot = ggplot(Xk, aes(x=gene_expression_mean, y=SON_Sucrose_gene_mean, color=conserve_staVScs)) +
#    geom_point()+  geom_smooth(method=lm)
#print(tsa.plot)
#tsa.plot = ggplot(Xk, aes(x=SE_num_20kb_per_Mb, y=SON_Sucrose_gene_mean, color=conserve_staVScs)) +
#    geom_point()+  geom_smooth()
#print(tsa.plot)
#tsa.plot = ggplot(Xk, aes(x=cosmic_snp_frequency, y=SON_Sucrose_gene_mean)) +
#    geom_point()+  geom_smooth()
#print(tsa.plot)
#tsa.plot = ggplot(Xk, aes(x=TF_number, y=SON_Sucrose_bin_mean)) +
#    geom_point()+  geom_smooth()
#print(tsa.plot)