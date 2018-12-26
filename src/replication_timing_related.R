col_list = c('Yes' = "#2b83ba", 'No' = "#d7191c")
pdf('main_fig/repli_seq_gene_vs_gene.pdf', width=1.5, height = 2.5)
p = ggplot(data=subset(G, cell=='k562'), aes(x=HIM, y=repli, color=HIM)) + 
    geom_violin() + geom_boxplot(width=0.1, outlier.shape = NA) + 
    scale_colour_manual(name="",values = col_list) +
    xlab('Gene assignment to HIMs') + ylab('Replication timing') + 
    theme(legend.position = 'none') + 
    scale_x_discrete(breaks=c(0, 1), labels=c('No', 'Yes')) + 
    geom_text(data=subset(Gpval_table, cell=='k562'), aes(x=1.5, y=maxy, label=pval), size=3, inherit.aes = F)  + 
    theme(axis.title = element_text(size=8), axis.text=element_text(size=6), 
          axis.ticks.x=element_blank())
print(p)
dev.off()    
#pdf('sup_fig/repli_seq_gene_vs_gene.pdf', width=6.5, height = 2)
#p = ggplot(data=G, aes(x=HIM, y=repli, color=HIM)) + 
#    geom_violin() + geom_boxplot(width=0.05, outlier.shape = NA) + 
#    scale_colour_manual(name="",values = col_list) +
#    facet_grid(~cell, labeller = as_labeller(rename_features(c('A'), nl))) + 
#    geom_text(data=Gpval_table, aes(x=1.5, y=maxy, label=pval), size=3, inherit.aes = F)  + 
#    xlab('Gene assignment to HIMs') + ylab('Replication timing') + 
#    theme(legend.position = 'none') + 
#    scale_x_discrete(breaks=c(0, 1), labels=c('No', 'Yes')) + 
#    theme(axis.title = element_text(size=10), axis.text=element_text(size=8), 
#          axis.ticks.x=element_blank())
#print(p)
#dev.off()    

feat = c('repli', 'repli')
ncol=5; nrow=length(feat)
pdf(file='sup_fig/repli_seq_gene_vs_gene.pdf', width=6.65, height=1.4 * nrow)
fi = 'HIM'
Gf = G[, c(feat, fi, 'cell')]
print(boxviolin2(df=Gf, fs='cell', fn=fi, g1='Yes', g2='No', ylabsize=1.1, pvalalter = F, psize=1, amaxy=1.05))
dev.off()
feat = c('repli_mean', 'repli_cv')
nrow=length(feat); ncol=5
pdf(file='sup_fig/repli_mhim_mnonhim.pdf', width=6.65, height=1.4 * nrow, useDingbats = FALSE)
fi = 'source'
Xf = X[, c(feat, 'cell', fi)]
print(boxdot2(df=Xf, fs='cell', fn=fi, g1='merged him', g2='merged non-him', ylabsize=1.1, pvalalter = F, psize=2))
dev.off()