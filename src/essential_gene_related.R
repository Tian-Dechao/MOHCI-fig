# sup_fig for essential genes 
diffess = X[X$source =='merged him', 'ess_all_gene_p'] - X[X$source =='merged non-him', 'ess_all_gene_p']
fivenum(diffess)

#Xsub = subset(X, X$source %in% c('merged him', 'merged non-him'))
#ggplot(data=Xsub, aes(x=chrom, y=ess_all_gene_p, colour=source, fill=source)) + 
#    geom_col(position=position_dodge()) + 
#    facet_grid(cell~.)
Xsubk = X[X$cell=='k562' & X$source %in%c('merged him', 'merged non-him'), 
          c('chrom','source', 'ess_k562_gene_p')]
Xsubk$chrom = factor(Xsubk$chrom, levels=paste('chr', c(1:22, 'X'), sep=''), ordered = T)
pdf('sup_fig/sup_mhim_vs_nonhim_ese_k562.pdf', width=6.5, height=2)
col_list = c("#2b83ba", "#d7191c")
names(col_list) <- c('merged him', 'merged non-him')
p = ggplot(data=Xsubk, aes(x=chrom, y=ess_k562_gene_p, fill=source, color=source))  + 
    geom_col(position=position_dodge(), width=0.5) + 
    scale_fill_manual(name = '', values = col_list, labels=rename_features('A', nl) ) +
    scale_color_manual(name = '', values = col_list, labels=rename_features('A', nl) ) +
    #xlab('Chromosome') + ylab('% of genes that\nare essential in K562') + 
    xlab('') + ylab('% of genes that\nare essential in K562') + 
    theme(axis.title.y = element_text(size=10), axis.text.y=element_text(size=8)) + 
    theme(axis.text.x=element_text(angle=30, size=8)) + 
    theme(plot.margin = unit(c(0,0,0.1,-0.05), "cm")) + 
    theme(legend.position="bottom", legend.box.margin=margin(t=-20, r=0, b=-10, l=0),
                                      legend.text=element_text(size=10))
print(p)
dev.off()
pdf('main_fig/gene_gene_essential_k562.pdf', width=1.5, height=2.67)
col_list = c("#2b83ba", "#d7191c")
names(col_list) <- c('Yes', 'No')
p = ggplot(data=propEk, aes(x=HIM, y=essp, fill=HIM, color=HIM)) + 
    geom_col(width=0.5) + 
    scale_fill_manual(name = '', values = col_list, labels=rename_features('A', nl) ) +
    scale_color_manual(name = '', values = col_list, labels=rename_features('A', nl) ) +
    xlab('Assignment to HIMs') + ylab('% of genes that\nare essential in K562') + 
    geom_text(data=Ek.p.df, aes(x=x, y=y, label=pval), size=3, inherit.aes = F) + 
    theme(axis.title.y = element_text(size=10), axis.text.y=element_text(size=8)) + 
    theme(axis.text.x=element_text(size=8)) + 
    theme(plot.margin = unit(c(0,0,0.1,-0.05), "cm")) + 
    theme(legend.position="none")
print(p)
dev.off()
pdf('sup_fig/sup_gene_gene_essential_all.pdf', width=6.5, height=2)
col_list = c("#2b83ba", "#d7191c")
names(col_list) <- c('Yes', 'No')
p = ggplot(data=propE.long, aes(x=HIM, y=prop, fill=HIM, color=HIM)) + 
    geom_col(width=0.5) + 
    scale_fill_manual(name = '', values = col_list, labels=rename_features('A', nl) ) +
    scale_color_manual(name = '', values = col_list, labels=rename_features('A', nl) ) +
    facet_grid(.~cell, labeller = as_labeller(rename_features('A', nl))) + 
    xlab('Assignment to HIMs') + ylab('% of genes that\nare essential') + 
    geom_text(data=propE.p.df, aes(x=x, y=y, label=pval), size=3, inherit.aes = F) + 
    theme(axis.title.y = element_text(size=10), axis.text.y=element_text(size=8)) + 
    theme(axis.text.x=element_text(size=8)) + 
    theme(plot.margin = unit(c(0,0,0.1,-0.05), "cm")) + 
    theme(legend.position="none")
print(p)
dev.off()
# Stop here
#pdf('sup_fig/sup_mhim_vs_nonhim_ese_all.pdf', width=6.5, height=6)
#Xsub = X[X$source %in%c('merged him', 'merged non-him'), 
#          c('chrom','source', 'cell', 'ess_all_gene_p')]
#Xsub$chrom = factor(Xsub$chrom, levels=paste('chr', c(1:22, 'X'), sep=''), ordered = T)
#col_list = c("#2b83ba", "#d7191c")
#names(col_list) <- c('merged him', 'merged non-him')
#ggplot(data=Xsub, aes(x=chrom, y=ess_all_gene_p, fill=source, color=source))  + 
#    geom_col(position=position_dodge(), width=0.5) + 
#    scale_fill_manual(name = '', values = col_list, labels=rename_features('A', nl) ) +
#    scale_color_manual(name = '', values = col_list, labels=rename_features('A', nl) ) +
#    facet_grid(cell~., labeller=as_labeller(rename_features('A', nl)) , shrink=T) + 
#    xlab('Chromosome') + ylab('% of genes that are essential') + 
#    theme(axis.title.y = element_text(size=10), axis.text.y=element_text(size=8)) + 
#    theme(axis.text.x=element_text(angle=30, size=8)) + 
#    theme(plot.margin = unit(c(0,0,0.1,-0.05), "cm")) + 
#    theme(legend.position="bottom", legend.box.margin=margin(t=-20, r=0, b=-10, l=0),
#                                      legend.text=element_text(size=10))
#dev.off()
## sup_fig for essential genes stops here
#