###### Barplot on A compartments
# 1. split by cutoffs, such as 50, 100
cutoffs = c(0, 50, 80,   99.99999999999,101)
labels = c('[0, 50)', '[50, 80)', '[80, 100)', '100' )
# A percent for HIMs
Apercent = extract_freq(df=subset(X, source=='him'), cutoffs=cutoffs, labels=labels)
featn = unique(X$cell)
Apercent[Apercent$category == '100', ]
Apercent[Apercent$category == '[0, 50)', ]
# main figure
pdf('main_fig/Apercent_k562.pdf', width=1.4, height = 2.5/2)
ggplot(data=subset(Apercent, cell=='k562'), aes(category, proportion)) + geom_col(width=0.5, fill=col_list[1], alpha=0.5) + 
    ylab('% of HIMs') + xlab('% of genes in A compartments') + 
    scale_y_continuous(expand=c(0, 0.5)) + 
    theme(axis.title = element_text(size=8), axis.text=element_text(size=6), 
          axis.text.x=element_text(angle=25, hjust=1), axis.ticks.x=element_blank())
dev.off()
# sup_fig
pdf('sup_fig/Apercent.pdf', width=6.5, height=2)
ggplot(data=Apercent, aes(category, proportion)) + geom_col(width=0.5, fill=col_list[1], alpha=0.5) + 
    ylab('% of HIMs') + xlab('% of genes in A compartments') + 
    scale_y_continuous(expand=c(0, 0.5)) + 
    theme(axis.title = element_text(size=10), axis.text=element_text(size=8), 
          axis.text.x=element_text(angle=25, hjust=1), axis.ticks.x=element_blank()) + 
    facet_grid(~cell, labeller = as_labeller(rename_features(featn, nl)))
dev.off()
