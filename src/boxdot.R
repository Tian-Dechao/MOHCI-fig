library(ggplot2)
library(scales)
library(reshape2)
library(plyr)
#library(venn)

nl = c('edge_density' = "Hi-C edge density",
       'tfgene_density' = 'GRN edge density',
       'motif_density' ='Motif density',
       'A_percent' = "% genes in A compartments",
       'repli_mean' = 'Mean replication timing',
       'repli_cv' = 'CV replication timing',
       'ppi_density' = 'PPI density',
       'gene_expression_mean' = 'Mean gene expression',
       'gene_expression_cv' = 'CV gene expression',
       'ess_all_gene_p' = '% genes that are essential',
       'ess_con_gene_p' = '% genes that are\nconsensus essential genes',
       'ess_k562_gene_p'= "% gnees that\nare essential in K562",
       'ess_k562_specific_gene_p' = '% genes that are K562\nspecific essential genes',
       'TF_number' = '# TFs',
       'gene_num' = '# genes',
       'master_in' = '# master TFs',
       'master_in_p' = '% TFs that are master TFs',
       'cosmic_snp_frequency' = 'cosmic SNV density',
       'snp_single_frequency' = 'Substitution SNV density\n(dbSNP)',
       'snp_ins_frequency' = 'Insertion SNV density\n(dbSNP)',
       'snp_del_frequency' = 'Deletion SNV density\n(dbSNP)',
       'ChIA_num' = '# ChIA-PET interactions',
       'ChIA_num_p' = '# ChIA-PET\ninteractions per Mb', 
       'csgp' = '% genes that are\ncell type-specific',
       'cstp' = '% tfs that are\ncell type-specific',
       'SE_num_50kb_per_Mb' = '# super-enhancers per Mb',
       'SE_num_100kb' = '# super-enhancers\ninteracting genes in a HIM (100kb)',
       'SE_density_100kb' = 'Density of super-enhancer--gene\ninteractions (100kb)',
       'SE_num_500kb' = '# super-enhancers\ninteracting genes in a HIM (500kb)',
       'SE_density_500kb' = 'Density of super-enhancer--gene\ninteractions (500kb)',
       'SE_num_1Mb' = '# super-enhancers\ninteracting genes in a HIM (1Mb)',
       'SE_density_1Mb' = 'Density of super-enhancer--gene\ninteractions (1Mb)',
       'SE_num_1Mb_per_Mb' = '# super-enhancers per Mb\n(window size=1Mb)',
       'SE_num_500kb_per_Mb' = '# super-enhancers per Mb\n(window size=500kb)',
       'SE_num_100kb_per_Mb' = '# super-enhancers per Mb\n(window size=100kb)',
       'SE_num_20kb_per_Mb' = '# super-enhancers per Mb\n(window size=20kb)',
       'gm12878' = 'GM12878',
       'hela' = 'HeLa',
       'huvec' = 'HUVEC',
       'k562' = 'K562',
       'nhek' = 'NHEK', 
       'merged him' = 'Merged-HIM',
       'merged non-him' = 'Non-HIM',
       'him' = 'HIM',
       #'stable' = 'Constitutive HIM',
       'stable' = 'More conserved HIM',
       'cs' = 'Cell type-specific HIM',
       '4-node' = 'HIM by the 4-node motif M',
       'triangle' = 'HIM by the triangle motif',
       'bifan' = 'HIM by the bifan motif',
       'tri_density' = 'Triangle density',
       'bifan_density' = 'Bifan density',
       'gene_distance' = '1D length (Mb) of the\nregion covered by\nthe genes in a HIM',
       'SON_Sucrose_gene_mean' = 'Mean SON TSA-seq',
       'LaminB_gene_mean' = 'Mean Lamin B TSA-seq',
       'SON_Sucrose' = 'SON TSA-seq',
       'LaminB' = 'Lamin B TSA-seq',
       'repli' = 'Replication timing'
       )

rename_features = function(x, nl){
    y = x
    names(y) = x
    y[names(nl)] = nl
    return(y)
}

rename_features_v2 = function(x, nl){
    y = x
    names(y) = x
    ysub = y[y%in%names(nl)]
    y[y%in%names(nl)] = nl[y[y%in%names(nl)]]
    return(y)
}

compare2vect = function(x, y){
    if(sum(!is.na(x) & is.finite(x))>5 & sum(!is.na(y) & is.finite(y))>5){
        pg = wilcox.test(x, y, alternative = 'g')$p.value
        pl = wilcox.test(x, y, alternative = 'l')$p.value
        if(pg < pl){
            pval = pg
            alter = '>'
        } else {
            pval = pl
            alter = '<'
        }
    }
    else {
        pval = NA
        alter = NA
    }
    return(data.frame(pval=pval, alter=alter, stringsAsFactors = F))
}

boxdot = function(df, fn='source', g1='him', g2='merged non-him', cell='F!', fnrow=1, fncol=NULL, psize=1, ylabsize=1.2, pvalalter=T){
    col_list = c(g1 = "#2b83ba", g2 = "#d7191c")
    names(col_list) <- c(g1, g2)
    # filter out x
    ind = df[, fn] %in% c(g1, g2)
    df = df[ind, ]
    featn =  colnames(df)[colnames(df)!=fn]  
    # calculate pvals 
    pval = c()
    alter = c()
    for(fi in 1:length(featn)){
        x = df[df[, fn] == g1, featn[fi]]
        y = df[df[, fn] == g2, featn[fi]]
        ptmp = compare2vect(x, y)
        pval[fi] = ptmp$pval
        alter[fi] = ptmp$alter
    }
    miny = apply(df[, featn], 2, min, na.rm=T)
    maxy = apply(df[, featn], 2, max, na.rm=T)
    maxy = maxy + 0.05 * (maxy - miny)
    pval_table = data.frame(feature=featn, pval=pval, alter=alter, maxy=maxy )
    pval_table$pval = format(pval_table$pval, scientific=T, digits=3)
    if(pvalalter){
        pval_table$pval = paste('P=', pval_table$pval, '\n', pval_table$alter, sep='')
    } else {
        pval_table$pval = paste('P=', pval_table$pval, sep='')
    }
    indna = grepl('NA', pval_table[, 'pval'])
    pval_table[indna, 'pval'] = ''
    
    df <- melt(df, id.vars = fn, measure.vars = featn, variable.name = "feature", value.name="value")
    df[, fn] = factor(df[, fn], levels=c(g1, g2), ordered=T)
    # calculate pvals
    #pval_table = data.frame(feature=as.factor(featn),  pval_single=letters[1:4], x=rep(1.5, 4), y=rep(0.1, 4), value=fessg)
    
    plot_compare <- ggplot(df, aes_string(x = fn, y = "value", fill = fn)) +
      geom_jitter(alpha = 0.7, color="black", stroke = 0.0, shape=21, size=psize,  width = 0.2, height = 0.0) +
      geom_boxplot(data=df, aes_string(x = fn, y = "value", color = fn), outlier.shape = NA, width=0.6, color = "black", fill="white") +
      scale_fill_manual(name="",values = col_list, labels=rename_features(featn, nl)) +
      xlab(label='') +
      ylab("") +
      facet_wrap(~feature, labeller=as_labeller(rename_features(featn, nl)), nrow = fnrow, ncol=fncol, scales = "free_y", strip.position = 'left', shrink=T) +
      theme_bw() +
      theme(plot.title = element_text(lineheight=.8, face="bold", size=rel(ylabsize))) +
      theme(panel.spacing = unit(0.5, 'lines')) + 
      theme(panel.border = element_blank()) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.title.y = element_text(size=rel(ylabsize * 0.9),margin=margin(0,1,0,0))) +
      theme(axis.title.x = element_text(size=rel(1.2),margin=margin(1,0,0,0))) +
      theme(axis.text.y = element_text(color='black', size=7)) +
      theme(axis.line = element_line(color = "black")) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      theme(strip.text.x = element_text(size=rel(1.2), face="bold")) +
      theme(strip.text.y = element_text(size=rel(ylabsize))) +
      theme(strip.background = element_blank())  +
      theme(strip.placement = 'outside')  + 
      theme(plot.margin = unit(c(0,0,0.1,-0.5), "cm"))
      plot_compare = plot_compare + theme(legend.position="bottom", legend.box.margin=margin(t=-22, r=0, b=-10, l=0))
      plot_compare <- plot_compare +
      #geom_text(data=pval_table, aes(x=1.5, y=maxy, label=pval), size=3, position=position_dodge(0.9), inherit.aes = F) 
      geom_text(data=pval_table, aes(x=1.5, y=maxy, label=pval), size=3, inherit.aes = F) 
          
      
    return(plot_compare)
}

boxdot_1panel = function(df, fn='source', g1='him', g2='merged non-him', cell='F!', fnrow=1, fncol=NULL, psize=1, ylabsize=1.2, pvalalter=T){
    print(fnrow)
    col_list = c(g1 = "#2b83ba", g2 = "#d7191c")
    names(col_list) <- c(g1, g2)
    # filter out x
    ind = df[, fn] %in% c(g1, g2)
    df = df[ind, ]
    featn =  colnames(df)[colnames(df)!=fn]  
    # calculate pvals 
    pval = c()
    alter = c()
    for(fi in 1:length(featn)){
        x = df[df[, fn] == g1, featn[fi]]
        y = df[df[, fn] == g2, featn[fi]]
        ptmp = compare2vect(x, y)
        pval[fi] = ptmp$pval
        alter[fi] = ptmp$alter
    }
    miny = min(df[, featn], na.rm=T)
    maxy = max(df[, featn], na.rm=T)
    maxy = maxy + 0.05 * (maxy - miny)
    pval_table = data.frame(feature=featn, pval=pval, alter=alter, maxy=maxy )
    pval_table$pval = format(pval_table$pval, scientific=T, digits=3)
    if(pvalalter){
        pval_table$pval = paste('P=', pval_table$pval, '\n', pval_table$alter, sep='')
    } else {
        pval_table$pval = paste('P=', pval_table$pval, sep='')
    }
    indna = grepl('NA', pval_table[, 'pval'])
    pval_table[indna, 'pval'] = ''
    
    df <- melt(df, id.vars = fn, measure.vars = featn, variable.name = "feature", value.name="value")
    df[, fn] = factor(df[, fn], levels=c(g1, g2), ordered=T)
    
    plot_compare <- ggplot(df, aes_string(x = fn, y = "value", fill = fn)) +
      geom_jitter(alpha = 0.7, color="black", stroke = 0.0, shape=21, size=psize,  width = 0.2, height = 0.0) +
      geom_boxplot(data=df, aes_string(x = fn, y = "value", color = fn), outlier.shape = NA, width=0.6, color = "black", fill="white") +
      xlab(label='') +
      ylab(label=rename_features(featn, nl)) +
      theme(axis.text.y = element_text(color='black', size=7)) +
      theme(axis.text.x = element_text(color='black', size=10, angle=45, hjust=1, vjust=1)) +
      scale_x_discrete(labels=rename_features_v2(c(g1, g2), nl)) +
      theme(plot.margin = unit(c(0,0,-0.3,0.1), "cm")) +
        guides(fill=F)
      plot_compare <- plot_compare +
      geom_text(data=pval_table, aes(x=1.5, y=maxy, label=pval), size=3, inherit.aes = F) 
      
    return(plot_compare)
}


boxdot2 = function(df, fs='cell', fn='source', g1='him', g2='merged non-him', fnrow=1, fncol=NULL, psize=1,ylabsize=1.2, pvalalter=T, amaxy=1){
    col_list = c(g1 = "#2b83ba", g2 = "#d7191c")
    names(col_list) <- c(g1, g2)
    # filter out x
    ind = df[, fn] %in% c(g1, g2)
    df = df[ind, ]
    featn =  colnames(df)[!colnames(df) %in% c(fs, fn)]  
    # cels
    cells = sort(unique(df[, fs]))
    # calculate pvals 
    pval = c()
    alter = c()
    for(ic in 1:length(cells)){
        dfc = df[df[, fs]==cells[ic], ]
        for(fi in 1:length(featn)){
            x = dfc[dfc[, fn] == g1, featn[fi]]
            y = dfc[dfc[, fn] == g2, featn[fi]]
            ptmp = compare2vect(x, y)
            pval = c(pval, ptmp$pval)
            alter = c(alter, ptmp$alter)
        }
    }
    miny = apply(df[, featn], 2, min, na.rm=T)
    maxy = apply(df[, featn], 2, max, na.rm=T)
    maxy = maxy + 0.05 * (maxy - miny)
    pval_table = data.frame(cell=rep(cells, rep(length(featn), length(cells))), feature=rep(featn, length(cells)), pval=pval, alter=alter, maxy=maxy )
    pval_tmp = pval_table$pval
    pval_table$pval = format(pval_table$pval, scientific=T, digits=3)
    pval_chr1 = rep('P<2.22e-16', nrow(pval_table))
    pval_chr2 = paste('P=', pval_table$pval, sep='')
    pval_table$pval = ifelse(pval_tmp<2.22e-16, pval_chr1, pval_chr2)
    if(pvalalter){
        pval_table$pval = paste(pval_table$pval, '\n', pval_table$alter, sep='')
    }
    #if(pvalalter){
    #    pval_table$pval = paste('P=', pval_table$pval, '\n', pval_table$alter, sep='')
    #} else {
    #    pval_table$pval = paste('P=', pval_table$pval, sep='')
    #}
    
    indna = grepl('NA', pval_table[, 'pval'])
    pval_table[indna, 'pval'] = ''
    
    df <- melt(df, id.vars = c(fs, fn), measure.vars = featn, variable.name = "feature", value.name="value")
    df[, fn] = factor(df[, fn], levels=c(g1, g2), ordered=T)
    df[, fs] = factor(df[, fs])
    
    plot_compare <- ggplot(df, aes_string(x = fn, y = "value", fill = fn)) +
      geom_jitter(alpha = 0.7, color="black", stroke = 0.0, shape=21, size=psize,  width = 0.2, height = 0.0) +
      geom_boxplot(data=df, aes_string(x = fn, y = "value", color = fn), outlier.shape = NA, width=0.6, color = "black", fill="white") +
      scale_fill_manual(name="",values = col_list, labels=rename_features(featn, nl)) +
      xlab("") +
      ylab("") +
      facet_grid(as.formula(paste('feature', '~', fs)), labeller=as_labeller(rename_features(featn, nl)), scales = "free_y", switch = 'y') +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
      theme(axis.title.y = element_text(size=10), axis.text.y=element_text(size=8)) + 
      theme(plot.margin = unit(c(0,0,0.1,-0.05), "cm"))
      plot_compare = plot_compare + theme(legend.position="bottom", legend.box.margin=margin(t=-20, r=0, b=-10, l=0), legend.text=element_text(size=10))
    #plot_compare = plot_compare + guides(colour=guide_legend(override.aes = list(size=3)))
      plot_compare <- plot_compare +
      geom_text(data=pval_table, aes(x=1.5, y=maxy*amaxy, label=pval), size=3, inherit.aes = F, vjust=1) 
    return(plot_compare)
}

# plot ppi comparison
boxdot_ppi = function(df, fm, fs){
    matrix_comp_0 <- melt(subset(df, TF_number >=0), id.vars = c('cell','source'), measure.vars = c("ppi_density"), variable.name = "feature", value.name = "value")
    matrix_comp_5 <- melt(subset(df, TF_number >=5), id.vars = c('cell','source'), measure.vars = c("ppi_density"), variable.name = "feature", value.name = "value")
    matrix_comp_10 <- melt(subset(df, TF_number >=10), id.vars = c('cell','source'), measure.vars = c("ppi_density"), variable.name = "feature", value.name = "value")
    #matrix_comp_15 <- melt(subset(df, TF_number >=15), id.vars = c('cell','source'), measure.vars = c("ppi_density"), variable.name = "feature", value.name = "value")
    matrix_comp_0$group <- "All"
    matrix_comp_5$group <- "#TF>=5"
    matrix_comp_10$group <- "#TF>=10"
    #matrix_comp_15$group <- "#TF>=15"
    #matrix_comp <- rbind(matrix_comp_0, matrix_comp_5, matrix_comp_10, matrix_comp_15)
    matrix_comp <- rbind(matrix_comp_0, matrix_comp_5, matrix_comp_10)
    #matrix_comp$group <- factor(matrix_comp$group, levels = c("All", "#TF>=5", "#TF>=10", "#TF>=15"))
    matrix_comp$group <- factor(matrix_comp$group, levels = c("All", "#TF>=5", "#TF>=10"))
    # get global mean
    global_mean = mean(unlist(subset(matrix_comp, source == "merged non-him")['value']))
    # gm12878 for the main figure
    plot_compare_ppi_v2 <- ggplot(subset(matrix_comp, source=='him' & cell=='gm12878'), aes(x = group, y = value)) +
    geom_jitter(alpha = 0.7, shape=21, color="black", fill = "#2b83ba", stroke = 0.0, width = 0.2, height = 0.0) +
    geom_boxplot(outlier.shape = NA, width=0.6, color = "black", fill = "white") +
    geom_hline(yintercept = global_mean, color = "#d7191c", size = 1) +
    coord_cartesian(ylim=c(0,0.65)) +
    xlab("") +
    ylab("PPI density") +
    theme_bw() +
    theme(plot.title = element_text(lineheight=.8, face="bold", size=rel(1.2))) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.text.y = element_text(color='black', size=7)) +
    theme(axis.text.x = element_text(color='black', size=8, angle=45, hjust=1, vjust=1)) +
    theme(axis.title.y = element_text(size=rel(1.2),margin=margin(0,10,0,0))) +
    theme(axis.title.x = element_text(size=rel(1.2),margin=margin(10,0,0,0))) +
    theme(axis.line = element_line(colour = "black")) +
    theme(strip.text.x = element_text(size=rel(1.2), face="bold")) +
    theme(strip.background = element_blank()) +
    theme(plot.margin = unit(c(0,0,-0.3,0.1), "cm"))
    
    pdf(file=fm, width=3, height=3)
    print(plot_compare_ppi_v2)
    dev.off()
        
    plot_compare_ppi_v3 <- ggplot(subset(matrix_comp, source=='him'), aes(x = group, y = value)) +
    geom_jitter(alpha = 0.7, shape=21, color="black", fill = "#2b83ba", stroke = 0.0, width = 0.2, height = 0.0) +
    geom_boxplot(outlier.shape = NA, width=0.6, color = "black", fill = "white") +
    geom_hline(yintercept = global_mean, color = "#d7191c", size = 1) +
    coord_cartesian(ylim=c(0,0.65)) +
    facet_grid(~cell, labeller = as_labeller(rename_features('ppi_density', nl))) + 
    xlab("") +
    ylab("PPI density") +
    theme_bw() +
    theme(plot.title = element_text(lineheight=.8, face="bold", size=rel(1.2))) +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.text.y = element_text(color='black', size=7)) +
    theme(axis.text.x = element_text(color='black', size=8, angle=45, hjust=1, vjust=1)) +
    theme(axis.title.y = element_text(size=rel(1.2),margin=margin(0,10,0,0))) +
    theme(axis.title.x = element_text(size=rel(1.2),margin=margin(10,0,0,0))) +
    theme(axis.line = element_line(colour = "black")) +
    theme(strip.text.x = element_text(size=rel(1.2), face="bold")) +
    theme(strip.background = element_blank()) + 
    theme(plot.margin = unit(c(0,0,-0.3,0.1), "cm"))
    
    pdf(file=fs, width=7.5, height=3)
    print(plot_compare_ppi_v3)
    dev.off()
}

boxdot3 = function(df, fs='cell', fn='source', g1='him', g2='non-him', g3='mhim', fnrow=1, fncol=NULL, psize=1,ylabsize=1.2, pvalalter=T, amaxy=1){
    col_list = c(g1 = "#d7191c", g2 = "#2b83ba", g3='#66c2a5')
    names(col_list) <- c(g1, g2, g3)
    # filter out x
    ind = df[, fn] %in% c(g1, g2, g3)
    df = df[ind, ]
    featn =  colnames(df)[!colnames(df) %in% c(fs, fn)]  
    # cels
    cells = sort(unique(df[, fs]))
    # calculate pvals 
    pval = c()
    alter = c()
    pval2 = c()
    alter2 = c()
    for(ic in 1:length(cells)){
        dfc = df[df[, fs]==cells[ic], ]
        for(fi in 1:length(featn)){
            x = dfc[dfc[, fn] == g2, featn[fi]]
            y = dfc[dfc[, fn] == g1, featn[fi]]
            z = dfc[dfc[, fn] == g3, featn[fi]]
            ptmp = compare2vect(x, y)
            pval = c(pval, ptmp$pval)
            alter = c(alter, ptmp$alter)
            ptmp2 = compare2vect(x, z)
            pval2 = c(pval2, ptmp2$pval)
            alter2 = c(alter2, ptmp2$alter)
        }
    }
    # hard coding switch < to >, vice versa in x vs y
    tmp = alter
    alter[tmp == '>'] = '<'
    alter[tmp == '<'] = '>'
    pval_table = data.frame(cell=rep(cells, rep(length(featn), length(cells))), feature=rep(featn, length(cells)), pval=pval, alter=alter,pval2=pval2, alter2=alter2 )
    pval_table$pval = format(pval_table$pval, scientific=T, digits=3)
    pval_table$pval2 = format(pval_table$pval2, scientific=T, digits=3)
    if(pvalalter){
        pval_table$pval = paste('P=\n', pval_table$pval, '\n', pval_table$alter, sep='')
        pval_table$pval2 = paste('P=\n', pval_table$pval2, '\n', pval_table$alter2, sep='')
    } else {
        pval_table$pval = paste('P=\n', pval_table$pval, sep='')
        pval_table$pval2 = paste('P=\n', pval_table$pval2, sep='')
    }
    indna = grepl('NA', pval_table[, 'pval'])
    pval_table[indna, 'pval'] = ''
    indna2 = grepl('NA', pval_table[, 'pval2'])
    pval_table[indna, 'pval2'] = ''
    
    df <- melt(df, id.vars = c(fs, fn), measure.vars = featn, variable.name = "feature", value.name="value")
    df[, fn] = factor(df[, fn], levels=c(g1, g2, g3), ordered=T)
    df[, fs] = factor(df[, fs])
    # manually remove outliers 
    if('gene_num' %in%featn){
        ind = df$feature == 'gene_num' & df$value >= 50
        df = df[!ind, ]
    }
    if('TF_number' %in% featn){
        ind = df$feature == 'TF_number' & df$value >= 75
        df = df[!ind, ]
    }
    miny = aggregate(x=df[, 'value'], by=list(df[, 'feature']), min, na.rm=T)[, 'x']
    maxy = aggregate(x=df[, 'value'], by=list(df[, 'feature']), max, na.rm=T)[, 'x']
    maxy = maxy + 0.05 * (maxy - miny)
    pval_table[, 'maxy'] = maxy
    
    plot_compare <- ggplot(df, aes_string(x = fn, y = "value", fill = fn)) +
      geom_jitter(alpha = 0.7, color="black", stroke = 0.0, shape=21, size=psize,  width = 0.2, height = 0.0) +
      geom_boxplot(data=df, aes_string(x = fn, y = "value", color = fn), outlier.shape = NA, width=0.6, color = "black", fill="white") +
      scale_fill_manual(name="",values = col_list, labels=rename_features(featn, nl)) +
      xlab("") +
      ylab("") +
      facet_grid(as.formula(paste('feature', '~', fs)), labeller=as_labeller(rename_features(featn, nl)), scales = "free_y", switch = 'y') +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
      theme(plot.margin = unit(c(0,0,0.1,-0.05), "cm"))
      plot_compare = plot_compare + theme(legend.position="bottom", legend.box.margin=margin(t=-20, r=0, b=-10, l=0))
    #plot_compare = plot_compare + guides(colour=guide_legend(override.aes = list(size=3)))
      plot_compare <- plot_compare +
      geom_text(data=pval_table, aes(x=1.5, y=maxy*amaxy, label=pval), size=1.5, inherit.aes = F, vjust=1) 
      plot_compare <- plot_compare +
      geom_text(data=pval_table, aes(x=2.5, y=maxy*amaxy, label=pval2), size=1.5, inherit.aes = F, vjust=1) 
    return(plot_compare)
}

extract_freq = function(df, cutoffs, labels){
    Apercent = c()
    cells = unique(X$cell)
    for( ic in 1:length(cells)){
        x = df[df$cell == cells[ic], 'A_percent']
        y=cut(x, breaks=cutoffs, labels=labels, include.lowest = T, right=F)
        z = as.data.frame(table(y))
        z[, 2] = z[, 2] / sum(z[, 2]) * 100
        z = data.frame(z, cells[ic], stringsAsFactors = F)
        Apercent = rbind(Apercent, z)
    }
    colnames(Apercent) = c('category', 'proportion', 'cell')
    return(Apercent)
}

boxviolin2 = function(df, fs='cell', fn='source', g1='him', g2='merged non-him', fnrow=1, fncol=NULL, psize=1,ylabsize=1.2, pvalalter=T, amaxy=1, facet.switch=NULL){
    col_list = c(g1 = "#2b83ba", g2 = "#d7191c")
    names(col_list) <- c(g1, g2)
    # filter out x
    ind = df[, fn] %in% c(g1, g2)
    df = df[ind, ]
    featn =  colnames(df)[!colnames(df) %in% c(fs, fn)]  
    # cels
    cells = sort(unique(df[, fs]))
    # calculate pvals 
    pval = c()
    alter = c()
    for(ic in 1:length(cells)){
        dfc = df[df[, fs]==cells[ic], ]
        for(fi in 1:length(featn)){
            x = dfc[dfc[, fn] == g1, featn[fi]]
            y = dfc[dfc[, fn] == g2, featn[fi]]
            ptmp = compare2vect(x, y)
            pval = c(pval, ptmp$pval)
            alter = c(alter, ptmp$alter)
        }
    }
    miny = apply(df[, featn], 2, min, na.rm=T)
    maxy = apply(df[, featn], 2, max, na.rm=T)
    maxy = maxy + 0.05 * (maxy - miny)
    pval_table = data.frame(cell=rep(cells, rep(length(featn), length(cells))), feature=rep(featn, length(cells)), pval=pval, alter=alter, maxy=maxy )
    indPval_small = pval_table$pval < 2.22e-16 & !is.na(pval_table$pval)
    pval_table$pval = format(pval_table$pval, scientific=T, digits=3)
    pval_table[!indPval_small, 'pval'] = paste('P=', pval_table[!indPval_small, 'pval'], sep='')
    pval_table[indPval_small, 'pval'] = 'P<2.22e-16'
    if(pvalalter){
        pval_table$pval = paste(pval_table$pval, '\n', pval_table$alter, sep='')
    }
    indna = grepl('NA', pval_table[, 'pval'])
    pval_table[indna, 'pval'] = ''
    
    df <- melt(df, id.vars = c(fs, fn), measure.vars = featn, variable.name = "feature", value.name="value")
    df[, fn] = factor(df[, fn], levels=c(g1, g2), ordered=T)
    df[, fs] = factor(df[, fs])
    
    plot_compare <- ggplot(df, aes_string(x = fn, y = "value", color = fn)) +
      geom_violin() + geom_boxplot(width=0.05, outlier.shape = NA) + 
      xlab("") + ylab("") +
      facet_grid(as.formula(paste('feature', '~', fs)), labeller=as_labeller(rename_features(featn, nl)), scales = "free_y", switch=facet.switch) +
      scale_color_manual(name = 'Gene assignment to HIMs', values = col_list, labels=rename_features(featn, nl) ) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
      theme(axis.title.y = element_text(size=10), axis.text.y=element_text(size=8)) + 
      theme(plot.margin = unit(c(0,0,0.1,-0.05), "cm"))
      plot_compare = plot_compare + theme(legend.position="bottom", legend.box.margin=margin(t=-20, r=0, b=-10, l=0),
                                      legend.text=element_text(size=10))
      plot_compare <- plot_compare +
      geom_text(data=pval_table, aes(x=1.5, y=maxy*amaxy, label=pval), size=3, inherit.aes = F, vjust=1) 
    return(plot_compare)
}

boxviolin = function(df, fn='source', g1='him', g2='merged non-him', cell='F!', fnrow=1, fncol=NULL, psize=1, ylabsize=1.2, pvalalter=T){
    col_list = c(g1 = "#2b83ba", g2 = "#d7191c")
    names(col_list) <- c(g1, g2)
    # filter out x
    ind = df[, fn] %in% c(g1, g2)
    df = df[ind, ]
    featn =  colnames(df)[colnames(df)!=fn]  
    # calculate pvals 
    pval = c()
    alter = c()
    for(fi in 1:length(featn)){
        x = df[df[, fn] == g1, featn[fi]]
        y = df[df[, fn] == g2, featn[fi]]
        ptmp = compare2vect(x, y)
        pval[fi] = ptmp$pval
        alter[fi] = ptmp$alter
    }
    miny = apply(df[, featn], 2, min, na.rm=T)
    maxy = apply(df[, featn], 2, max, na.rm=T)
    maxy = maxy + 0.05 * (maxy - miny)
    pval_table = data.frame(feature=featn, pval=pval, alter=alter, maxy=maxy )
    indPval_small = pval_table$pval < 2.22e-16 & !is.na(pval_table$pval)
    pval_table$pval = format(pval_table$pval, scientific=T, digits=3)
    pval_table[!indPval_small, 'pval'] = paste('P=', pval_table[!indPval_small, 'pval'], sep='')
    pval_table[indPval_small, 'pval'] = 'P<2.22e-16'
    if(pvalalter){
        pval_table$pval = paste(pval_table$pval, '\n', pval_table$alter, sep='')
    }
    indna = grepl('NA', pval_table[, 'pval'])
    pval_table[indna, 'pval'] = ''
    
    df <- melt(df, id.vars = fn, measure.vars = featn, variable.name = "feature", value.name="value")
    df[, fn] = factor(df[, fn], levels=c(g1, g2), ordered=T)
    # calculate pvals
    #pval_table = data.frame(feature=as.factor(featn),  pval_single=letters[1:4], x=rep(1.5, 4), y=rep(0.1, 4), value=fessg)
    
    plot_compare <- ggplot(df, aes_string(x = fn, y = "value", color = fn)) +
        geom_violin() + geom_boxplot(width=0.05, outlier.shape = NA) + 
       # geom_jitter(alpha = 0.7, color="black", stroke = 0.0, shape=21, size=psize,  width = 0.2, height = 0.0) +
        #geom_boxplot(data=df, aes_string(x = fn, y = "value", color = fn), outlier.shape = NA, width=0.6, color = "black", fill="white") +
        scale_color_manual(name="",values = col_list, labels=rename_features(featn, nl)) +
        xlab(label='') +
        ylab("") +
        facet_wrap(~feature, labeller=as_labeller(rename_features(featn, nl)), nrow = fnrow, ncol=fncol, scales = "free_y", strip.position = 'left', shrink=T) +
        theme_bw() +
        theme(plot.title = element_text(lineheight=.8, face="bold", size=rel(ylabsize))) +
        theme(panel.spacing = unit(0.5, 'lines')) + 
        theme(panel.border = element_blank()) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(axis.title.y = element_text(size=rel(ylabsize * 0.9),margin=margin(0,1,0,0))) +
        theme(axis.title.x = element_text(size=rel(1.2),margin=margin(1,0,0,0))) +
        theme(axis.text.y = element_text(color='black', size=7)) +
        theme(axis.line = element_line(color = "black")) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        theme(strip.text.x = element_text(size=rel(1.2), face="bold")) +
        theme(strip.text.y = element_text(size=rel(ylabsize))) +
        theme(strip.background = element_blank())  +
        theme(strip.placement = 'outside')  + 
        theme(plot.margin = unit(c(0,0,0.1,-0.5), "cm"))
    plot_compare = plot_compare + theme(legend.position="bottom", legend.box.margin=margin(t=-22, r=0, b=-10, l=0))
    plot_compare <- plot_compare +
        #geom_text(data=pval_table, aes(x=1.5, y=maxy, label=pval), size=3, position=position_dodge(0.9), inherit.aes = F) 
        geom_text(data=pval_table, aes(x=1.5, y=maxy, label=pval), size=3, inherit.aes = F) 
    return(plot_compare)
}

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}