extract_gene_tf_chr_him  = function(){
    hims = read.table('data/him_summary_allinone.txt', header=T, sep='\t', stringsAsFactors = F)
    hims = hims[hims$source == 'him',]
    hims_ids = rownames(hims)
    hims_chr = hims[, 'chrom']; names(hims_chr) = hims_ids
    hims_genes = list(); hims_tfs = list(); 
    for(i in 1:nrow(hims)){
        tmp_gene =  unlist(strsplit(hims[i, 'genes'], split=','))
        # A;B genes is treated as A as there is no big difference
        tmp_gene = gsub(';.*', '', tmp_gene)
        hims_genes[[ hims_ids[i] ]] = tmp_gene 
        hims_tfs[[ hims_ids[i] ]] =  unlist(strsplit(hims[i, 'TFs'], split=','))
    }
    result = list(genes=hims_genes, tfs=hims_tfs, chr=hims_chr, himids = hims_ids)
    return(result)
}
background_genes = function(){
    gene_bed = read.table('data/gene_chrom_bin_num_hg19_combined_sorted.bed', header=F, sep='\t', stringsAsFactors = F)
    rownames(gene_bed) = gene_bed[, 4]
    chrs = unique(gene_bed[,1]); chrs = chrs[chrs!='chrY']
    gene_bed_chr = lapply(chrs, function(z) gene_bed[gene_bed[, 1] == z, 4])
    names(gene_bed_chr) = chrs
    return(gene_bed_chr)
}
load_him_pairs_JI = function(){
    hp = read.table('data/him_dynamic_allinone.txt', header=T, stringsAsFactors = F, sep='\t')
    cn = c('himid1','cell1', 'himid2', 'cell2', 'jiTF', 'jiGene')
    hp = hp[, cn]
    return(hp)
}

#compute_ji_pval = function(){
#    library(doParallel)
#    registerDoParallel(cores=6)
#    system.time({ res = foreach(i=1:nrow(hims_ji), .combine=rbind) %dopar% {
#    #for(i in 1:500){
#        print(i)
#        himid1 = hims_ji[i, 'himid1']
#        himid2 = hims_ji[i, 'himid2']
#        chr = himinfo[['chr']][himid1]
#        gb = gene_bg_chr[[chr]]
#        xg = gb %in% himinfo[['genes']][[himid1]]
#        yg = gb %in% himinfo[['genes']][[himid2]]
#        ji_gene = jaccard(x=xg, y = yg)
#        ji_gene_test = unlist(jaccard.test(x=xg, y=yg, method='bootstrap',fix='x', verbose=F)[c('statistics', 'pvalue', 'expectation')])
#        xt = TFs_bg %in% himinfo[['tfs']][[himid1]]
#        yt = TFs_bg %in% himinfo[['tfs']][[himid2]]
#        ji_tf = jaccard(x=xt, y=yt)
#        ji_tf_test = unlist(jaccard.test(x=xt, y=yt, method='bootstrap', fix='x', verbose=F)[c('statistics', 'pvalue', 'expectation')])
#        tmp = c(ji=ji_gene, ji_gene_test, ji=ji_tf, ji_tf_test)
#        names(tmp) = paste(names(tmp), rep(c('_gene','_tf'), rep(4, 2)), sep='')
#        return(tmp)
#       #res = rbind(res, tmp)
#    #}
#    } })
#    return(res)
#}

hypergeom_vectors = function(v1, v2){
    # v2 is boolean vector of white balls
    # v3 is boolean vector of draw balls
    q = sum(v1 & v2)
    m = sum(v1)
    n = length(v1) - m
    k = sum(v2)
    pval = phyper(q=q, m=m, n=n, k=k, lower.tail = F, log.p = FALSE)
    names(pval) = 'pvalue'
    return(pval)
}

compute_ji_pval_hypergeo = function(){
    library(doParallel)
    registerDoParallel(cores=6)
    system.time({ res = foreach(i=1:nrow(hims_ji), .combine=rbind) %dopar% {
    #for(i in 1:500){
        print(i)
        himid1 = hims_ji[i, 'himid1']
        himid2 = hims_ji[i, 'himid2']
        chr = himinfo[['chr']][himid1]
        gb = gene_bg_chr[[chr]]
        xg = gb %in% himinfo[['genes']][[himid1]]
        yg = gb %in% himinfo[['genes']][[himid2]]
        pval1 = hypergeom_vectors(v1=xg, v2=yg)
        #ji_gene = jaccard(x=xg, y = yg)
        #ji_gene_test = unlist(jaccard.test(x=xg, y=yg, method='bootstrap',fix='x', verbose=F)[c('statistics', 'pvalue', 'expectation')])
        xt = TFs_bg %in% himinfo[['tfs']][[himid1]]
        yt = TFs_bg %in% himinfo[['tfs']][[himid2]]
        pval2 = hypergeom_vectors(v1=xt, v2=yt)
        #ji_tf = jaccard(x=xt, y=yt)
        #ji_tf_test = unlist(jaccard.test(x=xt, y=yt, method='bootstrap', fix='x', verbose=F)[c('statistics', 'pvalue', 'expectation')])
        tmp = c(pval1, pval2)
        names(tmp) = paste(names(tmp), c('_gene','_tf'), sep='')
        return(tmp)
       #res = rbind(res, tmp)
    #}
    } })
    return(res)
}

hims_stable_cs = function(){
    # stalbe hims are HIMs share significant number of genes with hims in at least one other cell types
    # the others are defined as cell type-specific HIMs
    ind_pval = hims_dyna$pvalue_gene <= 0.001
    #ind_fc = hims_dyna$jiGene >= 2 * hims_dyna$expectation_gene
    ind_abs = hims_dyna$jiGene >= 1/3
    table(ind_pval, ind_abs)
    fivenum(hims_dyna[ind_pval, 'jiGene'])
    hims_dyna_sub = hims_dyna[ind_pval & ind_abs, ]
    hims_dyna_sub = hims_dyna[ind_abs, ]
    hims_dyna_sub = hims_dyna[ind_pval, ]
    # stable hims 
    stable_hims = aggregate(himid1 ~ cell1, hims_dyna_sub, FUN=unique) 
    sapply(stable_hims$himid1, length)
    # test on TF set
    ind_pval = hims_dyna$pvalue_tf <= 0.001
    #ind_fc = hims_dyna$jiGene >= 2 * hims_dyna$expectation_gene
    ind_abs = hims_dyna$jiTF >= 1/3
    table(ind_pval, ind_abs)
    fivenum(hims_dyna[ind_pval, 'jiTF'])
    hims_dyna_sub = hims_dyna[ind_pval & ind_abs, ]
    hims_dyna_sub = hims_dyna[ind_abs, ]
    hims_dyna_sub = hims_dyna[ind_pval, ]
    # stable hims 
    stable_hims = aggregate(himid1 ~ cell1, hims_dyna_sub, FUN=unique) 
    sapply(stable_hims$himid1, length)
}