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

compute_ji_pval = function(){
    res = c()
    for(i in 1:nrow(hims_ji)){
        print(i)
        himid1 = hims_ji[i, 'himid1']
        himid2 = hims_ji[i, 'himid2']
        chr = himinfo[['chr']][himid1]
        gb = gene_bg_chr[[chr]]
        xg = gb %in% himinfo[['genes']][[himid1]]
        yg = gb %in% himinfo[['genes']][[himid2]]
        ji_gene = jaccard(x=xg, y = yg)
        ji_gene_test = unlist(jaccard.test(x=xg, y=yg, method='exact'))
        xt = TFs_bg %in% himinfo[['tfs']][[himid1]]
        yt = TFs_bg %in% himinfo[['tfs']][[himid2]]
        #ji_tf = jaccard(x=xt, y=yt)
        #ji_tf_test = unlist(jaccard.test(x=xt, y=yt, method='exact'))
        #tmp = c(ji=ji_gene, ji_gene_test, ji=ji_tf, ji_tf_test)
        #names(tmp) = paste(names(tmp), rep(c('_gene','_tf'), rep(4, 2)), sep='')
        #res = rbind(res, tmp)
    }
    return(res)
}
