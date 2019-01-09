library(splitstackshape)
library(ggplot2)
library(reshape2)
library(doParallel)

# fix the peaks and gene set; randomly pick new gene set while preserving the mean distance and sd 
# step 1. load the peak count near the tss table 
load_peak_table = function(cell, w, filter=T, chip_coverage=0.1){
    path = paste('data/chip-peaks/gene_TF_peak_combined_',cell,'_w_', w, '.txt', sep='')
    peak = read.table(path, header = T, row.names = 1, sep='\t')
    peak_coverage = apply(peak, 2, function(z) sum(z>0) / length(z))
    peak = peak[, peak_coverage>= chip_coverage]
    return(peak)
}
# step 2. load the genomic locations of genes in the GRN of a cell type
load_gene_grn_bed = function(cell){
    genes = read.table(paste('data/gene_symbol_list_', cell,'.txt', sep=''), stringsAsFactors = F)
    genes_bed = read.table('data/gene_chrom_bin_num_hg19_combined_sorted.bed', header=F, sep='\t', stringsAsFactors = F)
    genes_bed = genes_bed[genes_bed[, 4] %in% genes[, 1], ]
    return(genes_bed)
}
# step 3. compute the pairwise distance between genes per chromosome 
compute_pairwise_distance = function(df){
    # df is bed file 
    chroms = unique(df[, 1])
    chroms = chroms[chroms!= 'chrY']
    result = list()
    for(chr in chroms){
       df.chr = df[df[, 1] == chr, ] 
       tss = df.chr[, 2]
       names(tss) = df.chr[, 4]
       mat.dist = dist(tss, method='manhattan', upper=T)
       mat.dist = as.matrix(mat.dist)
       result[[chr]] = mat.dist
    }
    return(result)
}
# step 4. find the him genes per TFs
# step 4. find the chromosome of each him
extract_him_genes_per_TF = function(i, cell, tfs){
    hims = read.table(i, header=T, sep='\t', stringsAsFactors = F)
    hims = hims[hims$source == 'him', ]
    hims = hims[grepl(cell, rownames(hims)), ]
    hims = hims[, c('TFs', 'genes', 'chrom')]
    hims = data.frame(him_id = rownames(hims), hims, stringsAsFactors = F)
    result = list()
    result_chr = list()
    for(tf in tfs){
        ind = grepl(tf, hims$TFs)
        hims.sub = hims[ind, ]
        res2 = lapply(hims.sub[, 'genes'], function(z) unlist(strsplit(z, split = ';|,')))
        names(res2) = hims.sub[, 'him_id']
        result[[tf]] = res2
        
        tmp_chr = hims.sub[, 'chrom']
        names(tmp_chr) = hims.sub[, 'him_id']
        result_chr[[tf]] = tmp_chr
    }
    result_comb = list(genes=result, chr=result_chr)
    return(result_comb)
}
# step 5. generate new random gene sets
# this is vital, double check!
random_geneset_oneset = function(ngene, mat, mean.const, sigma.const, max.const){
    n1 = nrow(mat)
    genes = colnames(mat)
    # step 1. randomly pick two  suitable farthest apart genes 
    ind.max.iter = T
    max.try = 0
    while(ind.max.iter & max.try < 100){
    max.try = max.try + 1
    farest_dist = 0; ngene_cand = 0;
    g1.cand.index = 1:n1
    while(( (abs(farest_dist - max.const) >= 0.2*max.const) | (ngene_cand < ngene)) & length(g1.cand.index) > 0){
        g1.index = sample(g1.cand.index, 1)
        g1.cand.index = g1.cand.index[g1.cand.index != g1.index]
        g1 = genes[g1.index]
        g1.dist = mat[g1.index, ]
        g2.index = which.min(abs(g1.dist - max.const))
        g2 = genes[g2.index]
        farest_dist = mat[g1.index, g2.index]
        ngene_cand = abs(g1.index - g2.index)
    }
    if(length(g1.cand.index) > 0){
        if(g1.index > g2.index){
            tmp = g1.index
            g1.index = g2.index
            g2.index = tmp
        }
        # step 2; pick the rest ngene -2 genes
        # try brutial force 
        k = 1; 
        k.max = min(100, choose(ngene_cand, ngene-2))
        mean_initial = 0; sigma_initial=0
        g.index = c()
        while( (k<k.max) & (abs(mean.const - mean_initial)/mean.const >=0.2 | abs(sigma.const - sigma_initial)/sigma.const>=0.2)){
            k = k + 1
            gr.index = sample(g1.index:g2.index, ngene-2, replace = F)
            g.index = c(g1.index, g2.index, gr.index)
            mat_tmp = mat[g.index, g.index]
            mean_initial = mean(mat_tmp)
            sigma_initial = sd(mat_tmp)
        }
        
        if(k==k.max){
            ind.max.iter = T
        } else {
            ind.max.iter = F
        }
    
    } else {
        genesets=NULL
    }
    
    }
    if(max.try < 100){
        geneset = genes[g.index]
    } else {
        geneset = NULL
    }
    return(geneset)
}

random_geneset = function(genes, dist.mat, N, parallel=T){
    n = length(genes)
    ind = rownames(dist.mat) %in% genes 
    dist.mat.real = dist.mat[ind, ind]
    dist.mat.cad = dist.mat[!ind, !ind]
    
    # mean, sd, max are used as constraint
    dist.mean = mean(dist.mat.real)
    dist.sd = sd(dist.mat.real)
    dist.max = max(dist.mat.real)
    ### potentially parallel implementation
    ## parallel version 
    if(parallel){
        result = foreach(i=1:N, .combine=rbind) %dopar% random_geneset_oneset(ngene=n, mat=dist.mat.cad, mean.const = dist.mean, sigma.const = dist.sd, max.const=dist.max)
    } else {
        #result = matrix('', nrow=N, ncol=n)
        result = c()
        k.empty = 0
        for(i in 1:N){
            #if(i %% 100 == 0){
            #    print(paste(i, ' samples are geenrated'))
            #}
            tmp = random_geneset_oneset(ngene=n, mat=dist.mat.cad, mean.const = dist.mean, sigma.const = dist.sd, max.const=dist.max)
            if(! is.null(tmp)){
                #result[i, ] = sort(tmp)
                result = rbind(result, sort(tmp))
            } else {
                k.empty = k.empty + 1
            }
            
            if(k.empty >=50){
               result = NULL
               break 
            }
        }
    }
    return(result)
}

random_geneset_output = function(cell, N){
    peak_gene = load_peak_table(cell=cell, w='5000', filter=F, chip_coverage=0.05)
    tfs = colnames(peak_gene)
    genes_bed = load_gene_grn_bed(cell=cell)
    gene_dist = compute_pairwise_distance(genes_bed)
    gs_tf_chr = extract_him_genes_per_TF(i='data/him_summary_allinone.txt', cell=cell, tfs=tfs)
    gs_tf = gs_tf_chr[['genes']]; gs_chr = gs_tf_chr[['chr']]
    
    ncomb = c()
    for(i in 1:length(tfs)){
        nh = length(gs_tf[[i]])
        tmp = cbind(i, 1:nh)
        ncomb = rbind(ncomb, tmp)
    }
    res = c()
    for(k in 1:nrow(ncomb)){
        i = ncomb[k, 1]; j = ncomb[k, 2]
        tf = tfs[i]
        him = names(gs_tf[[i]])[j]
        n1 = length(gs_tf[[i]][[j]])
        tmp = paste(gs_tf[[i]][[j]], collapse = ';')
        tmp = data.frame(tf=tf, him=him, genes=tmp, status='real', n=n1, stringsAsFactors = F)
        res = rbind(res, tmp)
        permu_res = random_geneset(genes=gs_tf[[i]][[j]], dist.mat=gene_dist[[ gs_chr[[i]][[j]] ]], N=N, parallel=F)
        if(!is.null(permu_res) & length(permu_res)>0 ){
            if(is.null(dim(permu_res))){
                permu_res = t(data.frame(permu_res, stringsAsFactors = F))
            } else {
                permu_res = permu_res[!duplicated(permu_res), ]
            }
            if( nrow(permu_res) > 0.1 * N){
                permu_res_char = apply(permu_res, 1, function(z) paste(z, collapse = ';'))
                tmp = data.frame(tf=tf, him=him, genes=permu_res_char, status='random', n=ncol(permu_res), stringsAsFactors = F)
                res = rbind(res, tmp)
            }
        }
    }
    ofile = paste('inter_results/', 'random_genes_', cell, '.txt', sep='')
    write.table(res, ofile, col.names = T, row.names = F, sep='\t', quote=F)
    
}

pval_norm = function(df){
    ind_rand = df[, 1] == 'random'
    x = df[!ind_rand, 2]
    y = df[ind_rand, 2]
    pval = sum(y>=x) / length(y) 
    pval2 = pnorm(x, mean(y), sd(y), lower.tail = F)
    tmp = c(x, mean(y), pval, pval2)
    print(tmp)
    return(tmp)
}
random_geneset_comppute_pval = function(cell, w, binary=T){
    peak_gene = load_peak_table(cell=cell, w=w, filter=F, chip_coverage=0.05)
    if(binary){
        peak_gene = peak_gene >= 1
    }
    ifile = paste('inter_results/', 'random_genes_', cell, '.txt', sep='')
    rand_genes = read.table(ifile, header=T, sep='\t', stringsAsFactors = F) 
    npeak = c()
    for(i in 1:nrow(rand_genes)){
        tf = rand_genes[i, 'tf']
        genes = unlist(strsplit(rand_genes[i, 'genes'], split=';'))
        npeak[i] = sum(peak_gene[genes, tf])
    }
    rand_genes = cbind(rand_genes, npeak=npeak)
    # compute pval by norm distribution
    res = by(rand_genes[, c('status', 'npeak')], list(rand_genes$tf, rand_genes$him), FUN=pval_norm)
    print(str(res))
    stop()
    return(rand_genes)
} 

random_geneset_pval = function(i, j, gs_tf, gs_chr, peak_gene, gene_dist, N, parallel=T){
    print(c(i, j))
    permu_res = random_geneset(genes=gs_tf[[i]][[j]],  dist.mat=gene_dist[[ gs_chr[[i]][[j]] ]], N=N, parallel=parallel)
    if(!is.null(permu_res) & length(permu_res)>0 ){
        permu_res = permu_res[!duplicated(permu_res), ]
        if( nrow(permu_res) > 0.1 * N){
            # double check that there are no to litter overlaps between rows
            #ind = duplicated(permu_res)
            #if(sum(ind)> 0.1 * nrow(permu_res)){
            #    print(N)
            #    print(sum(ind))
            #   stop('Too many duplicated random sets') 
            #}
            
            npeak = sum(peak_gene[ gs_tf[[i]][[j]], tfs[i] ])
            npeak_random = apply(permu_res, 1, function(z) sum(peak_gene[z, tfs[i]]) )
            pval = sum(npeak_random >= npeak) / N
            pval2 = pnorm(npeak, mean(npeak_random), sd(npeak_random), lower.tail = F)
            tmp1 = c(length(gs_tf[[i]][[j]]), npeak, pval, pval2)
            tmp1 = as.character(tmp1)
            tmp = c(tfs[i], names(gs_tf[[i]])[j], tmp1)
        }
    } else {
        tmp= NULL
    }
    return(tmp)
}

# analysis the resutled enrichment results
library(ggplot2)
summarize_enrich = function(cell, w){
    ifile_master = paste('data/master_tfs_subset_grn_', cell, '.txt', sep='')
    tf_master = read.table(ifile_master, header=F, stringsAsFactors = F)[, 1]
    ifile = paste('inter_results/chip_enrichment_', cell, '_w_', w, '.txt',sep='')
    res = read.table(ifile, header = F, sep='\t', stringsAsFactors = F)
    colnames(res) = c('tf', 'him', 'ngene', 'ngene_w_peak', 'pval', 'pval_norm', 'prop')
    #### summarize per TF
    res_pval = aggregate(pval~tf, data=res, FUN=function(z) sum(z<=0.05)/length(z))
    res_prop = aggregate(prop~tf, data=res, FUN=function(z) sum(z>=0.5)/length(z))
    res2 = data.frame(res_pval, prop=res_prop$prop, stringsAsFactors = F)
    res2 = data.frame(res2, master=as.integer(res2$tf %in% tf_master), stringsAsFactors = F)
    res2 = res2[order(res2$master, res2$prop), ]
    print(res2)
    #### summarize by comparing master TFs and non-master TFs
    #res[, 'master'] = as.integer(res[, 'tf'] %in% tf_master)
    #res_tab = table(res[,'pval_norm']<=0.05, res[, 'prop']>=50, res[,'master'])
    #print(res_tab)
    #print(prop.table(res_tab))
    #####
    #res[, 'pval'] = -1*log10(res[, 'pval'])
    #res[, 'pval_norm'] = -1*log10(res[, 'pval_norm'])
    #p = ggplot(res, aes(x=prop, y=pval_norm, color=tf)) + geom_point() + geom_hline(yintercept = -1*log10(0.05))
    #print(p)
    return(res)
}
## observations 
# pval has weak correlations with prop
# 40% HIMs regualted by master TFs have master TFs'peaks in majority of HIMs genes 
# which is expected higher than non-master TFs because master TFs have higher number of peaks 
# still only small proportion of HIMs have pval <= 0.05; again it is not a good idea to test when sample size is small
# even HIMs with 100% genes having peaks are still non-sig (p~0.3 for master TF RUNX3)

######## method overhaul again 
# step1. find the him gnes 
gs_from_hims = function(i, cell){
    # load genes in the GRN of the cell type
    genes = read.table(paste('data/gene_symbol_list_', cell,'.txt', sep=''), stringsAsFactors = F)
    # load genes and TFs in each HIM
    hims = read.table(i, header=T, sep='\t', stringsAsFactors = F)
    hims = hims[hims$source == 'him', ]
    hims = hims[grepl(cell, rownames(hims)), ]
    hims = hims[, c('TFs', 'genes')]
    hims = data.frame(gs_id = rownames(hims), hims, stringsAsFactors = F)
    hims = cSplit(hims, 'genes', sep=';', direction='long')
    hims = cSplit(hims, 'genes', sep=',', direction='long')
    hims = as.data.frame(hims)
    hims$genes = as.character(hims[, 'genes'])
    # find the TF list
    TFs = unique(unlist(sapply(hims$TFs, function(z) unlist(strsplit(z, ";|,")))))
    result = list()
    for(i in TFs){
        ind = grepl(i, hims$TFs)
        if(sum(ind)>5){
            tmp = hims[ind, ]
            result[[i]] = tmp
        }
    }
    return(result)
}
#### load the summarize table for peaks near TSS
#load_peak_table = function(cell, w){
#    path = paste('data/chip-peaks/gene_TF_peak_combined_',cell,'_w_', w, '.txt', sep='')
#    peak = read.table(path, header = T, row.names = 1, sep='\t')
#    return(peak)
#}
#### compute enrichment per him per TF, using fisher exact test 
## output; ngene, ngene w peak, prop,  porp of peaks / gene, pvalue
chipenrich_him = function(df1, df2, method.test='individual'){
    # df1 saves the hims and their geens that are regulated by a TF
    # df2 saves the nunmber of TF peaks near the TSS of each gene
    #propreocess, df2 to yes or no
    df2 = df2 >= 1
    p_back = sum(df2) / nrow(df2)
    if(method.test == 'individual'){
        him_u = unique(df1[, 'gs_id'])
        res = c()
        for(i in 1:length(him_u)){
            genes = df1[df1[,'gs_id'] == him_u[i], 'genes']
            n = length(genes)
            gene_bool = rownames(df2) %in% genes
            n1 = sum(df2[gene_bool, ])
            p_him = n1 / n
            pval = binom.test(x=n1, n=n, p=p_back, alternative = 'g')$p.value
            tmp = c(n, n1, pval, p_him)
            res = rbind(res, tmp)
        }
        res = cbind(res, p_back)
    } else if(method.test=='population'){
        him_u = 'Population'
        genes = df1[, 'genes']
        n = length(genes)
        gene_bool = rownames(df2) %in% genes
        n1 = sum(df2[gene_bool, ])
        p_him = n1 / n
        pval = binom.test(x=n1, n=n, p=p_back, alternative = 'g')$p.value
        res = c(n, n1, pval, p_him, p_back)
        res = t(as.data.frame(res))
    }
    colnames(res) = c('N.gene', 'N.gene.w.peak', 'pval', 'Prop', 'Prop.background')
    res = data.frame(him=him_u, res, stringsAsFactors = F)
    return(res)
}

chipenrich_him_combine = function(tfs, hims_gs, peak_gene, method.test='individual'){
    res = c()
    for(i in 1:length(tfs)){
        tmp = chipenrich_him(df1=hims_gs[[tfs[i]]], df2=subset(peak_gene, select = tfs[i]), method.test=method.test)
        tmp = data.frame(tf=tfs[i],tmp, stringsAsFactors = F)
        res = rbind(res, tmp)
    }
    return(res)
}
chipenrich_result_summarize = function(df, method.test='individual', tf_group=F){
    sum.method = c('prop', 'pval', 'background')
    ind = data.frame(prop=df$Prop >= 0.8, pval=df$pval <= 0.05, background=df$Prop > df$Prop.background)
    ind = apply(ind, 2, as.integer)
    df = data.frame(df, ind, stringsAsFactors = F)
    return(df)
}
chipenrich_him_main = function(him_gs, cell, w, method.test, tf_group=F, chip_coverage=0.05){
    peak_gene = load_peak_table(cell=cell, w=w)
    ## filter out TFs with <= 5% genes having their peaks; these TFs might not co-operate with other TFs to regulate multiple genes
    peak_coverage = apply(peak_gene, 2, function(z) sum(z>0) / length(z))
    peak_gene = peak_gene[, peak_coverage>= chip_coverage]
    tfs_wc = colnames(peak_gene)
    res = chipenrich_him_combine(tfs=tfs_wc, hims_gs=hims_gs, peak_gene = peak_gene, method.test=method.test)
    res.sum = chipenrich_result_summarize(df=res, method.test=method.test, tf_group=tf_group)
    return(res.sum)
}
#### method  overhaul stops here 

#library(chipenrich)
# it is a model-based method to test enrichment of TF Chip-seq peaks in a given gene sets
# it is a competive model which compare the gene set and the other genes
#gs_from_hims = function(i, cell){
#    # load genes in the GRN of the cell type
#    genes = read.table(paste('data/gene_symbol_list_', cell,'.txt', sep=''), stringsAsFactors = F)
#    # map genes to Entrez
#    library("org.Hs.eg.db")
#    ENTREZ <- as.list(org.Hs.egALIAS2EG); ENTREZ <- ENTREZ[!is.na(ENTREZ)]
#    sym2entre = ENTREZ[genes[,1]]; sym2entre = unlist(sapply(sym2entre, function(z) z[1]))
#    sym2entre_u = unique(sym2entre)
#    
#    hims = read.table(i, header=T, sep='\t', stringsAsFactors = F)
#    hims = hims[hims$source == 'him', ]
#    hims = hims[grepl(cell, rownames(hims)), ]
#    hims = hims[, c('index', 'TFs', 'genes')]
#    hims = data.frame(gs_id = rownames(hims), hims, stringsAsFactors = F)
#    hims = cSplit(hims, 'genes', sep=';', direction='long')
#    hims = cSplit(hims, 'genes', sep=',', direction='long')
#    hims$genes = as.character(as.vector(hims[, 'genes'][[1]]))
#    # remove genes that do not have Entrez ID 
#    hims = hims[hims$genes %in% names(sym2entre), ]
#    hims$entrez = sym2entre[hims$genes]
#    # find the TF list
#    TFs = unique(unlist(sapply(hims$TFs, function(z) unlist(strsplit(z, ";|,")))))
#    result_comb = list()
#    for(i in TFs){
#        ind = grepl(i, hims$TFs)
#        if(sum(ind)>10){
#            entrez_i = hims[ind, 'entrez'][[1]]
#            entrez_i = unique(entrez_i)
#            entrez_io = sym2entre_u[!sym2entre_u %in% entrez_i]
#            tmp = data.frame(gs_id=i, entrez=entrez_i, stringsAsFactors = F)
#            tmp2 = data.frame(gs_id = 'Other', entrez=entrez_io, stringsAsFactors = F)
#            result_comb[[i]] = rbind(tmp, tmp2)
#        }
#    }
#    result = list(comb=result_comb)
#    return(result)
#}

#TFs_w_chip_seq = function(cell){
#    meta = read.table(paste('data/', cell, '/metadata_slim_from_all.tsv', sep=''), header=F, stringsAsFactors = F, sep='\t')
#    meta = meta[, c(1, 13)]
#    meta[, 'filepath'] = paste('data/', cell, '/', meta[,2], '_', meta[,1], '.bed.gz', sep='')
#    meta = meta[, c(2, 3)]
#    colnames(meta) = c('TF', 'filepath')
#    meta$TF = gsub('-human','', meta$TF)
#    rownames(meta) = meta$TF
#    return(meta)
#}
# copy bed files from the server in case they do not exist in the local directory
#scp -r cmu:/hive/dechaot/Data/ENCODE/chip-seq/gm12878 data/
#chip_enrichment = function(cell){
#    meta_tf = TFs_w_chip_seq(cell=cell)
#    hims_gs = gs_from_hims('data/him_summary_allinone.txt', cell=cell)
#    l1 = names(hims_gs[['comb']]); l2 = meta_tf$TF
#    # some TFs with Chip-seq data are not in the GRNs such as JUND
#    meta_tf = meta_tf[l2 %in% l1, ]
#    meta_tf['NFE2', ]
#    pval = list()
#    i = which(meta_tf$TF == 'NFE2')
#    for(i in 1:nrow(meta_tf)){
#        tf = meta_tf[i, 1]
#        gs_path = paste('inter_results/gs_', tf,'.txt', sep='')
#        write.table(hims_gs[['comb']][[tf]], file=gs_path, col.names = T, row.names = F, sep='\t', quote=F)
#        peaks_tf = read.table(gzfile(meta_tf[i, 2]), header=F, sep='\t', stringsAsFactors = T)
#        peaks_tf = peaks_tf[, 1:3]
#        colnames(peaks_tf) = c('chrom', 'start', 'end')
#        tryCatch({
#            # the number of GO terms in the genesets matters
#            results = chipenrich(peaks = peaks_tf, genome = 'hg19', genesets = gs_path,
#                          locusdef = '10kb', num_peak_threshold=1, qc_plots = FALSE, out_name = NULL, n_cores=1)
#            results.be = results$results
#            pval[[tf]] = results.be[, c('Geneset.ID', 'P.value', 'Odds.Ratio', 'Status', 'N.Geneset.Genes', 'N.Geneset.Peak.Genes', 'Geneset.Avg.Gene.Length')]
#        }, error=function(e) {conditionMessage(e)})
#    }
#    pval = do.call(rbind.data.frame, pval) 
#    ind = meta_tf$TF %in% pval[, 1]
#    ngene = sapply(hims_gs[['comb']][meta_tf$TF], nrow)
#    print(table(ind))
#    print(fivenum(ngene[ind]))
#    print(fivenum(ngene[!ind]))
#    return(pval)
#}
#
#pval = chip_enrichment(cell='gm12878')
#table(pval[, 'Status'], pval[, 'P.value']<=0.05)
#pval = chip_enrichment(cell='k562')
#table(pval[, 'Status'], pval[, 'P.value']<=0.05)
#pval = chip_enrichment(cell='hela')
#table(pval[, 'Status'], pval[, 'P.value']<=0.05)
## errors for some TFs & characterize the list
#    # could be 0 peaks
## refine the peaks by motif enrichment 
## window size 
## m
#RUNX3 NFATC1 ZBTB33  ESRRA   ELF1 
#20     53    193    334    640 
#PAX8 NR2F1  NFIC  IRF5  NFE2 
#14.0  82.5 161.0 345.0 869.0 # k562
#ZBED1  E4F1  NFIC   MYC  NFE2 
#42   137   249   467  1638 
#E2F8 HMBOX1  NR2F1   ARNT    SP1 
#80    225    326    516   1911 
## hela
#CTCF   CTCF  GABPA   REST NFE2L2 
#276    301    403    822   1164 
#EP300 EP300 EP300  MAFK  MAFK 
#230   230   513   796   796 