rm(list=ls())
library(splitstackshape)
#library(chipenrich)
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
load_peak_table = function(cell, w){
    path = paste('data/chip-peaks/gene_TF_peak_combined_',cell,'_w_', w, '.txt', sep='')
    peak = read.table(path, header = T, row.names = 1, sep='\t')
    return(peak)
}
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
        rownames(res) = him_u
    } else if(method.test=='population'){
        genes = df1[, 'genes']
        n = length(genes)
        gene_bool = rownames(df2) %in% genes
        n1 = sum(df2[gene_bool, ])
        p_him = n1 / n
        pval = binom.test(x=n1, n=n, p=p_back, alternative = 'g')$p.value
        res = c(n, n1, pval, p_him, p_back)
        res = t(as.data.frame(res))
        rownames(res) = 'Population'
    }
    colnames(res) = c('N.gene', 'N.gene.w.peak', 'pval', 'Prop', 'Prop.background')
    #res = data.frame(him=him_u, res, stringsAsFactors = F)
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
chipenrich_result_summarize = function(df, method, method.test='individual'){
    if(method=='prop'){
        ind = df$Prop >= 0.8
    } else if (method=='pval'){
        ind = df$pval <= 0.05
    } else if (method=='background'){
        ind = df$Prop > df$Prop.background 
    } else {
        stop("Double chcek the method")
    }
    if(method.test == 'individual'){
        tmp0 = as.data.frame(table(df[, 'tf'])); tmp0[, 1] = as.character(tmp0[, 1])
        tmp1 = as.data.frame(table(df[ind, 'tf'])); tmp1[, 1] = as.character(tmp1[, 1])
        rownames(tmp0) = tmp0[, 1]; rownames(tmp1) = tmp1[, 1]
        tmp2 = cbind(tmp0[rownames(tmp1), ], tmp1)
        ind = tmp2[, 4] / tmp2[, 2] >= 0.5
    } 
    
    ind = as.integer(ind)
    tmp2.sum = as.data.frame(table(ind))
    tmp2.sum[,1] = as.character(tmp2.sum[, 1])
    tmp3 = c(0, 0); names(tmp3) = as.character(0:1)
    tmp3[tmp2.sum[, 1]] = tmp2.sum[, 2]
    return(tmp3)
}
chipenrich_him_main = function(cell, w, method.test){
    hims_gs = gs_from_hims('data/him_summary_allinone.txt', cell=cell)
    peak_gene = load_peak_table(cell=cell, w=w)
    tfs_wc = colnames(peak_gene)
    res = chipenrich_him_combine(tfs=tfs_wc, hims_gs=hims_gs, peak_gene = peak_gene, method.test=method.test)
    sum.method = c('prop', 'pval', 'background')
    res.sum = sapply(sum.method, function(z) chipenrich_result_summarize(df=res, method=z, method.test=method.test) )
    return(res.sum)
}
#cell='gm12878';w='10000'
#chipenrich_him_main(cell=cell, w=w, method.test='individual')
#chipenrich_him_main(cell=cell, w=w, method.test='population')
cells = c('gm12878', 'hela', 'k562')
ws = c('5000', '10000', '50000', '100000')
result = list()
for(cell in cells){
    result[[cell]] = list()
    for(w in ws){
        print(c(cell, w))
        result[[cell]][[w]] = chipenrich_him_main(cell=cell, w=w, method.test='population')
    }
}
result
result_individual = list()
for(cell in cells){
    result_individual[[cell]] = list()
    for(w in ws){
        print(c(cell, w))
        result_individual[[cell]][[w]] = chipenrich_him_main(cell=cell, w=w, method.test='population')
    }
}
result_individual
# observation 1. window size does not matter much for pval & backgroud, but matters a lot on absolute proportion.
# thought, reconsider the definition of enrichment
# try different statistical test; hypergeo, binomial, fisher exact
# number of targets per him; it could be a confounding variable 
#### method  overhaul stops here 

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

TFs_w_chip_seq = function(cell){
    meta = read.table(paste('data/', cell, '/metadata_slim_from_all.tsv', sep=''), header=F, stringsAsFactors = F, sep='\t')
    meta = meta[, c(1, 13)]
    meta[, 'filepath'] = paste('data/', cell, '/', meta[,2], '_', meta[,1], '.bed.gz', sep='')
    meta = meta[, c(2, 3)]
    colnames(meta) = c('TF', 'filepath')
    meta$TF = gsub('-human','', meta$TF)
    rownames(meta) = meta$TF
    return(meta)
}
# copy bed files from the server in case they do not exist in the local directory
#scp -r cmu:/hive/dechaot/Data/ENCODE/chip-seq/gm12878 data/
chip_enrichment = function(cell){
    meta_tf = TFs_w_chip_seq(cell=cell)
    hims_gs = gs_from_hims('data/him_summary_allinone.txt', cell=cell)
    l1 = names(hims_gs[['comb']]); l2 = meta_tf$TF
    # some TFs with Chip-seq data are not in the GRNs such as JUND
    meta_tf = meta_tf[l2 %in% l1, ]
    meta_tf['NFE2', ]
    pval = list()
    i = which(meta_tf$TF == 'NFE2')
    for(i in 1:nrow(meta_tf)){
        tf = meta_tf[i, 1]
        gs_path = paste('inter_results/gs_', tf,'.txt', sep='')
        write.table(hims_gs[['comb']][[tf]], file=gs_path, col.names = T, row.names = F, sep='\t', quote=F)
        peaks_tf = read.table(gzfile(meta_tf[i, 2]), header=F, sep='\t', stringsAsFactors = T)
        peaks_tf = peaks_tf[, 1:3]
        colnames(peaks_tf) = c('chrom', 'start', 'end')
        tryCatch({
            # the number of GO terms in the genesets matters
            results = chipenrich(peaks = peaks_tf, genome = 'hg19', genesets = gs_path,
                          locusdef = '10kb', num_peak_threshold=1, qc_plots = FALSE, out_name = NULL, n_cores=1)
            results.be = results$results
            pval[[tf]] = results.be[, c('Geneset.ID', 'P.value', 'Odds.Ratio', 'Status', 'N.Geneset.Genes', 'N.Geneset.Peak.Genes', 'Geneset.Avg.Gene.Length')]
        }, error=function(e) {conditionMessage(e)})
    }
    pval = do.call(rbind.data.frame, pval) 
    ind = meta_tf$TF %in% pval[, 1]
    ngene = sapply(hims_gs[['comb']][meta_tf$TF], nrow)
    print(table(ind))
    print(fivenum(ngene[ind]))
    print(fivenum(ngene[!ind]))
    return(pval)
}

pval = chip_enrichment(cell='gm12878')
table(pval[, 'Status'], pval[, 'P.value']<=0.05)
pval = chip_enrichment(cell='k562')
table(pval[, 'Status'], pval[, 'P.value']<=0.05)
pval = chip_enrichment(cell='hela')
table(pval[, 'Status'], pval[, 'P.value']<=0.05)
# errors for some TFs & characterize the list
    # could be 0 peaks
# refine the peaks by motif enrichment 
# window size 
# m
RUNX3 NFATC1 ZBTB33  ESRRA   ELF1 
20     53    193    334    640 
PAX8 NR2F1  NFIC  IRF5  NFE2 
14.0  82.5 161.0 345.0 869.0 # k562
ZBED1  E4F1  NFIC   MYC  NFE2 
42   137   249   467  1638 
E2F8 HMBOX1  NR2F1   ARNT    SP1 
80    225    326    516   1911 
# hela
CTCF   CTCF  GABPA   REST NFE2L2 
276    301    403    822   1164 
EP300 EP300 EP300  MAFK  MAFK 
230   230   513   796   796 