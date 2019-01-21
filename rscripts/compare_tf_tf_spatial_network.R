rm(list=ls())
setwd('~/Documents/graph_cluster/')
# use this one
# load the TF-TF spatial interaction data 
motifs = function(k=3){
    if(k == 3){
       res = c('ETS1', 'GABPA', 'ELK1',
               'SP1', 'ETS1','ELF1',
               'SP1', 'GABPA','ELF1',
               'E2F4', 'ETS1','SP1',
               'E2F4', 'GABPA','SP1',
               #'SP1', 'MAZ','MAX',  # MAZ is not there 
               'ETS1', 'GABPA','E2F4',
               'ETS1', 'GABPA','ELF1',
               'ETS1', 'GABPA','SP1'
               ) 
       res[ !res %in% hims$TF]
       res = matrix(res, ncol=k, byrow=T)
    } else if (k == 4) {
        res = c('ETS1', 'GABPA', 'SP1', 'E2F4', 'ETS1', 'GABPA', 'SP1', 'ELF1')
        res = matrix(res, ncol=k, byrow=T)
    }
    
   return(res) 
}
load_TF_network = function(filter_list=c()){
    g1 = read.table('data_cano/gm12878_TF_TF.txt', header=T, stringsAsFactors = F, sep='\t')
    if(length(filter) > 0){
        ind1 = g1[, 1] %in% filter_list
        ind2 = g1[, 2] %in% filter_list
        g1 = g1[ ind1 & ind2, ]
    }
    rn1 = paste(g1[, 1], g1[, 2], sep='_')
    g1 = g1[, c(1, 2, 5)]
    rownames(g1) = apply(g1[, 1:2],1, function(z) paste(sort(z), collapse = '_'))
    ind1 = g1[, 3] == 1
    TF_g1 = sort(unique(c(g1[ind1, 1], g1[ind1, 2])))
    TF_g1 = data.frame(TF = TF_g1, group = 1, stringsAsFactors = F)
    TF_g2 = sort(unique(c(g1[!ind1, 1], g1[!ind1, 2])))
    TF_g2 = data.frame(TF = TF_g2, group = 2, stringsAsFactors = F)
    TF_g = rbind(TF_g1, TF_g2)
    rownames(TF_g) = TF_g[, 'TF']
    res = list(G=g1, TF = TF_g)
    return(res)
}
load_hims = function(cell='gm12878'){
    hims = read.table('him_summary_allinone.txt', header=T, sep='\t', stringsAsFactors = F)
    hims = hims[hims$cell==cell & hims$source =='him', ]
    TFs_list = lapply(hims$TFs, function(z) unique(unlist(strsplit(z, ';|,'))))
    names(TFs_list) = hims$index
    TFs = sort(unique(unlist(TFs_list)))
    res = list(hims=hims, TF_list = TFs_list, TF = TFs)
    return(res)
}
comp_freq = function(TFs, k=2, TF_list){
    TF_pair = combn(TFs, k); TF_pair = t(TF_pair)
    freq = c()
    for(i in 1:nrow(TF_pair)){
        freq[i] = sum(sapply(TF_list, function(z) sum(TF_pair[i, ] %in% z) == k))
    }
    res = data.frame(TF_pair, freq, stringsAsFactors = F)
    colnames(res) = c(paste('TF', 1:k, sep=''), 'freq')
    rownames(res) = apply(res[, 1:k], 1, function(z) paste(sort(z), collapse = '_'))
    return(res)
}
TF_in_hims = function(v, l, k=0){
    res = lapply(l, function(z) v[v %in% z])
    names(res) = names(l)
    if(k==0){
        res = res[sapply(res, length) > 0]
    } else {
        res = res[sapply(res, length) == k]
    }
    return(res)
}
hims = load_hims(cell = 'gm12878')
net = load_TF_network(filter_list=hims$TF)
# double check that all TFs in net are in hims
TF_in_num_hims = comp_freq(TFs = net$TF[net$TF[, 'group'] == 1, 1], k=2, TF_list = hims$TF_list)
TF_in_num_hims = data.frame(TF_in_num_hims, connected = rownames(TF_in_num_hims) %in% rownames(net$G))
aggregate(TF_in_num_hims$freq, by=list(TF_in_num_hims$connected), function(z) return(c(sum(z>0), length(z))))
# motfs 
motifs_3node = motifs(k=3)
apply(motifs_3node, 1, function(z) TF_in_hims(v=z, l=hims$TF_list, k=3))
motifs_4node = motifs(k=4)
apply(motifs_4node, 1, function(z) TF_in_hims(v=z, l=hims$TF_list, k=4))
# CTCF
TF_in_hims(v=c('CTCF'), l=hims$TF_list, k=1)
TF_in_hims(v=c('YY1'), l=hims$TF_list, k=1)
TF_freq = sapply(hims$TF, function(z) length(TF_in_hims(v=c(z), l=hims$TF_list, k=1)))
names(TF_freq) = hims$TF
fivenum(TF_freq)


## number of HIMs in CTCF and YY1
cells = c('gm12878', 'hela', 'huvec', 'k562', 'nhek')
res = c()
for(cell in cells){
    hims = load_hims(cell=cell) 
    TF_freq = sapply(hims$TF, function(z) length(TF_in_hims(v=c(z), l=hims$TF_list, k=1)))
    names(TF_freq) = hims$TF
    tmp = c(TF_freq[c('CTCF', 'YY1', 'ZNF143', 'AP1')], fivenum(TF_freq))
    res = rbind(res, tmp)
    print(cell)
    print(sum(TF_freq['CTCF'] > TF_freq) / length(TF_freq) * 100)
    print(sum(TF_freq['YY1'] > TF_freq) / length(TF_freq) * 100)
    print(sum(TF_freq['ZNF143'] > TF_freq) / length(TF_freq) * 100)
    print(tmp)
}
rownames(res) = cells
res
# compare out deg vs number of hims
hims = load_hims(cell='k562') 
TF_freq = sapply(hims$TF, function(z) length(TF_in_hims(v=c(z), l=hims$TF_list, k=1)))
names(TF_freq) = hims$TF
g = load_grn('grn_k562_txStart.txt')
deg = node_deg(el=g)
res = data.frame(nhim=TF_freq, deg=deg[names(TF_freq)], stringsAsFactors = F)
library(ggplot2)
ggplot(res, aes(y=nhim, x=deg.Freq)) + geom_point() + geom_smooth()
TF_inter = c('CTCF', 'YY1', 'ZNF143')
res[TF_inter,]
