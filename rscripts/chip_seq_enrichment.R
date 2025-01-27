## doulbe check the gene set sampling  
## remove the whole region from the candidate region
## count the absolute number instead of summing the binary variable
## are there duplicated random samples? The answer is no.
rm(list=ls())
library(ggplot2)
source('src/chip_seq_enrichment_test.R')
####################################################
####### proportion report on HIMs with at least two TF having ChIP-seq data
####### proportion report on single TF is not strong
####################################################
cells = c('gm12878', 'k562')
ws = c('10000')
res = c()
for(cell in cells){
    for(w in ws){
        peak_gene = load_peak_table(cell=cell, w=w, filter=T, chip_coverage=0.15)
        peak_gene = peak_gene >=1
        tfs = colnames(peak_gene)
        gs_tf_chr = extract_him_genes_per_TF(i='data/him_summary_allinone.txt', cell=cell, tfs=tfs)
        gs_tf = gs_tf_chr[['genes']]
        tmp = compute_prop_gene_w_peaks(l=gs_tf, df=peak_gene)
        hims_dup = tmp[duplicated(tmp[, 'him']), 'him']
        tmp = tmp[tmp[, 'him'] %in% hims_dup, ]
        ## only keep HIMs with at least two TFs
        tmp = data.frame(cell=cell, w=w, tmp, stringsAsFactors = F)
        res = rbind(res, tmp)
    }
}
res = res[order(res[, 'cell'], res[, 'him']), ]
res$w = factor(res$w, levels=ws)
aggregate(prop~cell, res, FUN=function(z) sum(z>=50)/length(z))
round(sum(res$prop >= 50) / nrow(res) * 100 , 2)

####################################################
####### proportion report on single TF
####################################################
cells = c('gm12878', 'k562')
#ws = c('5000', '10000', '50000', '100000')
ws = c('10000')
# steop 0. report the number of TFs overall
tf_l = lapply(cells, function(z) colnames(load_peak_table(cell=z, w='5000', filter=T, chip_coverage=0.15)))
sapply(tf_l, length)
length(unique(unlist(tf_l)))
# step 1. compute the proportions 
res = c()
for(cell in cells){
    for(w in ws){
        peak_gene = load_peak_table(cell=cell, w=w, filter=T, chip_coverage=0.15)
        peak_gene = peak_gene >=1
        tfs = colnames(peak_gene)
        gs_tf_chr = extract_him_genes_per_TF(i='data/him_summary_allinone.txt', cell=cell, tfs=tfs)
        gs_tf = gs_tf_chr[['genes']]
        tmp = compute_prop_gene_w_peaks(l=gs_tf, df=peak_gene)
        tmp = data.frame(cell=cell, w=w, tmp, stringsAsFactors = F)
        res = rbind(res, tmp)
    }
}
head(res)
res$w = factor(res$w, levels=ws)
aggregate(prop~cell, res, FUN=function(z) sum(z>=50)/length(z))
round(sum(res$prop >= 50) / nrow(res) * 100 , 2)
fivenum(res$ngene)

cellnew  = c('gm12878' = 'GM12878', 'k562' = "K562")
nhim = sapply(cells, function(z) length(unique(res[res$cell == z, 'him'])))
cellnew = paste(cellnew,'\n(# HIMs=',  nhim, ')', sep='')
# results is not interesting at all.
pdf('sup_fig/chip_seq_peak.pdf', width=2, height=3)
ggplot(res, aes(x=cell, y=prop)) + geom_boxplot(width=0.2) + 
    ylab('% of genes in a HIM have TF ChIP-seq peaks') + 
    scale_x_discrete(labels=cellnew) + 
    theme_classic() + 
    theme(axis.title.y = element_text(size=8), axis.text=element_text(size=7)) + 
    theme(axis.title.x = element_blank())  
dev.off()
##############################
#### stops here 
##############################
stop()
table(res$cell)
# number of unique HIMs
# report this number in the figure
##
tmp1 = aggregate(prop~tf+cell, data=res, FUN=fivenum)
tmp2 = round(apply(peak_gene, 2, function(z) sum(z)/length(z)), 3)
tmp3 = aggregate(prop~tf, data=res, FUN=length)
tmp = cbind(tmp1, tmp2, tmp3)
tmp; nrow(tmp)
fivenum(res$prop)
stop()

####################################################
####### enrichment test
####################################################
cells = c('gm12878', 'k562')
ws = c('5000', '10000', '50000', '100000')
N=1000 # the number of random gene set per him
##### step 1. generate random sets
#library(doParallel); registerDoParallel(cores=2)
#foreach(cell=cells) %dopar%  random_geneset_output(cell=cell, N=N)
#q(save='no')
##### step 2. load random sets, compute pval, and save the output
res1 = c()
res2 = c()
for(cell in cells){
    for(w in ws){
        for(ty in c(T, F)){
            tmp = random_geneset_comppute_pval(cell=cell, w=w, binary=T)
            tmp1 = tmp[['individual']]
            tmp2 = tmp[['group']]
            tmp1 = data.frame(cell, w, ty, tmp1, stringsAsFactors = F)
            tmp2 = data.frame(cell, w, ty, tmp2, stringsAsFactors = F)
            res1 = rbind(res1, tmp1)
            res2 = rbind(res2, tmp2)
        }
    }
}
write.table(res1, 'inter_results/chip_enrichment_individual.txt', row.names = F, sep='\t', quote=F)
write.table(res2, 'inter_results/chip_enrichment_group.txt', row.names = F, sep='\t', quote=F)
q(save='no')

######## method overhual again
source('src/chip_seq_enrichment_test.R')
#cells = c('gm12878', 'k562')
cells = c('k562')
# window size does not matter much
#ws = c('5000', '10000', '50000', '100000')
#ws = c('5000', '10000')
ws = c('100000')
N=1000 # the number of random gene set per him
for(cell in cells){
    for(w in ws){
        peak_gene = load_peak_table(cell=cell, w=w, filter=F, chip_coverage=0.05)
        # convert to binary
        peak_gene = peak_gene >= 1
        tfs = colnames(peak_gene)
        # only consider master TFs with chip-seq data
        #tfs.master = read.table(paste('data/master_tfs_subset_grn_',cell,'.txt', sep='') , stringsAsFactors = F)[,1]
        #tfs = tfs[tfs%in% tfs.master]
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
            tmp = random_geneset_pval(i=ncomb[k, 1], j=ncomb[k, 2], gs_tf=gs_tf, gs_chr=gs_chr, 
                                    peak_gene=peak_gene, gene_dist=gene_dist, N=N, parallel=F)  
            if(!is.null(tmp)){
              res = rbind(res, tmp)
            }
        }
        res = data.frame(res, stringsAsFactors = F)
        for(i in 3:6){
            res[, i] = as.numeric(res[, i])
        }
        res[, 'prop'] = res[, 4] / res[, 3] * 100
        ofile = paste('inter_results/chip_enrichment_', cell, '_w_', w, '.txt',sep='')
        write.table(res, ofile, col.names = F, row.names = F, quote=F, sep='\t')
    }
}
# summarize results
res = summarize_enrich(cell='gm12878', w='5000')
res = summarize_enrich(cell='gm12878', w='10000')
res = summarize_enrich(cell='gm12878', w='50000')
res = summarize_enrich(cell='gm12878', w='100000')
##
res = summarize_enrich(cell='k562', w='5000')
res = summarize_enrich(cell='gm12878', w='10000')
#### method overhaul
result = c()
for(cell in cells){
    hims_gs = gs_from_hims('data/him_summary_allinone.txt', cell=cell)
    for(w in ws){
        for(ty in type){
            for(cc in chip_coverage_cutoffs){
                tmp = chipenrich_him_main(him_gs=him_gs, cell=cell, w=w, method.test=ty, tf_group=F, chip_coverage = cc)
                tmp = data.frame(cell, w, ty,cc, tmp, stringsAsFactors = F)
                result = rbind(result, tmp)
            }
        }
    }
}
result$w = factor(result$w, levels=ws)
result$cc = factor(result$cc, levels=chip_coverage_cutoffs)
result_pop = subset(result, ty=='population')
result_pop_long = melt(result_pop, measure.vars = c('prop', 'pval.1', 'background'))
ggplot(result_pop_long, aes(y=Prop.background, x=as.factor(value), color=variable)) + geom_violin() +
    facet_grid(cell+cc~w)
result_pop_sum = aggregate(cbind(prop, pval.1, background)~cell+ w + cc , data=result_pop,
                        FUN=function(z) sum(z) / length(z))
result_pop_sum = result_pop_sum[order(result_pop_sum$cell, result_pop_sum$w),]
result_pop_sum
# figure for paper
ggplot(subset(result_pop, cc==0.1), aes(x=Prop, y= -1*log10(pval), color=pval<=0.05)) + geom_point(size=0.5) + 
    geom_vline(xintercept = 0.8, linetype='dotted') +  facet_grid(cell~w) 
# ignore pvalue 
library(ggrepel)
ggplot(subset(result_pop, cc==0.1), aes(x=1, y=Prop, label=tf)) + geom_point(size=0.5) + 
    geom_text_repel() + facet_grid(cell~w) 

# aggregate TFs in invidual type, x-axis is prop.background
result_indi = aggregate(cbind(prop, pval.1, background) ~ cell+ w + tf + Prop.background, data=result,
                        FUN=function(z) sum(z) / length(z),
                        subset=result$ty == 'individual')
result_indi_long = melt(result_indi, measure.vars=c('prop', 'pval.1', 'background'))
ggplot(result_indi_long, aes(x=Prop.background, y=value, color=variable)) + geom_point() +
    facet_grid(cell~w)
result_indi_sum = aggregate(cbind(prop, pval.1, background) ~ cell+ w, data=result_indi,
                        FUN=function(z) sum(z>=0.5) / length(z))
result_indi_sum = result_indi_sum[order(result_indi_sum$cell, result_indi_sum$w), ]
result_indi_sum
# do not aggregate; plot him_ids and TFs
result_indi_him = subset(result, ty=='individual')
ggplot(result_indi_him, aes(x=as.factor(tf), y=Prop)) + geom_point(size=0.5)+
    facet_grid(cell~w)

# observation 1. window size does not matter much for pval & backgroud, but matters a lot on absolute proportion.
               
# thought, reconsider the definition of enrichment
# pval is not suitable for individual hims due to small gene number in hims; this is done
# try different statistical test; hypergeo, binomial, fisher exact; not necessary
# number of targets per him; it could be a confounding variable 
# output the number for each TF to test whether there are significant variation among TFs, and the factors contribute to the variations 
# change the value of majority (0.5)