rm(list=ls())
source('src/chip_seq_enrichment_test.R')
cells = c('gm12878', 'k562')
ws = c('5000', '10000', '50000', '100000')
type = c('individual', 'population')
chip_coverage_cutoffs=c(0.05, 0.1, 0.2)
cell='gm12878'; w='10000'
peak_gene = load_peak_table(cell=cell, w=w)
tfs = colnames(peak_gene)
genes_bed = load_gene_grn_bed(cell=cell)
gene_dist = compute_pairwise_distance(genes_bed)
gs_tf_chr = extract_him_genes_per_TF(i='data/him_summary_allinone.txt', cell=cell, tfs=tfs)
gs_tf = gs_tf_chr[['genes']]; gs_chr = gs_tf_chr[['chr']]
source('src/chip_seq_enrichment_test.R')
gene_dist[['chr1']][1:5, 1:5]
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