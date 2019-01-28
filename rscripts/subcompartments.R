rm(list=ls())
source('src/boxdot.R')
load_subcompartment  = function(){
    sub = read.table('data/subcompartment_inter_genes.txt', header=F, stringsAsFactors = F)
    sub = sub[ order(sub[, 14], -sub[, 15]), ]
    sub = sub[!duplicated(sub[, 14]), ]
    sub_vec = sub[, 4]
    names(sub_vec) = sub[, 14]
    return(sub_vec)
}
load_hims_from_allinone = function(cell='gm12878'){
    hims = read.table('data/him_summary_allinone.txt', header=T, sep='\t', stringsAsFactors = F)
    hims = hims[hims$cell==cell & hims$source =='him', ]
    hims[, 'conserve_staVScs'] = 'stable'
    hims[hims[, 'conserve_All'] == "none", 'conserve_staVScs'] = 'cs'
    hims[hims[, 'source'] != "him", 'conserve_staVScs'] = NA
    return(hims)
}
# compute the proportion of genes in each subcompartments
subcompartment_distribution = function(genes, gene2sub, cutoff=0.5, dominate = T){
    if(dominate){
        if(sum(genes %in% names(gene2sub)) >= 5){
            freq = table(gene2sub[genes])
            freq = freq / length(genes)
            freq_m = max(freq)
            if(is.infinite(freq_m)){
                print(freq)
                stop()
            }
            if(freq_m > cutoff){
                sub_dom = names(freq)[which(freq == freq_m)]
            }else {
                sub_dom = 'NA'
            }
        } else {
            sub_dom = 'NA'
        }
        return(sub_dom)
    } else {
        subs = sort(unique(gene2sub))
        freq = table(gene2sub[genes])
        freq = freq / length(genes)
        res = rep(0, length(subs))
        names(res) = subs
        res[names(freq)] = as.vector(freq)
        return(res)
    }
}

extract_sub_freq = function(x, cutoffs, labels){
        y=cut(x, breaks=cutoffs, labels=labels, include.lowest = T, right=F)
        z = as.data.frame(table(y))
        z[, 2] = z[, 2] / sum(z[, 2]) * 100
        w = z[, 2]
        names(w) = z[, 1]
    return(w)
}
## main starts here
sub = load_subcompartment()
# proportion of subcompartments
prop.table(table(sub))
hims = load_hims_from_allinone(cell='gm12878')
# only consider hims with all their genes in A compartments
hims = hims[hims$A_percent == 100,]
hims_gene = lapply(hims$genes, function(z) unlist(strsplit(z, ';|,')))
hims_tf = lapply(hims$TFs, function(z) unlist(strsplit(z, ';|,')))
hims_sub = sapply(hims_gene, function(z) subcompartment_distribution(genes=z, gene2sub=sub, dominate =F))
hims_sub = hims_sub * 100
hims_sub = t(hims_sub)
rownames(hims_sub) = hims$index
fivenum(rowSums(hims_sub))
#sum(rowSums(hims_sub) ==0)
# 27 HIMs with all/majority their genes are not in any subcompartments
# optional
#hims_sub = hims_sub[ rowSums(hims_sub) > 0, ]
# create categories 
cutoffs = c(0, 50, 80,   99.99999999999,101)
labels = c('[0, 50)', '[50, 80)', '[80, 100)', '100' )
hims_sub_freq = apply(hims_sub[, 1:2], 2, function(z) extract_sub_freq(x=z, cutoffs = cutoffs, labels = labels))
hims_sub_freq
sum(hims_sub_freq[2:4, ])
sum(hims_sub_freq[3:4, ])
# prepare the table for figures 
hims_sub_freq_long = c()
hims_sub_freq_long[1:3] = hims_sub_freq[4:2, 1]
hims_sub_freq_long[5:7] = hims_sub_freq[2:4, 2]
hims_sub_freq_long[4] = 100 - sum(hims_sub_freq_long, na.rm=T)
names1 = paste(labels[4:2], '\nA1', sep='')
names2 = paste(labels[2:4], '\nA2', sep='')
names3 = c(names1, 'Others', names2)
names3 = factor(names3, levels=names3)
hims_sub_freq_long = data.frame(cat=names3, freq=hims_sub_freq_long)
# choose 3 colors 
# side by side bar
col_list = c('stable' = "#2b83ba", 'cs' = "#d7191c")
pdf('sup_fig/Subcomp_percent_gm12878.pdf', width=3, height = 2.5/2)
ggplot(hims_sub_freq_long, aes(x=cat, y=freq)) + geom_col(width=0.5, fill=col_list[1], alpha=0.5) + 
    ylab('% of HIMs') + xlab('% of genes in subcompartments') + 
    theme(axis.title = element_text(size=8), axis.text=element_text(size=6), 
          axis.text.x=element_text(angle=25, hjust=1), axis.ticks.x=element_blank())
dev.off()
# what's the conclusions     
A/B compartments can be furthur divided into subcompartments~(\citep{raoxx}), which provide an oppourtunity to 
refine the spatial lcoations of HIMs.
#Among the GM12878 HIMs with all their genes in A compartments, 
69.92 % of them are dominantly (>= 80% of the genes in a HIM) either in A1 or A2 subcompartments. 

# stop here
(apply(hims_sub, 2, function(z) sum(z>=1)) / nrow(hims_sub))
# 0.3135593 0.0819209 0.0000000 0.0000000 0.0000000 0.0000000 
#0.184590690 0.064205457 0.001605136 0.000000000 0.000000000 0.003210273 
(apply(hims_sub[, 1:2], 2, function(z) sum(z>=0.9)) / nrow(hims_sub)) 
# 0.4322034 0.1073446 
(apply(hims_sub[, 1:2], 2, function(z) sum(z>=0.8)) / nrow(hims_sub)) 
# 0.5593220 0.1694915 
# 0.5521669 0.1894061 
(apply(hims_sub[, 1:2], 2, function(z) sum(z>=0.5)) / nrow(hims_sub)) 
# 0.6892655 0.2909605 
(apply(hims_sub, 2, function(z) sum(z>=0.6)) / nrow(hims_sub))
fivenum(rowSums(hims_sub[, 1:2]))
#hims_sub = data.frame(hims_sub, conserve_staVScs= hims$conserve_staVScs, stringsAsFactors = F)
apply(hims_sub[hims$conserve_staVScs == 'stable', ], 2, fivenum)
apply(hims_sub[hims$conserve_staVScs == 'cs', ], 2, fivenum)
table(hims[, 'conserve_staVScs'])
# stop here
hims_sub = sapply(hims_gene, function(z) subcompartment_distribution(genes=z, gene2sub=sub, cutoff = 0.95))
table(hims_sub)
tf_A1 =  unlist(hims_tf[hims_sub == 'A1'])
tf_A2 =  unlist(hims_tf[hims_sub == 'A2'])
tf_freq = compute_freq(x = tf_A1, y = tf_A2)
head(tf_freq)
ggplot(tf_freq, aes(X1, X2)) + geom_point()
tf_A1_dominately = rownames(tf_freq)[ tf_freq[, 2] > 2 * tf_freq[, 3]]
tf_A2_dominately = rownames(tf_freq)[ tf_freq[, 3] > 2 * tf_freq[, 2]]
write.table(tf_A2_dominately, file='TFs_dominately in hims in A2.txt', row.names = F, quote=F, col.names = F)
write.table(tf_A1_dominately, file='TFs_dominately in hims in A1.txt', row.names = F, quote=F, col.names = F)
tfs_all = unique(unlist(hims_tf))
write.table(tf_A1_dominately, file='TFs in gm12878.txt', row.names = F, quote=F, col.names = F)

# 
ggplot(hims, aes(gene_expression_mean, TF_number)) + geom_point()
# GRN
gm = read.table('grn_gm12878_txStart.txt', header=T, sep='\t', stringsAsFactors = F)
deg = as.data.frame(table(gm[, 1]))
deg[,1] = as.character(deg[, 1])
rownames(deg) = deg[, 1]
deg = deg[order(-deg[, 2]), ]
sum(deg[, 2] >= deg['CTCF', 2])
sum(deg[, 2] >= deg['YY1', 2])
# HIMs
hims = read.table('him_summary_allinone.txt', header=T, sep='\t', stringsAsFactors = F)
hims['gm12878_267', 'gene_expression_mean']
hims['gm12878_267', 'gene_expression_mean_ratio']
hims['gm12878_628', 'gene_expression_mean']
hims['gm12878_267', 'genes']
colnames(hims)
