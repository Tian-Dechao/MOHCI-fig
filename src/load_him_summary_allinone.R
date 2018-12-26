X = read.table('data/him_summary_allinone.txt', header=T, stringsAsFactors = F, sep='\t')
cn = colnames(X)
# load and combine tsa
tsa = read.table('data/k562_tsa.txt', header=T, row.names=NULL, sep='\t', stringsAsFactors = F)
rownames(tsa) = paste('k562_', tsa$index, sep='')
tsa = tsa[, !colnames(tsa) %in% c('index', 'source')]
X = data.frame(X, tsa[rownames(X),], stringsAsFactors = F)
# Inf with Max
for(i in 1:ncol(X)){
    X[is.infinite(X[, i]), i] = max(X[is.finite(X[, i]), i], na.rm = T)
}
# compute % genes of a HIM in A compartments; 
#X[, 'A_percent'] = X[, 'A_num'] /rowSums(X[, c('A_num', 'B_num')]) * 100
# convert CV to CV * 100
#cncv = cn[grepl('_cv$', cn)]
#for(i in cncv){
#    X[, i] = X[, i] * 100
#}
# create a new feature to compare all,tf,gene vs none.j
#X[, 'conserve_staVScs'] = 'stable'
#jX[X[, 'conserve_All'] == "none", 'conserve_staVScs'] = 'cs'
#X[X[, 'source'] != "him", 'conserve_staVScs'] = NA
sum(table(X[, 'conserve_All']))
sum(table(X[, 'conserve_staVScs']))
table(X$conserve_staVScs, X$source)
# master proportion 
#X[, 'master_in_p'] = X[, 'master_in'] / X[, 'TF_number'] * 100
# ChIA-PET per MB
#X[, 'ChIA_num_p'] = X[, 'ChIA_num'] / X[, 'gene_distance'] * 100
# group the features
cells = sort(unique(X$cell))
fchara = c('edge_density', 'tfgene_density', 'motif_density', 'A_percent', 'repli_mean', 'repli_cv', 'ppi_density')
fsnp = c('cosmic_snp_frequency', 'snp_single_frequency', 'snp_ins_frequency', 'snp_del_frequency')
fessg = c('ess_all_gene_p','ess_con_gene_p', 'ess_k562_gene_p', 'ess_k562_specific_gene_p')
fesst = c('ess_all_tf_p', 'ess_con_tf_p', 'ess_k562_tf_p', 'ess_k562_specific_tf_p')
fcsg = c('csgp', 'cstp')
fexpg = c('gene_expression_mean', 'gene_expression_cv')
fmas = c('master_in', 'master_in_p', 'TF_number', 'gene_num', 'gene_distance')
fhmm = cn[grepl('ratio', cn) & grepl('X', cn)]
featratio = cn[grepl('ratio', cn) & !grepl('X', cn)]
fchia = c('ChIA_num', 'ChIA_num_p')
fse = c('SE_num_100kb', 'SE_num_100kb_per_Mb', 'SE_density_100kb','SE_num_500kb', 'SE_num_500kb_per_Mb', 'SE_density_500kb', 'SE_num_1Mb', 'SE_num_1Mb_per_Mb','SE_density_1Mb')
fsef = cn[grepl('SE', cn)]
