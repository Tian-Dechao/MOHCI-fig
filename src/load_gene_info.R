G = read.table('data/gene_info_allinone.txt', sep='\t', header=T)
cells = unique(G$cell)
indHIM = G$HIM == 1
G[indHIM, 'HIM'] = 'Yes'
G[!indHIM, 'HIM'] = 'No'
G$HIM = factor(G$HIM)
table(G$HIM)
# log transform gene expression 
G$expression = log2(1+ G$expression)
# compute A percent 
table(G$HIM, G$compartment, G$cell)
colnames(G)
# A genes are white balls; B genes are black balls
vennA = c()
for(i in 1:length(cells)){
  indA = G$compartment ==  'A' & G$cell == cells[i]
  indH = G$HIM == 'Yes' & G$cell == cells[i]
  nwdraw = sum(indA & indH)
  nw = sum(indA)
  nb = sum(!indA)
  ndraw = sum(indH)
  nwdraw;nw;nb;ndraw 
  logp = phyper(q=nwdraw, m=nw, n=nb, k=ndraw, lower.tail = F, log.p=T)  
  result = c(ndraw - nwdraw, nwdraw, nw - nwdraw, logp)
  vennA = rbind(vennA, result)
}
rownames(vennA) = cells
colnames(vennA) = c('HIM_B', 'HIM_A', 'nonHIM_A', 'Pval')
vennA
# essential genes are white balls; B genes are black balls
ess_comb = read.table('data/essential_combined.txt', stringsAsFactors = F)[, 1]
ess_k562 = read.table('data/essential_genes_k562.txt', stringsAsFactors = F)[, 1]
G[, 'ess_comb'] = G$gene_name %in% ess_comb
indE = G$ess_comb
G[indE, 'ess_comb'] = 'Yes'
G[!indE, 'ess_comb'] = 'No'
G[, 'ess_k562'] = G$gene_name %in% ess_k562
G[G$cell != 'k562', 'ess_k562'] = NA
indE = G$ess_k562 & !is.na(G$ess_k562)
indnE = !G$ess_k562 & !is.na(G$ess_k562)
G[indE, 'ess_k562'] = 'Yes'
G[indnE, 'ess_k562'] = 'No'
table(G$ess_k562)
vennE = c()
for(i in 1:length(cells)){
  indE = G$ess_comb ==  'Yes' & G$cell == cells[i]
  indH = G$HIM == 'Yes' & G$cell == cells[i]
  nwdraw = sum(indE & indH)
  nw = sum(indE)
  nb = sum(!indE)
  ndraw = sum(indH)
  nwdraw;nw;nb;ndraw 
  logp = phyper(q=nwdraw, m=nw, n=nb, k=ndraw, lower.tail = F, log.p=T)  
  result = c(ndraw - nwdraw, nwdraw, nw - nwdraw, logp)
  vennE = rbind(vennE, result)
}
rownames(vennE) = cells
colnames(vennE) = c('HIM_ne', 'HIM_e', 'nonHIM_e', 'Pval')
vennE
propE = matrix(0, nrow=5, ncol=2)
propE.p = c()
for(i in 1:length(cells)){
    Gsub = G[G$cell == cells[i], ]
    indE = Gsub$ess_comb == 'Yes'
    indH = Gsub$HIM == 'Yes'
    n1  = sum(indH); n2=sum(indH & indE); p1 = round(n2 / n1 * 100, 3)
    n3 = sum(!indH); n4=sum(!indH & indE); p2 = round(n4/n3 * 100, 3)
    propE[i, ] = c(p1, p2)
    pi = chisq.test(matrix(c(n1-n2, n2, n3-n4, n4), byrow = T, ncol=2))$p.value
    propE.p[i] = pi
}
propE.p = format(propE.p, scientific = T, digits=3)
propE.p = paste('P=',propE.p, sep='')
propE.p.df = data.frame(cell=cells, pval=propE.p, x=1.5, y=24)
colnames(propE) = c('Yes', 'No')
propE.df = data.frame(cell=cells, propE)
propE.long = melt(data=propE.df, id.vars='cell', measure.vars=c('Yes', 'No'), variable.name='HIM', value.name='prop')
# report some numbers
(propE.df$Yes / propE.df$No -1) * 100 
propEk = c()
Gsub = G[G$cell == 'k562', ]
indE = Gsub$ess_k562 == 'Yes'
indH = Gsub$HIM == 'Yes'
n1 = sum(indH); n2=sum(indH &indE); p1=round(n2/n1 * 100, 3)
n3 = sum(!indH); n4=sum(!indH &indE); p2=round(n4/n3 * 100, 3)
propEk = data.frame(essp=c(p1, p2), HIM=c('Yes', 'No'))
propEk$HIM = factor(propEk$HIM, levels=c('Yes', 'No'), ordered=T)
Ekmat = matrix(c(n1-n2, n2, n3-n4, n4), byrow=T, nrow=2)
rownames(Ekmat) = c('genes in HIMs', 'genes not in HIMs')
colnames(Ekmat) = c('Non-essential', 'Essential')
Ekmat
Ek.p = chisq.test(Ekmat)$p.value
Ek.p = format(Ek.p, scientific = T, digits=3)
Ek.p = paste('P=', Ek.p, sep='')
Ek.p.df = data.frame(x=1.5, y=14, pval=Ek.p, stringsAsFactors = F)

# compute P values 


# compute P value for features 
# compute P value for repli
pval = c()
alter = c()
for( i in 1:length(cells)){
    Gsub = G[G$cell == cells[i], ]
    x = Gsub[Gsub$HIM=="Yes", 'repli']
    y = Gsub[Gsub$HIM=="No", 'repli']
    ptmp = compare2vect(x, y)
    pval[i] = ptmp$pval
    alter[i] = ptmp$alter
}
miny = min(G$repli, na.rm=T)
maxy = max(G$repli, na.rm=T)
maxy = maxy + 0.05 * (maxy - miny)
Gpval_table = data.frame(cell=cells, pval=pval, alter=alter, maxy=maxy )
Gpval_table$pval = format(Gpval_table$pval, scientific=T, digits=3)
Gpval_table$pval = paste('P=', Gpval_table$pval, sep='')
indna = grepl('NA', Gpval_table[, 'pval'])
Gpval_table[indna, 'pval'] = ''
