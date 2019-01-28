rm(list=ls())

rna = read.table('data/expression/k562_gene_expression.txt', header=T, row.names = 1, sep='\t')
genes = unlist(strsplit('NKIRAS1;RPL15,NR1D2,SATB1,TBC1D5,THRB,UBE2E1', ';|,'))
genes
rna[genes, 'expression']
