rm(list=ls())
library(splitstackshape)
library(chipenrich)
# it is a model-based method to test enrichment of TF Chip-seq peaks in a given gene sets
# prepare gene set files from genes assigned to HIMs
gs_from_hims = function(i, cell){
    i = 'data/him_summary_allinone.txt'
    cell='gm12878'
    hims = read.table(i, header=T, sep='\t', stringsAsFactors = F)
    hims = hims[hims$source == 'him', ]
    hims = hims[grepl(cell, rownames(hims)), ]
    hims = hims[, c('index', 'TFs', 'genes')]
    hims = data.frame(gs_id = rownames(hims), hims, stringsAsFactors = F)
    hims = cSplit(hims, 'genes', sep=';', direction='long')
    hims = cSplit(hims, 'genes', sep=',', direction='long')
    hims$genes = as.character(hims$genes)
    # gene name to Entrez ID
    library("org.Hs.eg.db")
    ENTREZ <- as.list(org.Hs.egALIAS2EG); ENTREZ <- ENTREZ[!is.na(ENTREZ)]
    sym2entre = ENTREZ[hims$genes]; sym2entre = sapply(sym2entre, function(z) paste(z, collapse = ';'))
    hims$entrez = sym2entre
    # some gene symbols has multiple ENTREZ gene id. As a results, 8k to 10k.
    # This should not have big impact on the final results
    # refine this part if the results are not appelling.
    hims = cSplit(hims, 'entrez', ';', direction='long')
    # find the TF list
    TFs = unique(unlist(sapply(hims$TFs, function(z) unlist(strsplit(z, ";|,")))))
    # modify this part
    result_comb = list(); result_individual = list()
    for(i in TFs){
        ind = grepl(i, hims$TFs)
        genes = hims[ind, 'entrez']
        result_comb[[i]] = data.frame(gs_id=i, entrez=genes, stringsAsFactors = F)
        result_individual[[i]] = hims[ind, c('gs_id', 'entrez')]
    }
    result = list(comb=result_comb, individual=result_individual)
    return(result)
}

## modify the input Chip-seq peaks
#peaks_xx = read.table('xx')
data(peaks_E2F4, package = 'chipenrich.data')

hims_gs = gs_from_hims('data/him_summary_allinone.txt', cell='gm12878')
gs_path = 'inter_results/gs_E2F4.txt'
write.table(hims_gs[['comb']][['E2F4']], file=gs_path, col.names = T, row.names = F, sep='\t', quote=F)
results = chipenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path,
                      locusdef = '10kb', qc_plots = FALSE, out_name = NULL, n_cores=1)
results.be = results$results
print(results.be[, c('Geneset.ID', 'P.value', 'Odds.Ratio', 'Status', 'N.Geneset.Genes', 'N.Geneset.Peak.Genes', 'Geneset.Avg.Gene.Length')])
#results = polyenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path,
#                      locusdef = '10kb', qc_plots = FALSE, out_name = NULL, n_cores=1)
#results.be = results$results
#print(results.be[, c('Geneset.ID', 'P.value', 'Odds.Ratio', 'Status', 'N.Geneset.Genes', 'N.Geneset.Peak.Genes', 'Geneset.Avg.Gene.Length')])
## use the hybrid test but there are more cases that this test  encounter bugs 
#results = hybridenrich(peaks = peaks_E2F4, genome = 'hg19', genesets = gs_path,
#                      locusdef = '10kb', qc_plots = FALSE, out_name = NULL, n_cores=1)
#results.be = results$results
#print(results.be[, c('Geneset.ID', 'Odds.Ratio', 'N.Geneset.Genes', 'N.Geneset.Peak.Genes', 'Geneset.Avg.Gene.Length', 
                     'P.value.x', 'Status.x', 'P.value.y', 'Status.y', 'Status.Hybrid')])
