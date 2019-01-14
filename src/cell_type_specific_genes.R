### create cell tyep specific genes from expression data 
fc  = function(y, nc=1){
    y1 = y[, nc]
    y2 = y[, (1:ncol(y))[(1:ncol(y)) != nc]]
    #y2c = apply(y2, 1, max)
    y2c = apply(y2, 1, median)
    # mean is much worse than median
    #y2c = apply(y2, 1, mean)
    fchange = y1 / y2c
    # no consistent pattern
    #g = rownames(y)[fchange > 2 & y[,nc]>=log(1+1)]
    # 1 is worse than 0.5
    #g = rownames(y)[fchange > 2 & y[,nc]>=1]
    # 0.5 . 2 cell types are not significant
    #g = rownames(y)[fchange > 2 & y[,nc]>=0.5]
    g = rownames(y)[fchange > 2 & y[,nc]>=0.1]
    # 2 is much bettern than 1.5
    #g = rownames(y)[fchange > 2]
    # one of the best
    #g = rownames(y)[fchange > 1.5]
    #g = rownames(y)[fchange > 1 & y[,nc]>=0.1]
    return(g)
}

load_rna_cs = function(file, cs){
    files = sapply(cs, function(z) gsub('PH', z, file))
    L = list()
    for(i in 1:length(cells)){
        tmp = read.table(files[i], header = T, stringsAsFactors = F, sep = '\t')
        rownames(tmp) = tmp[, 'name']
        L[[i]] = tmp
    }
    genes = rownames(L[[1]])
    x = matrix(0,nrow=length(genes), ncol=length(cells))
    rownames(x) = genes
    colnames(x) = cells
    for(i in 1:length(cells)){
        tmp = L[[i]]
        x[unlist(tmp[, 'name']), i] = unlist(tmp[, 'expression'])
    }
    # log transformation
    x = log(1 + x)
    cg = list()
    for(i in 1:ncol(x)){
        g = fc(y=x, nc=i)
        cg[[i]] = g
    }
    names(cg) = cs
    return(cg)
}