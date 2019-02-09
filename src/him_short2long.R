### input is a dataframe
## cn must have at least 3 columns: him genes; cell type; him id 
library(splitstackshape)
short2long = function(X, CN, cn='genes'){
    hims = X[, CN]
    hims = data.frame(him_id = rownames(hims), hims, stringsAsFactors = F)
    hims = cSplit(hims, cn, sep=';', direction='long')
    hims = cSplit(hims, cn, sep=',', direction='long')
    hims = as.data.frame(hims)
    hims$genes = as.character(hims[, cn])
    return(hims)
}
