# this function takes a cell annotation file, a peak annotation file, and a peak count matrix, and then return a 
# cds object with monole2
ATAC_cds_construct <- function (UMI, df_cell, df_gene) 
{
    
    pd = new("AnnotatedDataFrame", data = df_cell)
    fd = new("AnnotatedDataFrame", data = df_gene)
    colnames(UMI) = df_cell$sample
    row.names(UMI) = df_gene$peak
    row.names(pd) = colnames(UMI)
    row.names(fd) = row.names(UMI)
    cds = newCellDataSet(UMI, phenoData = pd, featureData = fd, 
        expressionFamily = negbinomial.size())
    return(cds)
}