
# combine cds for ATAC-seq data
ATAC_cds_combine <- function (cds_1, cds_2) 
{
    df_cell_1 = pData(cds_1)
    df_gene_1 = fData(cds_1)
    gene_count_1 = exprs(cds_1)
    df_cell_2 = pData(cds_2)
    df_gene_2 = fData(cds_2)
    gene_count_2 = exprs(cds_2)
    df_gene = df_gene_1 %>% filter(peak %in% df_gene_2$peak)
    cat("\n Number of unmatched peak: ", sum(df_gene_1$peak != df_gene_2$peak))
    cat("\nnumber of gene in cds 1: ", nrow(df_gene_1))
    cat("\nnumber of gene in cds 2: ", nrow(df_gene_2))
    cat("\nnumber of gene in combined cds: ", nrow(df_gene))
    gene_count = cbind(gene_count_1, gene_count_2)
    cat("\nnumber of cell in cds 1: ", nrow(df_cell_1))
    cat("\nnumber of cell in cds 2: ", nrow(df_cell_2))
    cat("\nnumber of cell in combined cds: ", ncol(gene_count))
    df_cell_names_1 = names(df_cell_1)
    df_cell_names_1 = names(df_cell_1)
    df_cell_names_2 = names(df_cell_2)
    common_names = df_cell_names_1[df_cell_names_1 %in% df_cell_names_2]
    df_cell_1 = df_cell_1 %>% select(common_names)
    df_cell_2 = df_cell_2 %>% select(common_names)
    df_cell = rbind(df_cell_1, df_cell_2)
    cds = cds_construct(gene_count, df_cell, df_gene)
    return(cds)
}