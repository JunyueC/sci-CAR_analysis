
### Aggregate nearby peaks (script from the Hannah)

# here I define a function that accept a vector of peak name in the regression output, then it remove the first "X"
# and then generate a data frame with the chr, start, end location
df_peak_coor2 <- function(peak_names)
    {
    df_peak_name = data.frame(id = peak_names)
    df_peak_name = ((df_peak_name
                     %>% select(id) 
                     %>% separate(id, into = c("chr", "start", "stop"), sep = "-")))
    return(df_peak_name)
}


# accept a cds, and aggregate nearby peak within distance

make_bin_col <- function(cds, distance) {
    coords_string_df <- df_peak_coor2((fData(cds))$peak)
    coords_string_df$start = as.numeric(coords_string_df$start)
    coords_string_df$stop = as.numeric(coords_string_df$stop)
    
    coords_ranges <- GenomicRanges::makeGRangesFromDataFrame(coords_string_df)
    coords_range_merge <- GenomicRanges::reduce(coords_ranges, min.gapwidth = distance)

    merge_df <- data.frame(seqnames=GenomicRanges::seqnames(coords_range_merge),
                           starts=GenomicRanges::start(coords_range_merge),
                           ends=GenomicRanges::end(coords_range_merge))
    
    merge_df$name <- paste(merge_df$seqnames, merge_df$starts, merge_df$ends, sep="-")

    overlaps <- GenomicRanges::findOverlaps(coords_ranges,
                                            coords_range_merge,
                                            select="first")
    overlaps <- as.data.frame(overlaps)

    merge_df <- merge_df[overlaps$overlaps,]
    merge_df$name
}


#' Converts sparse matrix to data.table
#' @import data.table
#' @export

sparse_to_datatable <- function(sparse) {
  dgt_mat <- as(Matrix::t(sparse), "dgTMatrix")
  dt <- data.table::data.table(cell = dgt_mat@Dimnames[[1]][dgt_mat@i+1], site=dgt_mat@Dimnames[[2]][dgt_mat@j+1], val = dgt_mat@x)
  data.table::setkey(dt, site, cell)
  dt
}


#' Make an aggregate count cds by collapsing nearby peaks
#'
#' @param cds A cds object.
#' @param distance The distance under which peaks should be collapsed.
#'
#' @return A cds object with aggregated peaks.
#'
#' @export
#'
#' @examples
#' agg_cds <- aggregate_nearby_peaks(sample_cds, distance = 1000)
aggregate_nearby_peaks <- function(cds, distance = 1000) {
    
    
  fData(cds)$bin <- make_bin_col(cds, distance)
  cds <- cds[!is.na(fData(cds)$bin),]

  exprs_dt <- sparse_to_datatable(Matrix(exprs(cds), sparse = TRUE))
  bin_info <- data.table::data.table(site = (fData(cds))$peak,
                                     bin = fData(cds)$bin)
  data.table::setkey(bin_info, site)
  data.table::setkey(exprs_dt, site)
  exprs_dt <- merge(exprs_dt, bin_info)

  data.table::setkey(exprs_dt, cell, bin)
  genomic_bins <- exprs_dt[,sum(val), by="cell,bin"]
  out <- Matrix::sparseMatrix(j=as.numeric(factor(genomic_bins$cell)),
                              i=as.numeric(factor(genomic_bins$bin)),
                              x=genomic_bins$V1)

  fdf <- data.frame(dhs = levels(factor(genomic_bins$bin)),
                    row.names = levels(factor(genomic_bins$bin)))
  pdf <- data.frame(cells = levels(factor(genomic_bins$cell)),
                    row.names = levels(factor(genomic_bins$cell)))
  fdf$bin <- NULL
  #pdf <- pdf[row.names(pData(cds)),]
  pdf <- cbind(pdf, pData(cds)[rownames(pdf), ])
  pdf$pdf <- NULL

  fd <- new("AnnotatedDataFrame", data = fdf)
  pd <- new("AnnotatedDataFrame", data = pdf)

  if (class(exprs(cds)) == "dgCMatrix") {
    compart_cds <-  newCellDataSet(as(out, "sparseMatrix"),
                                    phenoData = pd,
                                    featureData = fd,
                                    expressionFamily=negbinomial.size(),
                                    lowerDetectionLimit=1)
  } else {
    compart_cds <-  newCellDataSet(as.matrix(out),
                                   phenoData = pd,
                                   featureData = fd,
                                   expressionFamily=negbinomial.size(),
                                   lowerDetectionLimit=1)
  }

  return(compart_cds)
}