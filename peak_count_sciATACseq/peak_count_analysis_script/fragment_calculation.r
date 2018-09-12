
# Here I define a funciton that accept a report folder and sample name list, and then extract the fragment size 
# correlated with the samples
fragment_size_extraction_list <- function(folder, sample_list)
    {
    fragment_size_extraction_file <- function(folder, sample)
        {
        frag_file = paste(folder, "frag_size/", sample, ".csv", sep = "")
        frag_file = read.csv(frag_file, sep=" ", header=F)
        return(frag_file$V1)
    }
    result = c()
    frags = sapply(sample_list, function(sample) {fragment_size_extraction_file(folder, sample)})
    return(as.numeric(unlist(frags)))
}