
# Run MAC2 on the combined reference for human and mouse
# This script accept a input sam file, a output folder and call MACS2 for calling peaks
input_file="./data/bulk_ATAC/ATAC_293T.bed"
output_folder="./data/macs2/293T/"
gene_size="hs" # for mouse, use "mm"
mkdir -p $output_folder
macs2 callpeak -t $input_file --nomodel --keep-dup all --extsize 200 --shift -100 -q 0.1 -B --SPMR -f BED -g $gene_size --outdir $output_folder --call-summits