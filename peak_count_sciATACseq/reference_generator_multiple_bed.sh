
# in this script, I am going to accept a folder of mutliple bed files, and a bed file sample id,
# then it will call peaks for each bed file, and then merge all the bed files with the promoter bed files

merged_folder="./data/merged_sam/human.merge.bed"
sample_ID=""
macs_folder="./data/macs_output"
gtf_file="/net/shendure/vol1/home/cao1025/reference/gtf_reference/hg19_GRCh37.p13/rmchr.gencode.v19.chr_patch_hapl_scaff.annotation.gtf"
reference_folder="./data/reference_folder"

promoter_length=500
script_folder="/net/shendure/vol1/home/cao1025/analysis_script/ATAC_RNA_coassay_pipe/scATAC_seq/peak_count_sciATACseq/"
python_use="/net/shendure/vol1/home/cao1025/anaconda3/bin/python3.6"

mkdir -p $reference_folder

# accept a merged sam file, a gtf file, and a output folder, and then generate the reference bed file for ATAC-seq count
#######********########

# This script accept a input sam file, a output folder and call MACS2 for calling peaks
input_folder=$merged_folder
output_folder=$macs_folder
gene_size="hs"

mkdir -p $output_folder
for sample in $(cat $sample_ID); do echo calling peaks $sample; macs2 callpeak -t $input_folder/$sample.bed --nomodel --keep-dup all --extsize 200 --shift -100 -q 0.1 -B --SPMR -f BED -g $gene_size --outdir $output_folder --call-summits -n $sample; done

#######********########
# I am going accept a gtf file and then extract the n = 500bp
# upstream of the gene start site, and output to a bed file with chr, start site, end site, and intersected gene name.
promoter_file=$reference_folder/promoter_gtf.bed
script=$script_folder/gtf_to_promoter_bed.r
Rscript $script $gtf_file $promoter_file $promoter_length

#######********########
# merge promoter peak
# (1) generate a merged bed file.
# (2) Generate the intersection between the merged bed file 
# and the promoter bed file: for each 
#gene, define the merged promoter region.

promoter_bed=$promoter_file
output_folder=$reference_folder

merge_file=$output_folder/merged.bed
promoter_file=$output_folder/promoter.bed

# generate a combined bed file
mkdir -p $output_folder
echo Start merging the bed files
# generating the merged files
cat $promoter_bed $macs_folder/*peaks.narrowPeak | cut -f 1,2,3 | sortBed -i - |bedtools merge -i - >$merge_file

echo Start generating the promoter bed file
# generate the promoter and gene files intersection data file
bedtools intersect -a $merge_file -b $promoter_bed -wa -wb >$promoter_file

echo promoter bed file are generated.
echo Merged file generated.