
#!/bin/bash

# this script takes single cell ATAC-seq data and generate the peak count object

# From scATAC-seq data, I am going to accept the folder of the scATAC-seq data (sam files), sample name for each sample and then:

# (1) Combine the scATAC-seq data for the samples.

# (2) use MACS for peak calling.

# (3) Define promoter peaks (the union of the annotated transcription start site (TSS) (Gencode V17) minus 500 base pairs) (same with Hannah)

# (4) Merge the promoter peaks and the MACS peaks, and produce a bed file for the merged peaks,
# and it also include the information (promoter?) and gene name (if it is a promoter).

# (5) For each cell in the provided samples, calculate the number of reads falling into each peaks. Return the sparse matrix with the sample id, peak id and read number; the sample matrix, which is the sample id and sample name; the feature matrix, which is the peak id and peak information (peak location, promoter/not promoter)
main_folder="../../nobackup/170624/scATAC_output/"


sample_list=$main_folder/barcode_samples.txt
sam_folder=$main_folder/splitted_sam/
reads_report_folder=$main_folder/report/
all_output=$main_folder/peak_count
core=20 # define the number of cores for peak counting

merged_folder=$all_output/merged_sam/
macs_folder=$all_output/macs2_output/
reference_folder=$all_output/reference/
peak_count_folder=$all_output/peak_cout/
summary_folder=$all_output/summary_count/

mkdir -p $merged_folder
mkdir -p $macs_folder
mkdir -p $reference_folder
mkdir -p $peak_count_folder
mkdir -p $summary_folder

promoter_length=500
script_folder="/net/shendure/vol1/home/cao1025/analysis_script/ATAC_RNA_coassay_pipe/scATAC_seq/peak_count_sciATACseq/"
gtf_file="/net/shendure/vol1/home/cao1025/reference/gtf_reference/hg19_mm10/rmchr.gencode.v19.chr_patch_hapl_scaff.annotation.gencode.vM12.chr_patch_hapl_scaff.annotation.gtf.gz"
python_use="/net/shendure/vol1/home/cao1025/anaconda3/bin/python3.6"
dis=50 # define the length of entending the reads from the tagmentation site for peak counting
pure_limit=0.85 # the limit for define human and mouse cells

#######********########
# this script is used to accept a sample list, a sam folder, a output folder, and 
# then merge all the sam files in the sample list into a merged sam file
# load the samtools
module load samtools/1.3

merged_file=$merged_folder/merge.sam
echo >$merged_file

echo "start merging the files..."
echo 
# output the combined sam file into the merged folder
for sample in $(cat $sample_list); do echo merging $sample; samtools view $sam_folder/$sample.sam >>$merged_file;done

echo
echo "Merged file is generated~"

#######********########

# This script accept a input sam file, a output folder and call MACS2 for calling peaks
input_file=$merged_file
output_folder=$macs_folder
gene_size="hs"

macs2 callpeak -t $input_file --nomodel --keep-dup all --extsize 200 --shift -100 -B --SPMR -f SAM -g $gene_size --outdir $output_folder --call-summits

#######********########
# Generate promoter file from the gtf reference
#For promoter bed files from UCSC, I am going accept a gtf file and then extract the n = 500bp
#upstream of the gene start site, and output to a bed file with chr, start site, end site, and intersected gene name.
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
peak_bed=$macs_folder/NA_peaks.narrowPeak
output_folder=$reference_folder

merge_file=$output_folder/merged.bed
promoter_file=$output_folder/promoter.bed

# generate a combined bed file
mkdir -p $output_folder
echo Start merging the bed files
# generating the merged files
cat $promoter_bed $peak_bed|sort -k1,1 -k2,2n |cut -f 1,2,3 |bedtools merge -i - >$merge_file

echo Start generating the promoter bed file
# generate the promoter and gene files intersection data file
bedtools intersect -a $merge_file -b $promoter_bed -wa -wb >$promoter_file

echo promoter bed file are generated.
echo Merged file generated.


#######********########
# count reads in each peaks and return the peak count matrix
# Take a sample folder, sample list and a reference bed file, 
# and then generate a sparse matrix for each sample for the read 
# count and also report sparse matrix including row id (DHS site), 
# column id (cell id) and also the read count

ref_bed=$merge_file
out_folder=$peak_count_folder
script=$script_folder/peak_map.py

$python_use $script $sample_list $sam_folder $ref_bed $peak_count_folder $dis $core

#######********########
# combine the peaks count files

echo combining the cell count data...
echo>$out_folder/cell_count.MM
for sample in $(cat $sample_list); do echo combining $sample; cat $sam_folder/$sample.count >>$out_folder/cell_count.MM; rm $sam_folder/$sample.count; done

echo All files are combined.

#######********########
# process the data and generate the R object for further processing
# return a sparse matrix including the peak name (row), cell name (column) 
# and read count in each peak, a peak annotation file with peak name, 
# peak chromatin, start, end, promoter or not; a cell annotation file 
# including cell name, human reads, mouse reads, total reads, human peak 
# reads, mouse peak reads, total peak reads, total peak number

script=$script_folder/peak_count_summary.r
Rscript $script $peak_count_folder $reference_folder $reads_report_folder $summary_folder $pure_limit
echo peak count summary generated.