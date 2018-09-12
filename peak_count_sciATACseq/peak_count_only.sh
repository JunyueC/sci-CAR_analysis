
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
reference_folder=""

sample_list=$main_folder/barcode_samples.txt
sam_folder=$main_folder/rmdup_splitted/
reads_report_folder=$main_folder/report/
all_output=$main_folder/peak_count
peak_count_folder=$all_output/peak_cout/
summary_folder=$all_output/summary_count/
merge_file=$reference_folder/merged.bed

core=10

mkdir -p $peak_count_folder
mkdir -p $summary_folder

script_folder="/net/shendure/vol1/home/cao1025/analysis_script/ATAC_RNA_coassay_pipe/scATAC_seq/peak_count_sciATACseq/"
python_use="/net/shendure/vol1/home/cao1025/anaconda3/bin/python3.6"
dis=50 # define the length of entending the reads from the tagmentation site for peak counting
pure_limit=0.85 # the limit for define human and mouse cells


#######********########
# count reads in each peaks and return the peak count matrix
# Take a sample folder, sample list and a reference bed file, 
# and then generate a sparse matrix for each sample for the read 
# count and also report sparse matrix including row id (DHS site), 
# column id (cell id) and also the read count

ref_bed=$merge_file
out_folder=$peak_count_folder
script=$script_folder/peak_map.py
mkdir -p $out_folder
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