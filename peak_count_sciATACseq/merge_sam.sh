
# Next I am going to write a script, that accept a folder of single cell samfiles, a sample list file, 
# and a output_folder, then it combines all the sam files in the sample list file and generate the merged sam file in
# the output folder

#######********########
# this script is used to accept a sample list, a sam folder, a output folder, and 
# then merge all the sam files in the sample list into a merged sam file
# load the samtools
sam_folder="../../nobackup/170624/scATAC_output/rmdup_splitted/"
merged_folder="./data/merged_sam_1"
sample_list="../../nobackup/170624/scATAC_output/barcode_samples.txt"

mkdir -p $merged_folder
module load samtools/1.3
merged_file=$merged_folder/merge.sam
echo >$merged_file
echo "start merging the files..."
echo 
# output the combined sam file into the merged folder
for sample in $(cat $sample_list); do echo merging $sample; samtools view $sam_folder/$sample.bam>>$merged_file;done
echo
echo "Merged file is generated~"