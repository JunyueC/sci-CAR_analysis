#!/bin/bash

# this script accept a sam file folder, a sample list, a output folder, a barcode file, then it will run the sam_split.py on each
# sam file and output the splited sam file to the output folder

sam_folder=$1
sample_list=$2
output_folder=$3
barcode_file=$4
cutoff=$5
core=$6

#define the location of the python
python_path="/net/shendure/vol1/home/cao1025/anaconda2/bin/python2.7"
script_path="/net/shendure/vol1/home/cao1025/analysis_script/ATAC_RNA_coassay_pipe/scATAC_seq/"

echo
echo "Start splitting the sam file..."
echo samfile folder: $sam_folder
echo sample list: $sample_list
echo ouput folder: $output_folder
echo barcode file: $barcode_file
echo cutoff value: $cutoff

mkdir $output_folder
for sample in $(cat $sample_list); do echo Now splitting $sample; sem -j+$core $python_path $script_path/sam_split_permuted.py $sam_folder/$sample.sam $barcode_file $output_folder $cutoff; done
sem --semaphoretimeout 1200
cat $output_folder/*sample_list.txt>$output_folder/All_samples.txt
cp $output_folder/All_samples.txt $output_folder/../barcode_samples.txt

# output the report the report/barcode_read_distribution folder
mkdir -p $output_folder/../report/barcode_read_distribution
mv $output_folder/*.txt $output_folder/../report/barcode_read_distribution/
mv $output_folder/*.png $output_folder/../report/barcode_read_distribution/

echo
echo "All sam file splitted."