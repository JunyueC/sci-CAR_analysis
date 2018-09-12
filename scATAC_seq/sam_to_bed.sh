#!/bin/bash
# this sript accept a input folder including sam file, a sample ID file and a output folder, then tranform all sam file in the input folder to bed file, and output the 
# result to the output folder

input_folder=$1
sample_ID=$2
output_folder=$3

echo 
echo "Start sam to bed file transformation.."
module load samtools/1.3
module load bedtools/2.24.0
mkdir $output_folder
for sample in $(cat $sample_ID); do echo transforming $sample; samtools view -bh $input_folder/$sample.sam|bedtools bamtobed -split -i -|sort -k1,1 -k2,2n ->$output_folder/$sample.bed; done
echo "all bed file transformed"
