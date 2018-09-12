
input_folder=$1
sample_ID=$2
output_folder=$3
core=$4

mkdir -p $output_folder
for sample in $(cat $sample_ID); do echo Processing $sample; samtools view $input_folder/$sample.sam | awk -v pat=$sample '{if ($9 > 0) {print $9,pat}}' - >$output_folder/$sample.csv; done