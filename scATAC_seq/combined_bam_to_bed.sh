
input_folder=$1
sample_ID=$2
output_folder=$3
output_file=$output_folder/combined.bed

mkdir $output_folder
cat >$output_file
for sample in $(cat $sample_ID); do echo combine $sample; samtools view -bh $input_folder/$sample.bam | bedtools bamtobed -split -i - >>$output_file; done
echo combining done
cat $output_file |sort -k1,1 -k2,2n > $output_folder/combined.sorted.bed