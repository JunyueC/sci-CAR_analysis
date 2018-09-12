read1="../../nobackup/170701/fastq_ATAC_2/Undetermined_S0_R1_001.fastq.gz"
read2="../../nobackup/170701/fastq_ATAC_2/Undetermined_S0_R2_001.fastq.gz"
output_folder="../../nobackup/170701/scATAC_output_2/"
P7_index="./P7_index.txt"
P5_index="./P5_index.txt"
barcode_file="/net/shendure/vol1/home/cao1025/analysis_script/ATAC_RNA_coassay_pipe/N7_N5_barcode_nextseq.txt"
core=5
cutoff=1000

script_folder="/net/shendure/vol1/home/cao1025/analysis_script/ATAC_RNA_coassay_pipe/scATAC_seq/"
python="/net/shendure/vol1/home/cao1025/anaconda2/bin/python2.7"
index="/net/shendure/vol1/home/cao1025/reference/index/STAR/STAR_hs37d5_mm10/"

dhs_bed="/net/shendure/vol1/home/cao1025/reference/bed_reference/mm10_hg19/dnase_hot_spot/3T3_293T_combined_hotspot_simp.bed"



# read in the read1 and read2, and P5 barcode, P7 barcode and then split the read1 and read2 based on the P5 barcode
# and P7 barcode, and output the splitted reads to the output_folder/fastq files
fastq_folder=$output_folder/fastq
mkdir -p $fastq_folder
echo Start splitting the ATACseq reads
$python $script_folder/ATAC_split_fastq_nextseq.py $read1 $read2 $fastq_folder $P7_index $P5_index

# Trim the reads
echo
echo "Start trimming the files..."
echo $(date)
module load python/2.7.3
module load cutadapt/1.8.3
module load trim_galore/0.4.1
mkdir $output_folder/trimmed_fastq
trimmed_fastq=$output_folder/trimmed_fastq
sample_ID=$output_folder/fastq/sample_ID.txt
mkdir -p $trimmed_fastq
for sample in $(cat $sample_ID); do echo trimming $sample; sem -j $core trim_galore $fastq_folder/$sample*R1*.gz $fastq_folder/$sample*R2*.gz --paired -a CTGTCTCTTATA -a2 CTGTCTCTTATA --three_prime_clip_R1 1 --three_prime_clip_R2 1 -o $trimmed_fastq; done
sem --semaphoretimeout 600
echo "All trimmed file generated."
module unload python/2.7.3


# align the reads with STAR
# this script take the input folder, sample ID, output folder and index as input, then it runs the bowtie2 with local alignment, and align single end read to the index
input_folder=$trimmed_fastq
sample_ID=$output_folder/fastq/sample_ID.txt
STAR_folder=$output_folder/STAR_align
mkdir -p $STAR_folder
echo
echo "Start the alignment using bowtie2"
echo input folder: $input_folder
echo sample ID file: $sample_ID
echo output folder: $STAR_folder
echo index: $index

mkdir -p $output_folder
for sample in $(cat $sample_ID); do echo aligning $sample; STAR --runThreadN $core --outSAMstrandField intronMotif --genomeDir $index --readFilesCommand zcat --readFilesIn $input_folder/$sample*R1*.gz $input_folder/$sample*R2*.gz --outFileNamePrefix $STAR_folder/$sample --genomeLoad LoadAndKeep; done

#remove the index from the memory
STAR --genomeDir $index --genomeLoad Remove
echo "all alignment done"


# This function accept a input bowtie2  folder, a sample ID file and a output folder, and then filter the reads: remove the mitochondrial reads, samtools -F4 -q 30, sort the file and use picard to remove the dupllicate; and output the mapped reads into the mapped folder of the output folder, output the remove duplicates reads into the rm_dup folder of the output folder
input_folder=$STAR_folder
sample_ID=$output_folder/fastq/sample_ID.txt
filtered_folder=$output_folder/filtered_bam
# make the output folder
mkdir -p $filtered_folder
mkdir $filtered_folder/mapped_bam

# filter the files and generate the mapped bam file
for sample in $(cat $sample_ID); do echo "generating mapped bam file" $sample; sem -j $core samtools view -h -F 4 -q 30 $input_folder/$sample*.sam |awk '$3 != "MT" && $3 != "chrM"' -|samtools view -bh -|samtools sort -|samtools view -bh>$filtered_folder/mapped_bam/$sample.bam; done
sem --semaphoretimeout 600
echo All mapped file generated~

# Transform the bam files to sam files
bam_folder=$filtered_folder/mapped_bam
sample_list=$output_folder/fastq/sample_ID.txt
splitted_folder=$output_folder/splitted_sam
# first convert the bam files into sam files
echo convert bam files to sam files...
for sample in $(cat $sample_list); do echo converting $sample; sem -j $core samtools view -h $bam_folder/$sample.bam>$bam_folder/$sample.sam; done
sem --semaphoretimeout 600
echo All bam files are converted to sam files.

echo Start splitting the sam files based on the barcode...
mkdir -p $splitted_folder
bash $script_folder/samfile_split_permuted.sh $bam_folder $sample_list $splitted_folder $barcode_file $cutoff $core
echo split files done

# remove the duplicates with samtools rmdup
sample_list=$output_folder/barcode_samples.txt
rmdup_folder=$output_folder/rmdup_splitted

mkdir -p $rmdup_folder
for sample in $(cat $sample_list); do echo remove duplicates $sample; samtools rmdup $splitted_folder/$sample.sam $rmdup_folder/$sample.bam; done
echo removing duplicates done~

echo "Start calculating the human and mouse reads number..."
bash $script_folder/report_human_mouse_fraction_bam.sh $rmdup_folder $output_folder/barcode_samples.txt $output_folder/report/human_mouse_read_number

echo "Calculate the fragment size distribution...."
bash $script_folder/frag_size.sh $rmdup_folder $output_folder/barcode_samples.txt $output_folder/report/frag_size $core
echo "All fragment size calculation is done."


# calculate the reads number
# this script accept the input parental folder, create a report folder and a read_number sub folder, and for each sample, 
# calculate the read number in fastq, trimmed_fastq, UMI_attached_R2, STAR_alignment, filtered_sam, samfile. after filter
# barcode
parental_folder=$output_folder
sample_ID=$output_folder/fastq/sample_ID.txt
fastq_folder=$parental_folder/fastq
trimmed_folder=$parental_folder/trimmed_fastq
alignment=$STAR_folder
filtered_sam=$filtered_folder/mapped_bam

#split_sam=$parental_folder/splited_sam
report_folder=$parental_folder/report/read_number
echo
echo "Start calculating the reads number..."
#make the report folder
mkdir -p $report_folder
#calculate the read number and output the read number into the report folder
echo sample,total reads,after trimming,all mapped reads, uniquely aligned reads >$report_folder/read_number.csv
for sample in $(cat $sample_ID); do echo calculating $sample; echo $sample,$(expr $(zcat $fastq_folder/$sample*R2*.gz|wc -l) / 4),$(expr $(zcat $trimmed_folder/$sample*R2*.gz|wc -l) / 4), $(samtools view -F 4 -q 30 $alignment/$sample*.sam|wc -l),$(samtools view $filtered_sam/$sample.bam|wc -l) >>$report_folder/read_number.csv; done
echo "Read number calculation is done."

# report the duplication rate
report_folder=$parental_folder/report/duplication_report
sample_list=$output_folder/barcode_samples.txt
split_folder=$output_folder/splitted_sam
rmdup_folder=$output_folder/rmdup_splitted
mkdir -p $report_folder
for sample in $(cat $sample_list); do echo calculating $sample; echo $sample, $(samtools view $split_folder/$sample.sam |wc -l), $(samtools view $rmdup_folder/$sample.bam |wc -l) >>$report_folder/rmdup_report.csv; done

echo "all analysis is done"