This folder includes script for generating peak count matrix from single cell ATAC-seq bam files, and functions for downstream analysis.

reference_generator.sh generate reference bed file from multiple bed files for peak counting.

peak_count.sh generate single cell peak count matrix based on the reference bed file and the single cell bam files.

peak_count_summary_no_report.r summarize and clean the output from peak_count.sh

peak_count_analysis_script includes functions for sci-ATAC-seq analysis (still updating).

