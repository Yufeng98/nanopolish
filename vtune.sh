#!/bin/bash
# source /opt/intel/vtune_profiler_2020/amplxe-vars.sh

# vtune -collect hotspots -r result -- <command line>

# vtune -collect hotspots -r result_8_thread -- ~/nanopolish/nanopolish call-methylation -t 8 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:5,000,000-10,000,000" > /z/scratch7/yufenggu/methylation_example/methylation_calls.tsv

# vtune -start-paused -data-limit=0 -collect memory-access -r result_memory_access ~/nanopolish/nanopolish call-methylation -t 8 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:5,000,000-6,000,000" > /z/scratch7/yufenggu/methylation_example/methylation_calls.tsv

# vtune -start-paused -data-limit=0 -collect memory-consumption -r result_memory_consumption ~/nanopolish/nanopolish call-methylation -t 8 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:5,000,000-10,000,000" > /z/scratch7/yufenggu/methylation_example/methylation_calls.tsv

# vtune -start-paused -data-limit=0 -collect uarch-exploration -r result_uarch_exploration ~/nanopolish/nanopolish call-methylation -t 8 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:5,000,000-10,000,000" > /z/scratch7/yufenggu/methylation_example/methylation_calls.tsv

# vtune -start-paused -data-limit=0 -collect hpc-performance -r result_hрc_performance ~/nanopolish/nanopolish call-methylation -t 8 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:5,000,000-10,000,000" > /z/scratch7/yufenggu/methylation_example/methylation_calls.tsv

# vtune -start-paused -data-limit=0 -collect threading -r result_threading ~/nanopolish/nanopolish call-methylation -t 8 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:5,000,000-10,000,000" > /z/scratch7/yufenggu/methylation_example/methylation_calls.tsv

~/nanopolish/nanopolish call-methylation -t 8 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:5,000,000-10,000,000" > /z/scratch7/yufenggu/methylation_example/methylation_calls.tsv

~/nanopolish/nanopolish call-methylation  -t 1 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -s -w "chr20:5,000,000-10,000,000"

~/nanopolish/nanopolish call-methylation  -t 12 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:0-64,444,167" > /z/scratch7/yufenggu/methylation_example/methylation_calls_fp_0_64.tsv

~/nanopolish/nanopolish call-methylation  -t 24 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:0-64,444,167" > /z/scratch7/yufenggu/methylation_example/methylation_calls_fix_0_64_17_15_30.tsv

~/nanopolish/nanopolish call-methylation  -t 12 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:0-64,444,167" > /z/scratch7/yufenggu/methylation_example/methylation_calls_fix_0_64_18_14_30.tsv

~/nanopolish/nanopolish call-methylation  -t 12 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:0-64,444,167" > /z/scratch7/yufenggu/methylation_example/methylation_calls_fix_0_64_16_16_30.tsv

~/nanopolish/nanopolish call-methylation  -t 12 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:0-64,444,167" > /z/scratch7/yufenggu/methylation_example/methylation_calls_fix_0_64_14_18_30.tsv

~/nanopolish/nanopolish call-methylation  -t 12 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:0-64,444,167" > /z/scratch7/yufenggu/methylation_example/methylation_calls_fix_0_64_12_20_30.tsv