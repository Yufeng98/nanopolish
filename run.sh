#!/bash/bin
# ~/nanopolish/nanopolish call-methylation -t 48 -h 1e-4 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:0-64,444,167" 2> result_4.txt

# ~/nanopolish/nanopolish call-methylation -t 48 -h 1e-8 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:0-64,444,167" 2> result_8.txt

# ~/nanopolish/nanopolish call-methylation -t 48 -h 1e-12 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:0-64,444,167" 2> result_12.txt

# ~/nanopolish/nanopolish call-methylation -t 48 -h 1e-16 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:0-64,444,167" 2> result_16.txt

# ~/nanopolish/nanopolish call-methylation -t 48 -h 1e-20 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:0-64,444,167" 2> result_20.txt

# ~/nanopolish/nanopolish call-methylation -t 48 -h 1e-24 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:0-64,444,167" 2> result_24.txt

# ~/nanopolish/nanopolish call-methylation -t 48 -h 1e-28 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:0-64,444,167" 2> result_28.txt

~/nanopolish/nanopolish call-methylation -t 16 -h 40 -r /z/scratch7/yufenggu/methylation_example/output.fastq -b /z/scratch7/yufenggu/methylation_example/output.sorted.bam -g /z/scratch7/yufenggu/methylation_example/reference.fasta -w "chr20:0-64,444,167" > methylation_calls_fix_0_64_32_32_40.tsv

python3 data_process.py 0
python3 data_process.py 4
python3 data_process.py 8
python3 data_process.py 12
python3 data_process.py 16
python3 data_process.py 20
python3 data_process.py 24
python3 data_process.py 28
python3 data_process.py 32
