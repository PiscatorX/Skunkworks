# A Shark eDNA pipeline

## Read QC

First we assess the quality of our reads, this informs the downstreem filtering and trimming of reads.

```
rawread_dir="/new-home/andhlovu/Shark_eDNA/RawReads"
output_dir="/new-home/andhlovu/Shark_eDNA/Results"
threads=16

module load app/fastqc app/multiqc

mkdir -p ${output_dir}/rawreads_fastqc

fastqc \
    --threads ${threads} \
    -o ${output_dir}/rawreads_fastqc \
    ${rawread_dir}/*


multiqc \
  ${output_dir}/rawreads_fastqc \
 -o  ${output_dir}/rawreads_multiqc

```
