# A Shark eDNA pipeline

## Read QC

First we assess the quality of our reads, this can also be done under DADA2 but fastqc/multiqc has visuals. These results can be used to inform the downstreem filtering and trimming of reads. We set a few variables to simplify the pipeline outputs and load the necessary modules.

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

## [DADA2](https://bioconductor.org/packages/release/bioc/html/dada2.html) workflow


Divisive Amplicon Denoising Algorithm v2 (DADA) was provided to correct amplicon sequencing errors generating Amplicon sequence variants (ASV) where seen to be more powerful than  conventionl OTUs clusteres at 97\%.


