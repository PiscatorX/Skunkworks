# A Shark eDNA pipeline

## Read QC

First we assess the quality of our reads, this can also be done under DADA2 but fastqc/multiqc visuals are superior. These results can be used to inform the downstreem filtering and trimming of reads. We set a few variables to simplify the pipeline outputs and load the necessary modules.

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


Divisive Amplicon Denoising Algorithm v2 (DADA) was provided to correct amplicon sequencing errors generating amplicon sequence variants (ASV) which have been established to more powereful (better error handling resulting in better measures of microbial diversity) than conventionl OTUs clusteres at 97\%.

### R/DADA2 Installation on the cluster

DADA2 is an R package and to run it on the cluster we have to set it up propally. It also important to point out that DADA2 is implemented as part of the [Qiime2 pipeline](https://docs.qiime2.org/2024.5/plugins/available/dada2/denoise-paired/). However, we are going to run DADA2 directly throug R.

To get R to run on the cluster we have to set up a few items. First we need to load the R module.

```
module load app/R/R-4.0.2

```

Please take note of the version number, the two other versions of R at the time of writing this would not load due to errors.
We also need to tell R to insall packages in the local directory as we don't have permission to write to main R libraries. We run the following command on thcommand line.

```
export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER

```

We can now start the R console by using ``` R ``` command. Once we are in in the R console were can run our install our packages that we need. We can test ifthe packages we need are installed by running ``` > library(dada2) ```

