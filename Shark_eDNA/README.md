# A Shark eDNA pipeline

## Read QC

First we assess the quality of our reads, this can also be done under DADA2 but fastqc/multiqc visuals are superior. These results can be used to inform the downstream filtering and trimming of reads. We set a few variables to simplify the pipeline outputs and load the necessary modules.

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


Divisive Amplicon Denoising Algorithm v2 (DADA) was provided to correct amplicon sequencing errors generating amplicon sequence variants (ASV) which have been established to more powerful (better error handling resulting in better measures of microbial diversity) than conventional OTUs clusters at 97\%.

### R/DADA2 Installation on the cluster

DADA2 is an R package and to run it on the cluster we have to set it up properly. It also important to point out that DADA2 is implemented as part of the [Qiime2 pipeline](https://docs.qiime2.org/2024.5/plugins/available/dada2/denoise-paired/). However, we are going to run DADA2 directly through R.

To get R to run on the cluster we have to set up a few items. First we request to request to used the cluster interactively by issuing the command ```qsubi``` command. More on that and related HP2 commands [here](https://www0.sun.ac.za/hpc/index.php?title=HOWTO_submit_jobs). To use the R console, we can load the R module.

```
module load  app/R/4.3.2

```

I have noticed that sometimes this will generate an error 
We also need to tell R to install packages in the local directory as we don't have permission to write to main R libraries. We run the following command on the command line.

```
export R_LIBS_USER=$HOME/apps/R:$R_LIBS_USER

```

We can now start the R console by using ``` R ``` command. Once we are on the R console were can run our install our packages that we need. We can not run R using the console, for example, we can test if dada2 has been installed by running ``` library(dada2) ```

DADA2 documentation is available from [biconductor](https://bioconductor.org/packages/release/bioc/vignettes/dada2/inst/doc/dada2-intro.html)  and there is also a recommended [tutorial](https://benjjneb.github.io/dada2/tutorial.html).

A [Jupyter notebook](DADA2_workflow.ipynb)  is provided showing dada workflow up to ASV and taxa table. To run the pipeline on the command line on the cluster we invoke inside a script where we request resources using ```#PBS``` directives including the modules required. DADA2 auto detects the available number of cores/cpus on the clusters and takes advantage of these resources to complete the analysis fast. We use the ```Rscript ``` to interpret and run the R script, this can also be achieved using ```R CMD batch``` see the discussion on the differences between using these two [here](https://stackoverflow.com/questions/21969145/why-or-when-is-rscript-or-littler-better-than-r-cmd-batch). The [R script](DADA2_workflow.R) run on the cluster  with dada2 workflow up to th taxonomic assignment is also provided.

```
#!/bin/bash
#PBS -N DADA2-workflow
#PBS -l select=1:ncpus=24
#PBS -l walltime=48:00:00

module load  app/R/4.3.2

#ensures that output files go to the working directory of the script

cd $PBS_O_WORKDIR

Rscript DADA2_workflow.R
```









