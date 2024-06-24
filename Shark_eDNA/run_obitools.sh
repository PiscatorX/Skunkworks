#/bin/bash

set -eu

#set you the directory of your reads here
read_dir="/home/drewx/DevOps/raw_reads"
#This file will contain paired the filenames and sample ids 
manifest="manifest.tsv"
#This is where your outputs will go this
output="/home/drewx/devs/Shark_eDNA"

#use python script to generate paired data file manifest
#https://github.com/piscatorx/blue-carbon-microbiomes/blob/main/scripts/qiime_buildmanifest.py
#this creates a temporary file with headers
qiime_buildmanifest.py \
    -r "${read_dir}" \
    -m  /tmp/manifest.tmp

#this removes the headers in the temporary file 
awk "NR!=1 {print}" /tmp/manifest.tmp  > "${manifest}"

mkdir -p "${output}/merged"

while read sample_id fwd rev
do
    obi import --fastq-input  ${fwd}  obitools3/${sample_id}_1
    obi annotate -S sample:\'${sample_id}\' obitools3/${sample_id}_1  obitools3/${sample_id}_tag1  
    obi import --fastq-input  ${rev}  obitools3/${sample_id}_2
    obi annotate -S sample:\'${sample_id}\' obitools3/${sample_id}_2  obitools3/${sample_id}_tag2  
    obi alignpairedend -R obitools3/${sample_id}_2 obitools3/${sample_id}_1 obitools3/raw_merged_${sample_id}
    obi grep -p "sequence['score_norm'] > 0.8" obitools3/raw_merged_${sample_id} obitools3/filtered_merged_${sample_id}
    obi export --fastq-output obitools3/filtered_merged_${sample_id} > "${output}/merged/${sample_id}.fastq"
    break

done <  ${manifest}

# mkdir -p "${output}/fastqc"
# mkdir -p "${output}/multiqc"
#while read fwd rev
#obi import --fasta {input[0]} {params[0]}/demultiplexed

# obi uniq -m sample {params[0]}/demultiplexed {params[0]}/dereplicated_sequences
# obi annotate -k COUNT -k MERGED_sample {params[0]}/dereplicated_sequences {params[0]}/cleaned_metadata_sequences
# obi grep -p "sequence['COUNT']>=10" {params[0]}/cleaned_metadata_sequences {params[0]}/denoised_sequences
# obi clean -s MERGED_sample -r 0.1 -H {params[0]}/denoised_sequences {params[0]}/cleaned_sequences
# obi export --fasta-output {params[0]}/cleaned_sequences >{output[0]}
