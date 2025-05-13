#! /bin/bash

#SBATCH --job-name=fastqc
#SBATCH --time=40:00:00
#SBATCH --mem=30G
#SBATCH --mail-type=END,ERROR,FAIL
#SBATCH --mail-user=h.c.r.ariese@prinsesmaximacentrum.nl

cd /hpc/local/Rocky8/pmc_rios/FastQC 

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

for FILE in /hpc/pmc_rios/1.projects/RA1_scRNA/00.runs/s007/s007_fastq/outs/fastq_path/HWKVKDMXY/*; 
do 
    echo "Processing $FILE"
    ./fastqc -t 4 -f fastq -o /hpc/pmc_rios/1.projects/RA1_scRNA/00.runs/s007/s007_fastq "$FILE"; 
done

