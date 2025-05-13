#!/bin/bash

#SBATCH --job-name=cellranger_count_BRO_RA1_exp023_d014
#SBATCH --time=12:00:00
#SBATCH --mem=100G
#SBATCH --mail-type=END,ERROR,FAIL
#SBATCH --mail-user=h.c.r.ariese@prinsesmaximacentrum.nl

export PATH=/hpc/local/Rocky8/pmc_rios/CellRanger/cellranger-8.0.1:$PATH

cellranger count --id=BRO_RA1_exp023_d014_analysis --transcriptome=/hpc/local/Rocky8/pmc_rios/CellRanger/refdata-gex-GRCh38-2024-A --fastqs=/hpc/pmc_rios/1.projects/RA1_scRNA/00.runs/s007/s007_fastq/outs/fastq_path/HWKVKDMXY --sample=BRO_RA1_exp023_d014 --create-bam=true 

