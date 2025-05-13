#! /bin/bash

#SBATCH --job-name=fastq_dual
#SBATCH --time=12:00:00
#SBATCH --mem=240G
#SBATCH --mail-type=END,ERROR,FAIL
#SBATCH --mail-user=h.c.r.ariese@prinsesmaximacentrum.nl


export PATH=/hpc/local/CentOS7/pmc_rios/CellRanger/cellranger-7.0.1:$PATH
export PATH=/hpc/local/CentOS7/pmc_rios/bcl2fastq/bin:$PATH

cd /hpc/pmc_rios/1.projects/RA1_scRNA/00.runs/s007

cellranger mkfastq \
--id=s007_fastq_redone_for_submission \
--run="/hpc/pmc_rios/1.projects/RA1_scRNA/00.runs/s007/240809_A01230_0473_AHWKVKDMXY/240809_A01230_0473_AHWKVKDMXY" \
--csv="/hpc/pmc_rios/1.projects/RA1_scRNA/00.runs/s007/s007_samplesheet_DI.csv" \