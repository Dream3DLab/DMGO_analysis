#! /bin/bash

cd /hpc/pmc_rios/1.projects/RA1_scRNA/00.runs/s007

if [ -d /jobs/count_dual ]; then
    echo "../jobs/count_dual exists, overwriting"
else
    echo "creating ../jobs/count_dual"
    mkdir jobs/count_dual/
fi

JOB_DIR="/hpc/pmc_rios/1.projects/RA1_scRNA/00.runs/s007/jobs/count_dual"
echo "${JOB_DIR}"
EXP='ExpIDDI.txt'

while read -r ExpIDDI; do
    echo $ExpIDDI
    RUN=${ExpIDDI}
    JOBS="${JOB_DIR}/${ExpIDDI}_dual.sh"

    echo "#!/bin/bash

#SBATCH --job-name=cellranger_count_${ExpIDDI}
#SBATCH --time=12:00:00
#SBATCH --mem=100G
#SBATCH --mail-type=END,ERROR,FAIL
#SBATCH --mail-user=h.c.r.ariese@prinsesmaximacentrum.nl

export PATH=/hpc/local/Rocky8/pmc_rios/CellRanger/cellranger-8.0.1:\$PATH

cellranger count --id="${ExpIDDI}_analysis" --transcriptome="/hpc/local/Rocky8/pmc_rios/CellRanger/refdata-gex-GRCh38-2024-A" --fastqs="/hpc/pmc_rios/1.projects/RA1_scRNA/00.runs/s007/s007_fastq/outs/fastq_path/HWKVKDMXY" --sample="${ExpIDDI}" --create-bam=true 
" > ${JOBS}

sbatch ${JOBS}

done < $EXP