#! /bin/bash

#SBATCH --job-name=SCENIC_I
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --ntasks=4
#SBATCH --mail-type=END,ERROR,FAIL
#SBATCH --mail-user=h.c.r.ariese@prinsesmaximacentrum.nl

source /hpc/local/Rocky8/pmc_rios/miniconda3/etc/profile.d/conda.sh #activate your conda, if needed
conda activate pyscenic	

# fill directories below

cd path/to/working/directory

f_loom_path_scenic='/path/to/data/file.loom'

# run grnboost

pyscenic grn f_loom_path_scenic auxilliaries/allTFs_hg38.txt -o output/adj.csv --num_workers 4