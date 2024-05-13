#! /bin/bash

#SBATCH --job-name=SCENIC_II-III
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --ntasks=4
#SBATCH --mail-type=END,ERROR,FAIL
#SBATCH --mail-user=h.c.r.ariese@prinsesmaximacentrum.nl

source /hpc/local/Rocky8/pmc_rios/miniconda3/etc/profile.d/conda.sh #activate your conda, if needed
conda activate pyscenic	

cd /hpc/pmc_rios/2.personal/rariese/BRO_t1/Scenic/

f_db_names='/path/to/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /path/to/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /path/to/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather /path/to/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather'
f_motif_path='/path/to/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl'
f_loom_path_scenic='/path/to/data/file.loom'

pyscenic ctx output/adj.csv \
    $f_db_names \
    --annotations_fname $f_motif_path \
    --expression_mtx_fname $f_loom_path_scenic \
    --output output/reg.csv \
    --mask_dropouts \
    --num_workers 4