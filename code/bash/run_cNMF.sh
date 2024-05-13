#! /bin/bash

#SBATCH --job-name=cNMF
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=6
#SBATCH --mail-type=END,ERROR,FAIL
#SBATCH --mail-user=h.c.r.ariese@prinsesmaximacentrum.nl

source /hpc/local/Rocky8/pmc_rios/miniconda3/etc/profile.d/conda.sh #activate your conda, if needed

cd /hpc/pmc_rios/2.personal/rariese/BRO_t1/NMF

conda activate cnmf_env

OUT="Jessa/output_ngenes2000_niter100_h33k27m_all"

# Step 1 - normalize the input matrix and prepare the run parameters
echo "@ Step 1"
python cnmf.py prepare --output-dir . --name $OUT -c Jessa/h33k27m_expression.tsv.gz \
  -k 5 6 7 8 9 \
  --n-iter 100 \
  --total-workers 6 \
  --seed 100 \
  --numgenes 2000 \
  --beta-loss frobenius

# Step 2 - factorize the matrix
# Computationally intensive step -- run in parallel with two workers
echo "@ Step 2"
nohup parallel \
  python cnmf.py factorize --output-dir . --name $OUT \
  --worker-index {} ::: 0 1 2 3 4 5

# Step 3 - combine the results
echo "@ Step 3"
python cnmf.py combine --output-dir . --name $OUT

# remove intermediate files
rm $OUT/cnmf_tmp/$OUT.spectra.k_*.iter_*.df.npz

# Step 4 - select k
echo "@ Step 4"
python cnmf.py k_selection_plot --output-dir . --name $OUT

# Step 5 - get consensus estimates at chosen k
# this step is submitted manually for a chosen value K
#
# K=? 
#
# python ${CNMF}/cnmf.py consensus --output-dir . --name ${OUT} \
#  --components \$K --local-density-threshold 2.00 --show-clustering
#
# python ${CNMF}/cnmf.py consensus --output-dir . --name ${OUT} \
#  --components \$K --local-density-threshold 0.1 --show-clustering

