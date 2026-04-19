#!/bin/bash

#SBATCH --job-name=toh5ad
#SBATCH --mem=25GB
#SBATCH --account=tp_2616_fnom_183960
#SBATCH --partition fast
#SBATCH --output=log/toh5ad.out
#SBATCH --error=log/toh5ad.err

module load conda
#conda init
source /shared/software/miniconda/etc/profile.d/conda.sh
conda activate /shared/projects/tp_2616_fnom_183960/conda/envs/PP_r_seuratdisk
Rscript /shared/projects/tp_2616_fnom_183960/Puigdevall_P/Cajal-scRNAseq-2026/day2/3_Data_Integration_Batch_Correction/workdir/p3_export_to_h5ad.R

