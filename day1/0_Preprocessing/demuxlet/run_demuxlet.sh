#!/bin/bash

#SBATCH --job-name=demuxlet
#SBATCH --mem=20GB
#SBATCH --account=tp_2616_fnom_183960
#SBATCH --partition=fast
#SBATCH --cpus-per-task=2
#SBATCH --output=log/demuxlet_%A_%a.out
#SBATCH --error=log/demuxlet_%A_%a.err

module load conda
source /shared/software/miniconda/etc/profile.d/conda.sh
conda activate /shared/projects/tp_2616_fnom_183960/conda/envs/PP_demuxlet

INPUTLINE=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n1)
echo $INPUTLINE

INPUTBAM=$(echo $INPUTLINE | awk {'print $1'})
poolnum=$(echo $INPUTLINE | awk {'print $2'})
INPUTVCF=/shared/projects/tp_2616_fnom_183960/Puigdevall_P/midbrainDataset/DataS1_inVitroLines_hg38LO_modified2.vcf
rootfolder=$(echo $INPUTBAM | sed 's/possorted_genome_bam.bam//g')
INPUTBARCODES=$(echo ${rootfolder}filtered_feature_bc_matrix/barcodes.tsv.gz)
OUTPUT=$(echo ${rootfolder}output.demuxlet.doublet0.5.noFilter)
SAMPLELIST=${poolnum}_sample_list.txt

demuxlet --sam $INPUTBAM --vcf $INPUTVCF --field GT --doublet-prior 0.5 --group-list $INPUTBARCODES --out $OUTPUT --sm-list $SAMPLELIST

