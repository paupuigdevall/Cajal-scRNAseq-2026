#!/bin/bash

#SBATCH --job-name=mid_cranger
#SBATCH --mem=82GB
#SBATCH --account=tp_2616_fnom_183960
#SBATCH --partition=fast
#SBATCH --cpus-per-task=16
#SBATCH --output=log/cellRanger_%A_%a.out
#SBATCH --error=log/cellRanger_%A_%a.err


module load cellranger/9.0.1


OUTPUT_BASE=/shared/projects/tp_2616_fnom_183960/Puigdevall_P/cr_out

mkdir -p ${OUTPUT_BASE}
cd ${OUTPUT_BASE}


INPUTLINE=$(head -n ${SLURM_ARRAY_TASK_ID} $1 | tail -n1)
echo $INPUTLINE
RUN_COUNT_SAMPLE=$(echo $INPUTLINE | awk {'print $1'})
fastqs_folder=$(echo $INPUTLINE | awk {'print $2'})
TPATH=/shared/projects/tp_2616_fnom_183960/Puigdevall_P/cr_ref/refdata-gex-GRCh38-2024-A

cellranger count --id=cr_output_${RUN_COUNT_SAMPLE} \
   --fastqs=${fastqs_folder} \
   --transcriptome=${TPATH} \
   --create-bam true \
   --localmem=80 \
   --localcores=16

