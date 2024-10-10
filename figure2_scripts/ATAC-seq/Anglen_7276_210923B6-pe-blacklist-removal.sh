#!/bin/bash
#SBATCH --job-name=cwl_atac_seq
#SBATCH --output=/data/gersbachlab/tda17/Anglen_7276_210923B6/processing/atac_seq/logs/Anglen_7276_210923B6-pe-blacklist-removal-%a.out
#SBATCH --mail-user=tda17@duke.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16


export PATH="/data/reddylab/software/bin:$PATH"
export PATH="/data/common/shared_conda_envs/ucsc/bin:$PATH"
export PATH="/data/reddylab/software/cwl/bin:$PATH"
export PATH="/data/reddylab/software/preseq_v2.0:$PATH"
export PATH="/data/reddylab/software/rsem-1.2.21/:$PATH"
export PATH="/data/reddylab/software/phantompeakqualtools-1.2/:$PATH"
export PATH="/data/reddylab/software/miniconda2/envs/cwl10/bin:$PATH"

module load bedtools2
module load fastqc
module load samtools
module load bowtie2
module load java

# For Fastqc
export DISPLAY=:0.0

# Make sure temporary files and folders are created in a specific folder
mkdir -p /data/gersbachlab/tda17/Anglen_7276_210923B6/processing/atac_seq/tmpdirs/tmp-Anglen_7276_210923B6-pe-blacklist-removal-${SLURM_ARRAY_TASK_ID}-
export TMPDIR="/data/gersbachlab/tda17/Anglen_7276_210923B6/processing/atac_seq/tmpdirs/tmp-Anglen_7276_210923B6-pe-blacklist-removal-${SLURM_ARRAY_TASK_ID}-"

cwltool --debug \
    --non-strict \
    --preserve-environment PATH \
    --preserve-environment DISPLAY \
    --preserve-environment TMPDIR \
    --outdir /data/gersbachlab/tda17/Anglen_7276_210923B6/processing/atac_seq/Anglen_7276_210923B6-pe-blacklist-removal  \
    --no-container \
    /data/reddylab/software/cwl/GGR-cwl/v1.0/ATAC-seq_pipeline/pipeline-pe-blacklist-removal.cwl \
    /data/gersbachlab/tda17/Anglen_7276_210923B6/processing/atac_seq/jsons/atac_seq_download_metadata.Anglen_7276_210923B6-pe-blacklist-removal-${SLURM_ARRAY_TASK_ID}.json


# Delete any tmpdir not removed by cwltool
rm -rf /data/gersbachlab/tda17/Anglen_7276_210923B6/processing/atac_seq/tmpdirs/tmp-Anglen_7276_210923B6-pe-blacklist-removal-${SLURM_ARRAY_TASK_ID}-
