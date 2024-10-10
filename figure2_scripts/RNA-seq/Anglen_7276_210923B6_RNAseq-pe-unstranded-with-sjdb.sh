#!/bin/bash
#SBATCH --job-name=cwl_rna_seq
#SBATCH --output=/home/tda17/gershbach_tda17/Anglen_7276_210923B6_RNAseq/processing/rna_seq/logs/Anglen_7276_210923B6_RNAseq-pe-unstranded-with-sjdb-%a.out
#SBATCH --mail-user=tda17@duke.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --mem=48000
#SBATCH --cpus-per-task=24


export PATH="/data/reddylab/software/bin:$PATH"
export PATH="/data/common/shared_conda_envs/ucsc/bin:$PATH"
export PATH="/data/reddylab/software/cwl/bin:$PATH"
export PATH="/data/reddylab/software/preseq_v2.0:$PATH"
export PATH="/data/reddylab/software/rsem-1.2.21/:$PATH"
export PATH="/data/reddylab/software/STAR-STAR_2.4.1a/bin/Linux_x86_64/:$PATH"
export PATH="/data/reddylab/software/subread-1.4.6-p4-Linux-x86_64/bin/:$PATH"
export PATH="/data/reddylab/software/bamtools-2.2.3/bin/:$PATH"

export PATH="/data/reddylab/software/miniconda2/envs/cwl10/bin:$PATH"

module load bedtools2
module load fastqc
module load samtools
module load bowtie2
module load java

# For Fastqc
export DISPLAY=:0.0

# Make sure temporary files and folders are created in a specific folder
mkdir -p /home/tda17/gershbach_tda17/Anglen_7276_210923B6_RNAseq/processing/rna_seq/tmpdirs/tmp-Anglen_7276_210923B6_RNAseq-pe-unstranded-with-sjdb-${SLURM_ARRAY_TASK_ID}-
export TMPDIR="/home/tda17/gershbach_tda17/Anglen_7276_210923B6_RNAseq/processing/rna_seq/tmpdirs/tmp-Anglen_7276_210923B6_RNAseq-pe-unstranded-with-sjdb-${SLURM_ARRAY_TASK_ID}-"

cwltool --debug \
    --non-strict \
    --preserve-environment PATH \
    --preserve-environment DISPLAY \
    --preserve-environment TMPDIR \
    --outdir /home/tda17/gershbach_tda17/Anglen_7276_210923B6_RNAseq/processing/rna_seq/Anglen_7276_210923B6_RNAseq-pe-unstranded-with-sjdb  \
    --no-container \
    /data/reddylab/software/cwl/GGR-cwl/v1.0/RNA-seq_pipeline/pipeline-pe-unstranded-with-sjdb.cwl \
    /home/tda17/gershbach_tda17/Anglen_7276_210923B6_RNAseq/processing/rna_seq/jsons/rna_seq_download_metadata.Anglen_7276_210923B6_RNAseq-pe-unstranded-with-sjdb-${SLURM_ARRAY_TASK_ID}.json


# Delete any tmpdir not removed by cwltool
rm -rf /home/tda17/gershbach_tda17/Anglen_7276_210923B6_RNAseq/processing/rna_seq/tmpdirs/tmp-Anglen_7276_210923B6_RNAseq-pe-unstranded-with-sjdb-${SLURM_ARRAY_TASK_ID}-
