#!/bin/bash
#SBATCH --job-name=qc
#SBATCH --output=/data/gersbachlab/tda17/Anglen_7276_210923B6/processing/atac_seq/logs/qc_gen.Anglen_7276_210923B6-pe-blacklist-removal.out

source /data/reddylab/software/miniconda2/bin/activate alex
cd /data/gersbachlab/tda17/Anglen_7276_210923B6/processing/atac_seq/Anglen_7276_210923B6-pe-blacklist-removal

python /data/reddylab/software/cwl/bin/generate_stats_atacseq_paired_end.py ./ \
    -samples $(/bin/ls -1 *PBC.txt | sed 's@.PBC.txt@@') \
> qc.txt
