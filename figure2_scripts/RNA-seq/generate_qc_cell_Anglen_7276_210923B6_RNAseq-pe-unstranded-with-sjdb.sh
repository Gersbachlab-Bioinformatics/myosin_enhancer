#!/bin/bash
#SBATCH --job-name=qc
#SBATCH --output=/home/tda17/gershbach_tda17/Anglen_7276_210923B6_RNAseq/processing/rna_seq/logs/qc_gen.Anglen_7276_210923B6_RNAseq-pe-unstranded-with-sjdb.out

source /data/reddylab/software/miniconda2/bin/activate alex
cd /home/tda17/gershbach_tda17/Anglen_7276_210923B6_RNAseq/processing/rna_seq/Anglen_7276_210923B6_RNAseq-pe-unstranded-with-sjdb

python /data/reddylab/software/cwl/bin/generate_stats_rnaseq_paired_end.py ./ \
    -samples $(/bin/ls -1 *PBC.txt | sed 's@.PBC.txt@@') \
> qc.txt
