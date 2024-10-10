#!/bin/bash
ORDER=Anglen_7276_210923B6_RNAseq
PROCESSED_DATA_DIR=/home/tda17/gershbach_tda17/Anglen_7276_210923B6_RNAseq/data/rna_seq/processed_raw_reads/${ORDER}
METADATA=/home/tda17/gershbach_tda17/Anglen_7276_210923B6_RNAseq/data/rna_seq/metadata/rna_seq_download_metadata.Anglen_7276_210923B6_RNAseq.txt

python /data/reddylab/software/cwl/GGR-cwl/v1.0/json-generator/run.py \
    -m ${METADATA} \
    -d ${PROCESSED_DATA_DIR} \
    -o /home/tda17/gershbach_tda17/Anglen_7276_210923B6_RNAseq/processing/rna_seq/jsons \
    -t rna-seq \
    --fastq-gzipped \
    --mem 48000 \
    --nthreads 24 \
    --separate-jsons \
    --skip-star-2pass \
