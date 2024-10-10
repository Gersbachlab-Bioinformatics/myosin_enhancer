#!/bin/bash
ORDER=Anglen_7276_210923B6
PROCESSED_DATA_DIR=/data/gersbachlab/tda17/Anglen_7276_210923B6/data/atac_seq/processed_raw_reads/${ORDER}
METADATA=/data/gersbachlab/tda17/Anglen_7276_210923B6/data/atac_seq/metadata/atac_seq_download_metadata.Anglen_7276_210923B6.txt

python /data/reddylab/software/cwl/GGR-cwl/v1.0/json-generator/run.py \
    -m ${METADATA} \
    -d ${PROCESSED_DATA_DIR} \
    -o /data/gersbachlab/tda17/Anglen_7276_210923B6/processing/atac_seq/jsons \
    -t atac-seq \
    --fastq-gzipped \
    --mem 24000 \
    --nthreads 16 \
    --separate-jsons

