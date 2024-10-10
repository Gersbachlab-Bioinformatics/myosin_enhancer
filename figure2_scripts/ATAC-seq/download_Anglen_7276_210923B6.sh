#!/bin/bash
METADATA=/data/gersbachlab/tda17/Anglen_7276_210923B6/data/atac_seq/metadata/atac_seq_download_metadata.Anglen_7276_210923B6.txt
DATA_HOME=/data/gersbachlab/tda17/Anglen_7276_210923B6/data/atac_seq
mkdir -p ${DATA_HOME}/raw_reads/

cp -R /gpfs/fs1/data/gersbachlab/tda17/Anglen_7276_210923B6/CHIR.YAP5SAgfp.ATACseq/*fastq.gz /gpfs/fs1/data/gersbachlab/tda17/Anglen_7276_210923B6/data/atac_seq/raw_reads/Anglen_7276_210923B6/
