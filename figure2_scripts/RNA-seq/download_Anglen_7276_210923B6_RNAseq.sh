#!/bin/bash
METADATA=/home/tda17/gershbach_tda17/Anglen_7276_210923B6_RNAseq/data/rna_seq/metadata/rna_seq_download_metadata.Anglen_7276_210923B6_RNAseq.txt
DATA_HOME=/home/tda17/gershbach_tda17/Anglen_7276_210923B6_RNAseq/data/rna_seq
mkdir -p ${DATA_HOME}/raw_reads/

cp -R /home/tda17/gershbach_tda17/Anglen_7276_210923B6_RNAseq/Anglen_7276_210923B6_RNAseq/*fastq.gz /gpfs/fs1/data/gersbachlab/Anglen_7276_210923B6_RNAseq/Anglen_7276_210923B6_RNAseq/Anglen_7276_210923B6/
