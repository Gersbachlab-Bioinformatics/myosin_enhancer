#!/bin/bash
#SBATCH --array=0-9%20
ORDER=Anglen_7276_210923B6_RNAseq
RAW_DATA_DIR=/home/tda17/gershbach_tda17/Anglen_7276_210923B6_RNAseq/data/rna_seq/raw_reads/${ORDER}
PROCESSED_DATA_DIR=/home/tda17/gershbach_tda17/Anglen_7276_210923B6_RNAseq/data/rna_seq/processed_raw_reads/${ORDER}
METADATA=/home/tda17/gershbach_tda17/Anglen_7276_210923B6_RNAseq/data/rna_seq/metadata/rna_seq_download_metadata.Anglen_7276_210923B6_RNAseq.txt

mkdir -p ${PROCESSED_DATA_DIR}
cd ${PROCESSED_DATA_DIR}

seq_name_header=$(/bin/grep -Eoi "sequencing.?core.?library.?name" ${METADATA})
if [[ $? == 1 ]];
then
    echo -e "ERROR: Sequencing core library name not found in ${METADATA}"
    exit 1
fi

name_header=$(/bin/grep -Poi "\tname\t" ${METADATA})
if [[ $? == 1 ]];
then
    echo -e "ERROR: Library Name column not found in ${METADATA}"
    exit 1
fi
name_header=$(echo ${name_header} | cut -f2)

seq_type_header=$(head -1 ${METADATA} | /bin/grep -Poi "paired.?end.?or.?single.?end")
if [[ $? == 1 ]];
then
    echo -e "ERROR: Paired-end or single-end column not found in ${METADATA}"
    exit 1
fi

sample_seq_name=$(/data/reddylab/software/bin/print_tab_cols.awk -v cols="${seq_name_header}" ${METADATA} \
    | awk -v SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} 'NR==SLURM_ARRAY_TASK_ID+1{print}');
sample_name=$(/data/reddylab/software/bin/print_tab_cols.awk -v cols="${name_header}" ${METADATA} \
    | awk -v SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} 'NR==SLURM_ARRAY_TASK_ID+1{print}');
seq_type=$(/data/reddylab/software/bin/print_tab_cols.awk -v cols="${seq_type_header}" ${METADATA} \
    | awk -v SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} 'NR==SLURM_ARRAY_TASK_ID+1{print}');


for read_pair in R1 R2;
do
    sample_files=$(/bin/ls ${RAW_DATA_DIR}/${sample_seq_name/ /}_${read_pair}_[0-9][0-9][0-9]* 2> /dev/null)
    if [[ $? != 0 ]]; # If no samples found with that read_pair, continue
    then
        continue;
    fi
    if [[ ${read_pair} == "R1" || (${seq_type/ /} == "PE" || ${seq_type/ /} == "pe") ]];
    then
        # Merge all lanes
        merged=$(basename $(echo ${sample_files} | awk '{print $1}') | sed -e 's/_[0-9]\{3\}_/_/')
        cat ${sample_files} > ${merged};

        # Rename samples with our sample Names
        dest_filename=$(basename $(echo ${merged} | awk '{print $1}') | sed -r 's/\_S[0-9]+//; s/\_(R1|R2|UMI)\_/\.\1\./; s/\.[0-9]+\.fastq/\.fastq/')
        mv ${merged} ${dest_filename}

        cleaned_dest_filename=${dest_filename/${sample_seq_name/ /}/${sample_name/ /}}

        if [[ ${seq_type/ /} == "SE" || ${seq_type/ /} == "se" ]];
        then
            cleaned_dest_filename=${cleaned_dest_filename/.R1/}
        fi
        
        mv ${dest_filename} ${cleaned_dest_filename}
    fi
done




