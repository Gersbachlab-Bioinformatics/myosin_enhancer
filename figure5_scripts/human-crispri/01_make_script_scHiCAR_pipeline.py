#!/usr/bin/python

import os


cpu=12


def make_script(sample_id):
    output1=open('PBS_scHiCAR_pipeline_'+sample_id+'.pbs','w')
    output1.write('#PBS -N scHiCAR_pipeline_'+sample_id+'\n')
    output1.write('#PBS -q workq\n')
    output1.write('#PBS -l nodes=1:ppn='+str(cpu)+'\n')
    output1.write('#PBS -j oe\n')
    output1.write('\n')
    output1.write('# go workdir\n')
    output1.write('cd $PBS_O_WORKDIR\n')
    output1.write('\n')
    output1.write('# run command \n')
    output1.write('sleep 5\n')
    output1.write('\n')
    output1.write('echo -n \"I am on: \"\n')
    output1.write('hostname;\n')    
    output1.write('echo finding ssh-accessible nodes:\n')
    output1.write('echo -n \"running on: \"\n')
    output1.write('\n')

    output1.write("source activate Charlie_capture_HiCAR\n")
    output1.write('#pairtools parse2 -c hg38/Annotation/hg38.chrom.sizes --min-mapq 10 --max-insert-size 2000 --max-inter-align-gap 50 --report-position outer --add-pair-index --no-flip --drop-seq --expand --max-expansion-depth 6 --output-stats ./aligned_pairs/{}.pairsam.stat -o ./aligned_pairs/{}.pairsam.gz ./bam/{}.bam\n'.format(sample_id, sample_id, sample_id))
    output1.write('#pairtools select \'(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")\' --chrom-subset target_chromosome_list.txt -o ./aligned_pairs/{}.selected.pairs.gz --output-rest ./aligned_pairs/{}.unselected.pairs.gz ./aligned_pairs/{}.pairsam.gz\n'.format(sample_id, sample_id, sample_id))
    output1.write('#pairtools restrict -f hg38_CviQI_restricted.bed -o ./aligned_pairs/{}.restrict.pairs.gz ./aligned_pairs/{}.selected.pairs.gz\n'.format(sample_id, sample_id))
    output1.write('#pairtools select "(COLS[-6]==COLS[-3]) and (chrom1==chrom2)" -o ./aligned_pairs/{}.selected.pairs.gz --output-rest ./aligned_pairs/{}.unselected.pairs.gz ./aligned_pairs/{}.restrict.pairs.gz\n'.format(sample_id, sample_id,sample_id))
    output1.write('#pairtools flip -c hg38/Annotation/hg38.chrom.sizes --nproc-in 2 --nproc-out 2 -o ./aligned_pairs/{}.flip.gz ./aligned_pairs/{}.unselected.pairs.gz\n'.format(sample_id,sample_id))
    output1.write('#pairtools sort --tmpdir ./ --nproc 12 --memory 100G -o ./aligned_pairs/{}.sorted.pairs.gz ./aligned_pairs/{}.flip.gz\n'.format(sample_id,sample_id))
    output1.write('#pairtools dedup --max-mismatch 1 --method max -o ./aligned_pairs/{}.dedup.pairs.gz --output-stats ./aligned_pairs/{}.dedup.pairs.stat ./aligned_pairs/{}.sorted.pairs.gz\n'.format(sample_id,sample_id,sample_id))
    output1.write('#pairix ./aligned_pairs/{}.dedup.pairs.gz\n'.format(sample_id))
    output1.write('#/home/schloss/miniconda3/envs/4dn/bin/python ./bin/pairqc/pairsqc.py -p ./aligned_pairs/{}.dedup.pairs.gz -c hg38/Annotation/hg38.chrom.sizes -t P -O ./pairsqc/{} -s {}\n'.format(sample_id,sample_id,sample_id))
    output1.write('#cooler cload pairix --max-split 2 --nproc 12 hg38/Annotation/hg38.chrom.sizes:10000 ./aligned_pairs/{}.dedup.pairs.gz ./aggregated_cool/{}.10000.cool\n'.format(sample_id,sample_id))
    output1.write('#cooler cload pairix --max-split 2 --nproc 12 hg38/Annotation/hg38.chrom.sizes:2000 ./aligned_pairs/{}.dedup.pairs.gz ./aggregated_cool/{}.2000.cool\n'.format(sample_id,sample_id))
    output1.write('#cooler cload pairix --max-split 2 --nproc 12 hg38/Annotation/hg38.chrom.sizes:20000 ./aligned_pairs/{}.dedup.pairs.gz ./aggregated_cool/{}.20000.cool\n'.format(sample_id,sample_id))
    output1.write('cooler cload pairix --max-split 2 --nproc 12 hg38/Annotation/hg38.chrom.sizes:1000 ./aligned_pairs/{}.dedup.pairs.gz ./aggregated_cool/{}.1000.cool\n'.format(sample_id,sample_id))
    output1.write('#cooler balance --cis-only -p 12 ./aggregated_cool/{}.10000.cool\n'.format(sample_id))
    output1.write('#cooler balance --cis-only -p 12 ./aggregated_cool/{}.2000.cool\n'.format(sample_id))
    output1.write('#cooler balance --cis-only -p 12 ./aggregated_cool/{}.20000.cool\n'.format(sample_id))
    output1.write('cooler balance --cis-only -p 12 ./aggregated_cool/{}.1000.cool\n'.format(sample_id))
    output1.write('#cooler zoomify --balance -r 10000N -n 12 -o ./aggregated_cool/{}.10000.mcool ./aggregated_cool/{}.10000.cool\n'.format(sample_id,sample_id))
    output1.write('#cooler zoomify --balance -r 2000N -n 12 -o ./aggregated_cool/{}.2000.mcool ./aggregated_cool/{}.2000.cool\n'.format(sample_id,sample_id))
    output1.write('#cooler zoomify --balance -r 20000N -n 12 -o ./aggregated_cool/{}.20000.mcool ./aggregated_cool/{}.20000.cool\n'.format(sample_id,sample_id))
    output1.write('cooler zoomify --balance -r 1000N -n 12 -o ./aggregated_cool/{}.1000.mcool ./aggregated_cool/{}.1000.cool\n'.format(sample_id,sample_id))
    output1.write('\n')

    output1.write('sleep 30\n')
    output1.write('exit 0')

    output1.close()


sample_list = ['krabnontarget1','krabnontarget2','krabregion31','krabregion32']

for sample_id in sample_list:
    make_script(sample_id)
    os.system('qsub PBS_scHiCAR_pipeline_'+sample_id+'.pbs')



