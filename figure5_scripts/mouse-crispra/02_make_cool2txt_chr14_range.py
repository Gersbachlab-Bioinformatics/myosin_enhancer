#!/usr/bin/python

import os


cpu=1

def make_script(sample_id,range_set):
    output1=open('PBS_make_cool2txt_same_CRISPRi_'+sample_id+'.pbs','w')
    output1.write('#PBS -N cool2txt_same_CRISPRi_'+sample_id+'\n')
    output1.write('#PBS -q workq\n')
    output1.write('#PBS -l nodes=n13:ppn='+str(cpu)+'\n')
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
    #output1.write('cooler dump --header --join --no-balance --range {} -o ./cool2txt/{}.1kb.txt.gz ./aggregated_cool/{}.1000.mcool::/resolutions/1000\n'.format(range_set, sample_id, sample_id))
    output1.write('cooler dump --header --join --no-balance --range {} -o ./cool2txt/{}.2kb.txt.gz ./aggregated_cool/{}.1000.mcool::/resolutions/2000\n'.format(range_set, sample_id, sample_id))
    #output1.write('cooler dump --header --join --no-balance --range {} -o ./cool2txt/{}.5kb.txt.gz ./aggregated_cool/{}.1000.mcool::/resolutions/5000\n'.format(range_set, sample_id, sample_id))
    output1.write('cooler dump --header --join --no-balance --range {} -o ./cool2txt/{}.10kb.txt.gz ./aggregated_cool/{}.1000.mcool::/resolutions/10000\n'.format(range_set, sample_id, sample_id))
    #output1.write('cooler dump --header --join --no-balance --range {} -o ./cool2txt/{}.20kb.txt.gz ./aggregated_cool/{}.1000.mcool::/resolutions/20000\n'.format(range_set, sample_id, sample_id))
    #output1.write('cooler dump --header --join --no-balance --range {} -o ./cool2txt/{}.40kb.txt.gz ./aggregated_cool/{}.20000.mcool::/resolutions/40000\n'.format(range_set, sample_id, sample_id))
    output1.write('\n')
    
    output1.write('sleep 30\n')
    output1.write('exit 0')

    output1.close()



#sample_list = ['krabnontarget1','krabnontarget2','krabregion31','krabregion32']
sample_list = [
'HL1_TA432_NT_1_stim',
'HL1_TA432_NT_1_ctrl',
'HL1_TA432_NT_2_stim',
'HL1_TA432_NT_2_ctrl',
'HL1_TA432_NT_3_stim',
'HL1_TA432_NT_3_ctrl',
'HL1_TA446_CA_1_stim',
'HL1_TA446_CA_1_ctrl',
'HL1_TA446_CA_2_stim',
'HL1_TA446_CA_2_ctrl',
'HL1_TA446_CA_3_stim',
'HL1_TA446_CA_3_ctrl'
]

range1='chr14:55,111,190-55,258,690'
for sample_id in sample_list:
    make_script(sample_id,range1)
    os.system('qsub PBS_make_cool2txt_same_CRISPRi_'+sample_id+'.pbs')



