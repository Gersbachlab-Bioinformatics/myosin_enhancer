#~/anaconda2/envs/Charlie_capture_HiCAR/bin/python
import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import cooltools
import cooler
from matplotlib.ticker import EngFormatter
from matplotlib.colors import LogNorm
bp_formatter = EngFormatter('b')
norm = LogNorm(vmax=200)

def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)

def Draw_heatmap(sampleid, res):
	clr = cooler.Cooler("./aggregated_cool/"+sampleid+".1000.mcool::resolutions/"+str(res))
	print(f'chromosomes: {clr.chromnames}, binsize: {clr.binsize}')
	### to make a list of chromosome start/ends in bins:
	chromstarts = []
	for i in clr.chromnames:
		print(f'{i} : {clr.extent(i)}')
		chromstarts.append(clr.extent(i)[0])

	plt.clf()
	f, ax = plt.subplots(figsize=(7,6))
	start, end = 55_111_190, 55_258_690
	region = ('chr14', start, end)
	norm = LogNorm(vmax=np.quantile(clr.matrix(balance=False).fetch(region),q=0.995))
	im = ax.matshow(clr.matrix(balance=False).fetch(region),extent=(start, end, end, start),norm=norm, cmap="Reds");
	plt.colorbar(im ,ax=ax,fraction=0.046, pad=0.04, label='log10 raw counts '+str(np.quantile(clr.matrix(balance=False).fetch(region),q=0.995)))
	format_ticks(ax)
	plt.tight_layout()
	plt.savefig("./figures/"+sampleid+"_"+res+"b_cool_balanced_max_log_top_0.5_percent.png")


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

#res_list=['1000','2000','5000']
res_list=['10000','2000']
for res in res_list:
	for sampleid in sample_list:
		Draw_heatmap(sampleid,res)


