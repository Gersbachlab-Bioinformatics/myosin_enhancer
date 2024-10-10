#~/anaconda2/envs/Charlie_capture_HiCAR/bin/python
import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
from itertools import product

##############################
#target_coord = (23303865,23451365)
target_coord = (55111190, 55258690)
def make_coordinate(coord,res):
	scale_factor=res*1000
	coord_start, coord_end=coord
	start=math.floor(coord_start/scale_factor)*scale_factor
	end=math.ceil(coord_end/scale_factor)*scale_factor
	target_region=[]
	for i in range(int(start/scale_factor),int(end/scale_factor)):
		target_region.append("chr14."+str(i*scale_factor)+"."+str((i+1)*scale_factor))
	return(target_region)

def main(resol,name,stimulation):
	target_coord = (55111190, 55258690)
	draw_list=make_coordinate(target_coord,resol)
	tmp=pd.DataFrame(np.zeros((len(draw_list),len(draw_list))))
	tmp.index=draw_list
	tmp.columns=draw_list
	input_data=pd.read_csv("./cool2txt/HL1_Region3_active_vs_nontarget_"+stimulation+"_Union_"+str(resol)+"kb"+name+"DE_results.csv",index_col=0)
	#### log2Folchange
	for i in input_data.index:
		read1=".".join(i.split(".")[0:3])
		read2=".".join(i.split(".")[3:6])
		#if input_data.loc[i,'padj']<=0.05:
		tmp.loc[read1,read2]=input_data.loc[i,'log2FoldChange']
		tmp.loc[read2,read1]=input_data.loc[i,'log2FoldChange']

	plt.clf()
	fig, ax = plt.subplots(figsize=(10, 8))
	sns.heatmap(tmp,vmax=2,vmin=-2,square=True,cmap="RdBu_r")
	for i, j in product(range(tmp.shape[0]), range(tmp.shape[1])):
		if tmp.index[i]+"."+tmp.columns[j] in input_data.index:
			if input_data.loc[tmp.index[i]+"."+tmp.columns[j],'padj']<=0.05:
				ax.text(i + 0.5,j + 0.5,"*",ha="center",va="center",color="black",fontsize=20,fontweight="bold")
				ax.text(j + 0.5,i + 0.5,"*",ha="center",va="center",color="black",fontsize=20,fontweight="bold")

	plt.tight_layout()
	plt.savefig("./figures/HL1_Region3_active_vs_nontarget_"+stimulation+"_Union_"+str(resol)+"kb"+name+"DE_log2FC_heatmap.pdf")
	
	#### pval
	tmp=pd.DataFrame(np.zeros((len(draw_list),len(draw_list))))
	tmp.index=draw_list
	tmp.columns=draw_list
	for i in input_data.index:
		read1=".".join(i.split(".")[0:3])
		read2=".".join(i.split(".")[3:6])
		if input_data.loc[i,'log2FoldChange']>0:
			indicator=-1
		if input_data.loc[i,'log2FoldChange']<0:
			indicator=1
		tmp.loc[read1,read2]=(indicator*np.log10(input_data.loc[i,'padj']))
		tmp.loc[read2,read1]=(indicator*np.log10(input_data.loc[i,'padj']))

	plt.clf()
	fig, ax = plt.subplots(figsize=(10, 8))
	sns.heatmap(tmp,vmax=3,vmin=-3,square=True,cmap="RdBu_r")
	for i, j in product(range(tmp.shape[0]), range(tmp.shape[1])):
		if tmp.index[i]+"."+tmp.columns[j] in input_data.index:
			if abs(input_data.loc[tmp.index[i]+"."+tmp.columns[j],'log2FoldChange'])>=1:
				ax.text(i + 0.5,j + 0.5,"*",ha="center",va="center",color="black",fontsize=20,fontweight="bold")
				ax.text(j + 0.5,i + 0.5,"*",ha="center",va="center",color="black",fontsize=20,fontweight="bold")

	plt.tight_layout()
	plt.savefig("./figures/Region3_active_vs_nontarget_"+stimulation+"_Union_"+str(resol)+"kb"+name+"DE_adjpval_heatmap.pdf")

main(10,"_",'ctrl' )
main(10,"_",'stim')
