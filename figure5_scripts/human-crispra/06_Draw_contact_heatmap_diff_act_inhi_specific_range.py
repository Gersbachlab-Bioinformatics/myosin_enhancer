#~/anaconda2/envs/Charlie_capture_HiCAR/bin/python
import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
from itertools import product

##############################
target_coord = (23303865,23451365)
def make_coordinate(coord,res):
	scale_factor=res*1000
	coord_start, coord_end=coord
	start=math.floor(coord_start/scale_factor)*scale_factor
	end=math.ceil(coord_end/scale_factor)*scale_factor
	target_region=[]
	for i in range(int(start/scale_factor),int(end/scale_factor)):
		target_region.append("chr14."+str(i*scale_factor)+"."+str((i+1)*scale_factor))
	return(target_region)

def main(resol,name,compare, var_set):
	target_coord = (23303865,23451365)
	draw_list=make_coordinate(target_coord,resol)
	tmp=pd.DataFrame(np.zeros((len(draw_list),len(draw_list))))
	tmp.index=draw_list
	tmp.columns=draw_list
	pos, neg=var_set
	input_data_inh=pd.read_csv("./Read_count_matrix/"+compare+"_Union_in_"+pos+"_"+str(resol)+"kb"+name+"DE_results.csv",index_col=0)
	input_data_act=pd.read_csv("./Read_count_matrix/"+compare+"_Union_in_"+neg+"_"+str(resol)+"kb"+name+"DE_results.csv",index_col=0)
	#### add inhi log2Folchange
	for i in input_data_inh.index:
		read1=".".join(i.split(".")[0:3])
		read2=".".join(i.split(".")[3:6])
		#if input_data.loc[i,'padj']<=0.05:
		tmp.loc[read1,read2]=input_data_inh.loc[i,'log2FoldChange']
		tmp.loc[read2,read1]=input_data_inh.loc[i,'log2FoldChange']
	#### subtract act log2Foldchange
	for i in input_data_act.index:
		read1=".".join(i.split(".")[0:3])
		read2=".".join(i.split(".")[3:6])
		tmp.loc[read1,read2]=tmp.loc[read1,read2]-input_data_act.loc[i,'log2FoldChange']
		tmp.loc[read2,read1]=tmp.loc[read2,read1]-input_data_act.loc[i,'log2FoldChange']
	
	plt.clf()
	fig, ax = plt.subplots(figsize=(10, 8))
	sns.heatmap(tmp,vmax=2,vmin=-2,square=True,cmap="RdBu_r")
	"""
	for i, j in product(range(tmp.shape[0]), range(tmp.shape[1])):
		if tmp.index[i]+"."+tmp.columns[j] in input_data.index:
			if input_data.loc[tmp.index[i]+"."+tmp.columns[j],'padj']<=0.05:
				ax.text(i + 0.5,j + 0.5,"*",ha="center",va="center",color="black",fontsize=20,fontweight="bold")
				ax.text(j + 0.5,i + 0.5,"*",ha="center",va="center",color="black",fontsize=20,fontweight="bold")
	"""
	plt.tight_layout()
	plt.savefig("./figures/Difference_"+compare+"_in_"+pos+"_vs_"+compare+"_in_"+neg+"_Union_"+str(resol)+"kb"+name+"DE_log2FC_heatmap.pdf")
	

main(10,"_","TA448_vs_NT",("stim","ctrl"))
main(10,"_","stim_vs_ctrl",("TA448","NT"))
