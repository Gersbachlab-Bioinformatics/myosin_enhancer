#~/anaconda2/envs/Charlie_capture_HiCAR/bin/python
import os
import pandas as pd
import numpy as np
import math
import glob


def name_func(input_data):
	#if "chr" not in input_data['chrom1']:
	#	input_data['chrom1']="chr"+input_data['chrom1']
	#	input_data['chrom2']="chr"+input_data['chrom2']
	tmp1=input_data[input_data['chrom1']=='chr14']
	tmp2=tmp1[tmp1['chrom2']=='chr14']
	tmp2['region1']=tmp2['chrom1']+"."+tmp2['start1'].astype('str')+"."+tmp2['end1'].astype('str')
	tmp2['region2']=tmp2['chrom2']+"."+tmp2['start2'].astype('str')+"."+tmp2['end2'].astype('str')
	name_vec=tmp2['region1']+"."+tmp2['region2']
	tmp2.index=name_vec
	#tmp2_re=tmp2.loc[[x for x in tmp2.index if tmp2.loc[x,'region1']!=tmp2.loc[x,'region2']]]
	return(tmp2)


def coordinate_convert(coord,res,length):
	scale_factor=res*1000
	chr_num, coord_start, coord_end=coord
	start=(math.floor(coord_start/scale_factor)-length)*scale_factor
	end=(math.ceil(coord_end/scale_factor)+length)*scale_factor
	target_region=[]
	for i in range(int(start/scale_factor),int(end/scale_factor)):
		target_region.append(chr_num+"."+str(i*scale_factor)+"."+str((i+1)*scale_factor))
	return(target_region)


def make_matrix_for_target_region(res, var_set, fix):
	## read data
	feature_list=glob.glob("./cool2txt/*"+fix+"*"+str(res)+"kb"+"*"+"gz")
	Exp, ctrl=var_set
	Exp_list=[x for x in feature_list if Exp in x]
	print(len(Exp_list))
	ctrl_list=[x for x in feature_list if ctrl in x]
	print(len(ctrl_list))
	## save in dict
	Exp_data_dict={Exp+"_"+x.split("_")[2] : name_func(pd.read_csv(x,compression='gzip', header=0, sep='\t')) for x in Exp_list}
	ctrl_data_dict={ctrl+"_"+x.split("_")[2] : name_func(pd.read_csv(x,compression='gzip', header=0, sep='\t')) for x in ctrl_list}	
	### get nontarget and region3 union interaction
	Exp_tmp=[]
	for key in Exp_data_dict.keys():
		if len(Exp_tmp)==0:
			Exp_tmp=Exp_data_dict[key].index
		else:
			Exp_tmp=Exp_tmp.intersection(Exp_data_dict[key].index)	
	ctrl_tmp=[]
	for key in ctrl_data_dict.keys():
		if len(ctrl_tmp)==0:
			ctrl_tmp=ctrl_data_dict[key].index
		else:
			ctrl_tmp=ctrl_tmp.intersection(ctrl_data_dict[key].index)
	total_name=Exp_tmp.union(ctrl_tmp)
	### Construct matirx
	Exp_columns_dict={x:Exp+"_"+x.split("_")[1] for x in sorted(Exp_data_dict.keys())}
	ctrl_columns_dict={x:ctrl+"_"+x.split("_")[1] for x in sorted(ctrl_data_dict.keys())}
	column_list=list(Exp_columns_dict.values())+list(ctrl_columns_dict.values())
	output_Data=pd.DataFrame(np.zeros((len(total_name),len(column_list))),index=total_name,columns=column_list)
	for key in Exp_columns_dict.keys():
		output_Data.loc[total_name.intersection(Exp_data_dict[key].index),Exp_columns_dict[key]]=Exp_data_dict[key].loc[total_name.intersection(Exp_data_dict[key].index),'count']
	for key in ctrl_columns_dict.keys():
		output_Data.loc[total_name.intersection(ctrl_data_dict[key].index),ctrl_columns_dict[key]]=ctrl_data_dict[key].loc[total_name.intersection(ctrl_data_dict[key].index),'count']
	if not os.path.exists(os.getcwd()+"/Read_count_matrix/"+Exp+"_vs_"+ctrl+"/"):
		os.makedirs(os.getcwd()+"/Read_count_matrix/"+Exp+"_vs_"+ctrl+"/")
	output_Data.to_csv("./Read_count_matrix/"+Exp+"_vs_"+ctrl+"_Union_in_"+fix+"_"+str(res)+"kb_count_mat_CRISPRi_range_res_"+str(res)+"kb.txt",sep="\t")

def main():
	res=10
	combination_set=[(("TA448",'NT'),'ctrl'),(("TA448",'NT'),'stim'),(("stim",'ctrl'),'NT'),(("stim",'ctrl'),'TA448')]
	### define data list
	for var_set_fix in combination_set:
		print(var_set_fix)
		var_set, fix= var_set_fix
		make_matrix_for_target_region(res,var_set,fix)

main()
