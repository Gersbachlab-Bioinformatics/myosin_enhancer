#~/anaconda2/envs/Charlie_capture_HiCAR/bin/python

import pandas as pd
import numpy as np
import math


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

def name_func_filter(input_data,target_region):
	tmp1=input_data[input_data['chrom1']=='chr14']
	tmp2=tmp1[tmp1['chrom2']=='chr14']
	tmp2['region1']=tmp2['chrom1']+"."+tmp2['start1'].astype('str')+"."+tmp2['end1'].astype('str')
	tmp2['region2']=tmp2['chrom2']+"."+tmp2['start2'].astype('str')+"."+tmp2['end2'].astype('str')
	name_vec=tmp2['region1']+"."+tmp2['region2']
	tmp2.index=name_vec
	### filter
	tmp2_re=tmp2.loc[(tmp2["region1"].isin(target_region)) | (tmp2["region2"].isin(target_region))]
	#tmp2_re=tmp2.loc[[x for x in tmp2.index if tmp2.loc[x,'region1']!=tmp2.loc[x,'region2']]]
	return(tmp2_re)


def coordinate_convert(coord,res):
	scale_factor=res*1000
	coord_start, coord_end=coord
	start=math.floor(coord_start/scale_factor)*scale_factor
	end=math.ceil(coord_end/scale_factor)*scale_factor
	target_region=[]
	for i in range(int(start/scale_factor),int(end/scale_factor)):
		target_region.append("chr14."+str(i*scale_factor)+"."+str((i+1)*scale_factor))
	return(target_region)

def main(res):
	nontarget1 = pd.read_csv("./cool2txt/krabnontarget1."+str(res)+"kb.txt.gz", compression='gzip', header=0, sep='\t')
	nontarget2 = pd.read_csv("./cool2txt/krabnontarget2."+str(res)+"kb.txt.gz", compression='gzip', header=0, sep='\t')
	region31 = pd.read_csv("./cool2txt/krabregion31."+str(res)+"kb.txt.gz", compression='gzip', header=0, sep='\t')
	region32 = pd.read_csv("./cool2txt/krabregion32."+str(res)+"kb.txt.gz", compression='gzip', header=0, sep='\t')

	nontarget1_chr14=name_func(nontarget1)
	nontarget2_chr14=name_func(nontarget2)
	region31_chr14=name_func(region31)
	region32_chr14=name_func(region32)

	### get nontarget and region3 union interaction
	nontarget_unique=nontarget1_chr14.index.intersection(nontarget2_chr14.index)
	region3_unique=region31_chr14.index.intersection(region32_chr14.index)
	total_name=nontarget_unique.union(region3_unique)
	
	output_Data=pd.DataFrame(np.zeros((len(total_name),4)),index=total_name,columns=['region31','region32','nontarget1','nontarget2'])
	output_Data.loc[total_name.intersection(region31_chr14.index),'region31']=region31_chr14.loc[total_name.intersection(region31_chr14.index),'count']
	output_Data.loc[total_name.intersection(region32_chr14.index),'region32']=region32_chr14.loc[total_name.intersection(region32_chr14.index),'count']
	####
	output_Data.loc[total_name.intersection(nontarget1_chr14.index),'nontarget1']=nontarget1_chr14.loc[total_name.intersection(nontarget1_chr14.index),'count']
	output_Data.loc[total_name.intersection(nontarget2_chr14.index),'nontarget2']=nontarget2_chr14.loc[total_name.intersection(nontarget2_chr14.index),'count']
	output_Data.to_csv("./cool2txt/Region3_vs_nontarget_Union_"+str(res)+"kb_count_mat.txt",sep="\t")
	#output_Data.to_csv("./cool2txt/Region3_vs_nontarget_Union_"+str(res)+"kb_count_mat_wo_self.txt",sep="\t")

def main_only_target_region(res):
	## define target region
	PABPN1_CTCF=(23318750,23326998)	
	MYH6_gene_body=(23393234,23409273)
	MYH7_promoter=(23434660,23436660)
	Region3_enhancer=(23379863,23380961)
	## conversion	
	PABPN1_CTCF_res=coordinate_convert(PABPN1_CTCF,res)
	MYH6_gene_body_res=coordinate_convert(MYH6_gene_body,res)
	MYH7_promoter_res=coordinate_convert(MYH7_promoter,res)
	Region3_enhancer_res=coordinate_convert(Region3_enhancer,res)
	target_region=PABPN1_CTCF_res+MYH6_gene_body_res+MYH7_promoter_res+Region3_enhancer_res
	## read data
	nontarget1 = pd.read_csv("./cool2txt/krabnontarget1."+str(res)+"kb.txt.gz", compression='gzip', header=0, sep='\t')
	nontarget2 = pd.read_csv("./cool2txt/krabnontarget2."+str(res)+"kb.txt.gz", compression='gzip', header=0, sep='\t')
	region31 = pd.read_csv("./cool2txt/krabregion31."+str(res)+"kb.txt.gz", compression='gzip', header=0, sep='\t')
	region32 = pd.read_csv("./cool2txt/krabregion32."+str(res)+"kb.txt.gz", compression='gzip', header=0, sep='\t')
	
	nontarget1_chr14=name_func_filter(nontarget1,target_region)
	nontarget2_chr14=name_func_filter(nontarget2,target_region)
	region31_chr14=name_func_filter(region31, target_region)
	region32_chr14=name_func_filter(region32, target_region)

	### get nontarget and region3 union interaction
	nontarget_unique=nontarget1_chr14.index.intersection(nontarget2_chr14.index)
	region3_unique=region31_chr14.index.intersection(region32_chr14.index)
	total_name=nontarget_unique.union(region3_unique)

	output_Data=pd.DataFrame(np.zeros((len(total_name),4)),index=total_name,columns=['region31','region32','nontarget1','nontarget2'])
	output_Data.loc[total_name.intersection(region31_chr14.index),'region31']=region31_chr14.loc[total_name.intersection(region31_chr14.index),'count']
	output_Data.loc[total_name.intersection(region32_chr14.index),'region32']=region32_chr14.loc[total_name.intersection(region32_chr14.index),'count']
	####
	output_Data.loc[total_name.intersection(nontarget1_chr14.index),'nontarget1']=nontarget1_chr14.loc[total_name.intersection(nontarget1_chr14.index),'count']
	output_Data.loc[total_name.intersection(nontarget2_chr14.index),'nontarget2']=nontarget2_chr14.loc[total_name.intersection(nontarget2_chr14.index),'count']
	output_Data.to_csv("./cool2txt/Region3_vs_nontarget_Union_"+str(res)+"kb_count_mat_include_target_region.txt",sep="\t")


def main_only_PABPN1_CTCF(res):
	## define target region
	PABPN1_CTCF=(23318750,23326998)	
	## conversion	
	PABPN1_CTCF_res=coordinate_convert(PABPN1_CTCF,res)
	target_region=PABPN1_CTCF_res
	## read data
	nontarget1 = pd.read_csv("./cool2txt/krabnontarget1."+str(res)+"kb.txt.gz", compression='gzip', header=0, sep='\t')
	nontarget2 = pd.read_csv("./cool2txt/krabnontarget2."+str(res)+"kb.txt.gz", compression='gzip', header=0, sep='\t')
	region31 = pd.read_csv("./cool2txt/krabregion31."+str(res)+"kb.txt.gz", compression='gzip', header=0, sep='\t')
	region32 = pd.read_csv("./cool2txt/krabregion32."+str(res)+"kb.txt.gz", compression='gzip', header=0, sep='\t')
	
	nontarget1_chr14=name_func_filter(nontarget1,target_region)
	nontarget2_chr14=name_func_filter(nontarget2,target_region)
	region31_chr14=name_func_filter(region31, target_region)
	region32_chr14=name_func_filter(region32, target_region)

	### get nontarget and region3 union interaction
	nontarget_unique=nontarget1_chr14.index.intersection(nontarget2_chr14.index)
	region3_unique=region31_chr14.index.intersection(region32_chr14.index)
	total_name=nontarget_unique.union(region3_unique)

	output_Data=pd.DataFrame(np.zeros((len(total_name),4)),index=total_name,columns=['region31','region32','nontarget1','nontarget2'])
	output_Data.loc[total_name.intersection(region31_chr14.index),'region31']=region31_chr14.loc[total_name.intersection(region31_chr14.index),'count']
	output_Data.loc[total_name.intersection(region32_chr14.index),'region32']=region32_chr14.loc[total_name.intersection(region32_chr14.index),'count']
	####
	output_Data.loc[total_name.intersection(nontarget1_chr14.index),'nontarget1']=nontarget1_chr14.loc[total_name.intersection(nontarget1_chr14.index),'count']
	output_Data.loc[total_name.intersection(nontarget2_chr14.index),'nontarget2']=nontarget2_chr14.loc[total_name.intersection(nontarget2_chr14.index),'count']
	output_Data.to_csv("./cool2txt/Region3_vs_nontarget_Union_"+str(res)+"kb_count_mat_PABPN1_CTCF_res.txt",sep="\t")
	

res_list=[1,2,5,10,20,40]
for res in res_list:
	#main(res)
	#main_only_target_region(res)
	main_only_PABPN1_CTCF(res)
