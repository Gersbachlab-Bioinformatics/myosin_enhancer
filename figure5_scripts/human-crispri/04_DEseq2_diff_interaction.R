#~/anaconda2/envs/Charlie_capture_HiCAR/bin/R

library( "DESeq2" )
# Deseq2 version 3.17
library(ggplot2)

Run_DESeq2_target<-function(resol){
data=read.table(paste0("./cool2txt/Region3_vs_nontarget_Union_",resol,"kb_count_mat_include_target_region.txt"),header=T)
metadata=read.csv("metadata.csv",header=T,sep=",",row.names="id")
metadata['Group.']<-factor(metadata[,'Group.'])
dds <- DESeqDataSetFromMatrix(countData=data,colData=metadata,design=~Group.)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <- results(dds)
res[is.na(res['padj']),'padj']<-1
write.csv(res, file=paste0("./cool2txt/Region3_vs_nontarget_Union_",resol,"kb_include_target_region_DE_results.csv"))
#res1<-res[which(res$padj<0.1),]
#summary(res1)
#write.csv(res1[rev(order(res1$log2FoldChange)),], file=paste0("./cool2txt/Region3_vs_nontarget_Union_",resol,"kb_include_target_region_DE_results_pval_0.1.csv"))
}


Run_DESeq2<-function(resol){
data=read.table(paste0("./cool2txt/Region3_vs_nontarget_Union_",resol,"kb_count_mat.txt"),header=T)
metadata=read.csv("metadata.csv",header=T,sep=",",row.names="id")
metadata['Group.']<-factor(metadata[,'Group.'])

dds <- DESeqDataSetFromMatrix(countData=data,colData=metadata,design=~type)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <- results(dds)
res[is.na(res['padj']),'padj']<-1
write.csv(res, file=paste0("./cool2txt/Region3_vs_nontarget_Union_",resol,"kb_DE_results.csv"))
#res1<-res[which(res$padj<0.1),]
#summary(res1)
#write.csv(res1[rev(order(res1$log2FoldChange)),], file=paste0("./cool2txt/Region3_vs_nontarget_Union_",resol,"kb_DE_results_pval_0.1.csv"))
}

Run_DESeq2_PABPN1_CTCF<-function(resol){
data=read.table(paste0("./cool2txt/Region3_vs_nontarget_Union_",resol,"kb_count_mat_PABPN1_CTCF_res.txt"),header=T)
metadata=read.csv("metadata.csv",header=T,sep=",",row.names="id")
metadata['Group.']<-factor(metadata[,'Group.'])

dds <- DESeqDataSetFromMatrix(countData=data,colData=metadata,design=~type)
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <- results(dds)
res[is.na(res['padj']),'padj']<-1
write.csv(res, file=paste0("./cool2txt/Region3_vs_nontarget_Union_",resol,"kb_PABPN1_CTCF_res_DE_results.csv"))
#res1<-res[which(res$padj<0.1),]
#summary(res1)
#write.csv(res1[rev(order(res1$log2FoldChange)),], file=paste0("./cool2txt/Region3_vs_nontarget_Union_",resol,"kb_PABPN1_CTCF_res_DE_results_pval_0.1.csv"))
}



Run_DESeq2_target("10")
Run_DESeq2_target("20")
Run_DESeq2("10")
Run_DESeq2("20")
Run_DESeq2_PABPN1_CTCF("20")
Run_DESeq2_PABPN1_CTCF("10")

