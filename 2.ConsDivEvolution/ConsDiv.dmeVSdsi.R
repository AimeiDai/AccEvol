library(sva)
library(pamr)
library(ggplot2)
library(edgeR)
library(reshape)
library(ComplexHeatmap)
options(stringsAsFactors=F)
require(RColorBrewer)
library(circlize)
cols <- brewer.pal(12, "Set1")[1:6]
#setwd("C:/Users/Administrator/OneDrive/10.shangrui/remove_batch_effect")
setwd("/Users/aimeidai/OneDrive/10.shangrui/remove_batch_effect")
sample_table = read.table("sample_table_rc.txt",sep="\t",header=T)
rownames(sample_table) = sample_table[,3]
cpm_adj = read.table("adjusted_raw_count_5species.txt",header=T)
colnames(cpm_adj)=sample_table[,3]

cal_con_div <- function(counts,st,tissue1,tissue2,ref_species,file_name){
 # browser()
  st_ym = subset(st,(merged_tissue==tissue1|merged_tissue==tissue2)&(species=="dme"|species==ref_species))
  st_ym$species = factor(st_ym$species)
  st_ym$species = relevel(st_ym$species,ref=ref_species)
  group=paste(st_ym$species,st_ym$merged_tissue,sep=".")
  y <- DGEList(counts=counts[,st_ym[,3]],group = group)
  keep <- filterByExpr(y)
  table(keep)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  
  #write.table(rownames(y),file_name,sep="\t",row.names = F, col.names = F,quote = F)
  
  # tissue evolving: Dyak vs Dmel/Dsim fold change
  expr=cpm(y)
  st12 = subset(st_ym, merged_tissue==tissue1&species==ref_species)
  st13 = subset(st_ym, merged_tissue==tissue2&species==ref_species)
  st2 = subset(st_ym, merged_tissue==tissue1&species=="dme")
  st3 = subset(st_ym, merged_tissue==tissue2&species=="dme")
  fc1 = log2(apply(expr[,st2$merged_SRR],1,mean)/apply(expr[,st12$merged_SRR],1,mean))
  fc2 = log2(apply(expr[,st3$merged_SRR],1,mean)/apply(expr[,st13$merged_SRR],1,mean))
  fc = data.frame(fc1,fc2)
  
  #concordant evolution
  design = model.matrix(~0+group)
  y <- estimateGLMCommonDisp(y, design)
  fit <- glmQLFit(y, design)
  colnames(fit)
  qlf <- glmQLFTest(fit,contrast = c(0.5,0.5,-0.5,-0.5))
  con = topTags(qlf,nrow(y))$table
  con$ifSig=con$FDR<=0.05
  table(con$ifSig)
  g_con = rownames(con)[con$ifSig]
  
  #divergent evolution
  design = model.matrix(~merged_tissue+merged_tissue:species,st_ym)
  y <- estimateGLMCommonDisp(y, design)
  fit <- glmQLFit(y, design)
  colnames(fit)
  qlf <- glmQLFTest(fit,contrast=c(0,0,1,-1))
  div = topTags(qlf,nrow(y))$table
  div$ifSig=div$FDR<=0.05
  table(div$ifSig)
  g_div = rownames(div)[div$ifSig]
  
  fc$type=0
  fc[g_con,"type"]=fc[g_con,"type"]+1
  fc[g_div,"type"]=fc[g_div,"type"]-1
  print(table(fc$type))
  #fc[fc$type==3,"type"]=0
  fc$type = factor(fc$type,levels=c(1,0,-1))
  
  fc$direct=0
  fc$direct=ifelse(fc$type==1,ifelse(fc$fc1>0&fc$fc2>0,"1","2"),fc$direct)
  fc$direct=ifelse(fc$type==-1,ifelse(fc$fc1<0&fc$fc2>0,"3","4"),fc$direct)
  print(table(fc[,"direct"]))
  
  write.table(fc,paste(ref_species,file_name,sep=""),row.names = T,col.names = T,sep = "\t",quote=F)
  return(fc)
}
#setwd("C:/Users/Administrator/OneDrive/10.shangrui/ConsDivEvolution")
setwd("/Users/aimeidai/OneDrive/10.shangrui/ConsDivEvolution")
######################################
ref_species="dsi"
fc = cal_con_div(cpm_adj,sample_table,"AG","TE",ref_species,"_AGvsTE.ConsDiv.txt")
p<-ggplot()+ geom_abline(intercept = 0, slope = 1, color="grey",linetype="dashed")+ geom_abline(intercept = 0, slope = -1, color="grey",linetype="dashed") + geom_vline(xintercept = 0, color="grey",linetype="dashed")+ geom_hline(yintercept = 0, color="grey",linetype="dashed") + geom_point(data=subset(fc,type==0), aes(x=fc1,y=fc2,color=type),size=0.6)+ geom_point(data=subset(fc,type!=0), aes(x=fc1,y=fc2,color=type),size=0.6)  + xlim(-5,5)+ylim(-5,5) + theme_bw(10) + scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) + ylab(paste("Testis evolving log2(dme/",ref_species,")",sep="")) +  xlab(paste("Accessory gland evolving log2(dme/",ref_species,")",sep=""))
ggsave(paste("FC_AGvsTE_",ref_species,".pdf",sep=""), width=4.5,height = 3.5)

fc = cal_con_div(cpm_adj,sample_table,"TE","OV",ref_species,"_TEvsOV.ConsDiv.txt")
p<-ggplot()+ geom_abline(intercept = 0, slope = 1, color="grey",linetype="dashed")+ geom_abline(intercept = 0, slope = -1, color="grey",linetype="dashed") + geom_vline(xintercept = 0, color="grey",linetype="dashed")+ geom_hline(yintercept = 0, color="grey",linetype="dashed") + geom_point(data=subset(fc,type==0), aes(x=fc1,y=fc2,color=type),size=0.6)+ geom_point(data=subset(fc,type!=0), aes(x=fc1,y=fc2,color=type),size=0.6)  + xlim(-5,5)+ylim(-5,5) + theme_bw(10) + scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))  + xlab(paste("Testis evolving log2(dme/",ref_species,")",sep="")) + ylab(paste("Ovary evolving log2(dme/",ref_species,")",sep=""))
ggsave(paste("FC_TEvsOV_",ref_species,".pdf",sep=""),  width=4.5,height = 3.5)


fc = cal_con_div(cpm_adj,sample_table,"AG","OV",ref_species,"_AGvsOV.ConsDiv.txt")
p<-ggplot()+ geom_abline(intercept = 0, slope = 1, color="grey",linetype="dashed")+ geom_abline(intercept = 0, slope = -1, color="grey",linetype="dashed") + geom_vline(xintercept = 0, color="grey",linetype="dashed")+ geom_hline(yintercept = 0, color="grey",linetype="dashed") + geom_point(data=subset(fc,type==0), aes(x=fc1,y=fc2,color=type),size=0.6)+ geom_point(data=subset(fc,type!=0), aes(x=fc1,y=fc2,color=type),size=0.6)  + xlim(-5,5)+ylim(-5,5) + theme_bw(10) + scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))  + xlab(paste("Accessory gland evolving log2(dme/",ref_species,")",sep="")) + ylab(paste("Ovary evolving log2(dme/",ref_species,")",sep=""))
ggsave(paste("FC_AGvsOV_",ref_species,".pdf",sep=""),  width=4.5,height = 3.5)

############################################
##############################################
ref_species="dan"
fc = cal_con_div(cpm_adj,sample_table,"AG","TE",ref_species,"_AGvsTE.ConsDiv.txt")
p<-ggplot()+ geom_abline(intercept = 0, slope = 1, color="grey",linetype="dashed")+ geom_abline(intercept = 0, slope = -1, color="grey",linetype="dashed") + geom_vline(xintercept = 0, color="grey",linetype="dashed")+ geom_hline(yintercept = 0, color="grey",linetype="dashed") + geom_point(data=subset(fc,type==0), aes(x=fc1,y=fc2,color=type),size=0.6)+ geom_point(data=subset(fc,type!=0), aes(x=fc1,y=fc2,color=type),size=0.6)  + xlim(-5,5)+ylim(-5,5) + theme_bw(10) + scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))  + ylab(paste("Testis evolving log2(dme/",ref_species,")",sep="")) + xlab(paste("Accessory gland evolving log2(dme/",ref_species,")",sep=""))
ggsave(paste("FC_AGvsTE_",ref_species,".pdf",sep=""),  width=4.5,height = 3.5)

fc = cal_con_div(cpm_adj,sample_table,"TE","OV",ref_species,"_TEvsOV.ConsDiv.txt")
p<-ggplot()+ geom_abline(intercept = 0, slope = 1, color="grey",linetype="dashed")+ geom_abline(intercept = 0, slope = -1, color="grey",linetype="dashed") + geom_vline(xintercept = 0, color="grey",linetype="dashed")+ geom_hline(yintercept = 0, color="grey",linetype="dashed") + geom_point(data=subset(fc,type==0), aes(x=fc1,y=fc2,color=type),size=0.6)+ geom_point(data=subset(fc,type!=0), aes(x=fc1,y=fc2,color=type),size=0.6)  + xlim(-5,5)+ylim(-5,5) + theme_bw(10) + scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))  + xlab(paste("Testis evolving log2(dme/",ref_species,")",sep="")) + ylab(paste("Ovary evolving log2(dme/",ref_species,")",sep=""))
ggsave(paste("FC_TEvsOV_",ref_species,".pdf",sep=""),  width=4.5,height = 3.5)


fc = cal_con_div(cpm_adj,sample_table,"AG","OV",ref_species,"_AGvsOV.ConsDiv.txt")
p<-ggplot()+ geom_abline(intercept = 0, slope = 1, color="grey",linetype="dashed")+ geom_abline(intercept = 0, slope = -1, color="grey",linetype="dashed") + geom_vline(xintercept = 0, color="grey",linetype="dashed")+ geom_hline(yintercept = 0, color="grey",linetype="dashed") + geom_point(data=subset(fc,type==0), aes(x=fc1,y=fc2,color=type),size=0.6)+ geom_point(data=subset(fc,type!=0), aes(x=fc1,y=fc2,color=type),size=0.6)  + xlim(-5,5)+ylim(-5,5) + theme_bw(10) + scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))  + xlab(paste("Accessory gland evolving log2(dme/",ref_species,")",sep="")) + ylab(paste("Ovary evolving log2(dme/",ref_species,")",sep=""))
ggsave(paste("FC_AGvsOV_",ref_species,".pdf",sep=""),  width=4.5,height = 3.5)

########################################
ref_species="dya"
fc = cal_con_div(cpm_adj,sample_table,"AG","TE",ref_species,"_AGvsTE.ConsDiv.txt")
p<-ggplot()+ geom_abline(intercept = 0, slope = 1, color="grey",linetype="dashed")+ geom_abline(intercept = 0, slope = -1, color="grey",linetype="dashed") + geom_vline(xintercept = 0, color="grey",linetype="dashed")+ geom_hline(yintercept = 0, color="grey",linetype="dashed") + geom_point(data=subset(fc,type==0), aes(x=fc1,y=fc2,color=type),size=0.6)+ geom_point(data=subset(fc,type!=0), aes(x=fc1,y=fc2,color=type),size=0.6)  + xlim(-5,5)+ylim(-5,5) + theme_bw(10) + scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) + ylab(paste("Testis evolving log2(dme/",ref_species,")",sep="")) + xlab(paste("Accessory gland evolving log2(dme/",ref_species,")",sep=""))
ggsave(paste("FC_AGvsTE_",ref_species,".pdf",sep=""),  width=4.5,height = 3.5)

fc = cal_con_div(cpm_adj,sample_table,"TE","OV",ref_species,"_TEvsOV.ConsDiv.txt")
p<-ggplot()+ geom_abline(intercept = 0, slope = 1, color="grey",linetype="dashed")+ geom_abline(intercept = 0, slope = -1, color="grey",linetype="dashed") + geom_vline(xintercept = 0, color="grey",linetype="dashed")+ geom_hline(yintercept = 0, color="grey",linetype="dashed") + geom_point(data=subset(fc,type==0), aes(x=fc1,y=fc2,color=type),size=0.6)+ geom_point(data=subset(fc,type!=0), aes(x=fc1,y=fc2,color=type),size=0.6)  + xlim(-5,5)+ylim(-5,5) + theme_bw(10) + scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) + xlab(paste("Testis evolving log2(dme/",ref_species,")",sep="")) + ylab(paste("Ovary evolving log2(dme/",ref_species,")",sep=""))
ggsave(paste("FC_TEvsOV_",ref_species,".pdf",sep=""),  width=4.5,height = 3.5)


fc = cal_con_div(cpm_adj,sample_table,"AG","OV",ref_species,"_AGvsOV.ConsDiv.txt")
p<-ggplot()+ geom_abline(intercept = 0, slope = 1, color="grey",linetype="dashed")+ geom_abline(intercept = 0, slope = -1, color="grey",linetype="dashed") + geom_vline(xintercept = 0, color="grey",linetype="dashed")+ geom_hline(yintercept = 0, color="grey",linetype="dashed") + geom_point(data=subset(fc,type==0), aes(x=fc1,y=fc2,color=type),size=0.6)+ geom_point(data=subset(fc,type!=0), aes(x=fc1,y=fc2,color=type),size=0.6)  + xlim(-5,5)+ylim(-5,5) + theme_bw(10) + scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) + xlab(paste("Accessory gland evolving log2(dme/",ref_species,")",sep="")) + ylab(paste("Ovary evolving log2(dme/",ref_species,")",sep=""))
ggsave(paste("FC_AGvsOV_",ref_species,".pdf",sep=""),  width=4.5,height = 3.5)

########################################
ref_species="dps"
fc = cal_con_div(cpm_adj,sample_table,"AG","TE",ref_species,"_AGvsTE.ConsDiv.txt")
p<-ggplot()+ geom_abline(intercept = 0, slope = 1, color="grey",linetype="dashed")+ geom_abline(intercept = 0, slope = -1, color="grey",linetype="dashed") + geom_vline(xintercept = 0, color="grey",linetype="dashed")+ geom_hline(yintercept = 0, color="grey",linetype="dashed") + geom_point(data=subset(fc,type==0), aes(x=fc1,y=fc2,color=type),size=0.6)+ geom_point(data=subset(fc,type!=0), aes(x=fc1,y=fc2,color=type),size=0.6)  + xlim(-10,10)+ylim(-10,10) + theme_bw(10) + scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))  + ylab(paste("Testis evolving log2(dme/",ref_species,")",sep="")) + xlab(paste("Accessory gland evolving log2(dme/",ref_species,")",sep=""))
ggsave(paste("FC_AGvsTE_",ref_species,".pdf",sep=""),  width=3.8,height = 3)

fc = cal_con_div(cpm_adj,sample_table,"TE","OV",ref_species,"_TEvsOV.ConsDiv.txt")
p<-ggplot()+ geom_abline(intercept = 0, slope = 1, color="grey",linetype="dashed")+ geom_abline(intercept = 0, slope = -1, color="grey",linetype="dashed") + geom_vline(xintercept = 0, color="grey",linetype="dashed")+ geom_hline(yintercept = 0, color="grey",linetype="dashed") + geom_point(data=subset(fc,type==0), aes(x=fc1,y=fc2,color=type),size=0.6)+ geom_point(data=subset(fc,type!=0), aes(x=fc1,y=fc2,color=type),size=0.6)  + xlim(-10,10)+ylim(-10,10) + theme_bw(10) + scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))  + xlab(paste("Testis evolving log2(dme/",ref_species,")",sep="")) + ylab(paste("Ovary evolving log2(dme/",ref_species,")",sep=""))
ggsave(paste("FC_TEvsOV_",ref_species,".pdf",sep=""),  width=3.8,height = 3)


fc = cal_con_div(cpm_adj,sample_table,"AG","OV",ref_species,"_AGvsOV.ConsDiv.txt")
p<-ggplot()+ geom_abline(intercept = 0, slope = 1, color="grey",linetype="dashed")+ geom_abline(intercept = 0, slope = -1, color="grey",linetype="dashed") + geom_vline(xintercept = 0, color="grey",linetype="dashed")+ geom_hline(yintercept = 0, color="grey",linetype="dashed") + geom_point(data=subset(fc,type==0), aes(x=fc1,y=fc2,color=type),size=0.6)+ geom_point(data=subset(fc,type!=0), aes(x=fc1,y=fc2,color=type),size=0.6)  + xlim(-10,10)+ylim(-10,10) + theme_bw(10) + scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) + xlab(paste("Accessory gland evolving log2(dme/",ref_species,")",sep="")) + ylab(paste("Ovary evolving log2(dme/",ref_species,")",sep=""))
ggsave(paste("FC_AGvsOV_",ref_species,".pdf",sep=""),  width=3.8,height = 3)




