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
##################################
setwd("C:/Users/Administrator/OneDrive/10.shangrui/remove_batch_effect")
sample_table = read.table("sample_table_rc.txt",sep="\t",header=T)
rownames(sample_table) = sample_table[,3]

count_matrix = as.matrix(read.table("raw_counts_5_species.txt",header=T, sep="\t"))
colnames(count_matrix) = sample_table[,3]
################################
setwd("C:/Users/Administrator/OneDrive/10.shangrui/coevolution")
sample_table = subset(sample_table,species=="dme"|species=="dsi")
count_matrix = count_matrix[,sample_table[,3]]
st = subset(sample_table,merged_tissue=="AG"|merged_tissue=="OV"|merged_tissue=="TE")
counts=count_matrix[,st[,3]]
write.table(st,"sample_table_dmds.txt",row.names = F,col.names = T,quote=F,sep="\t")
adjusted <- ComBat_seq(counts, batch=st$batch, group=NULL, covar_mod = st[,c("merged_tissue2","species")])

write.table(adjusted,"adjusted_rawCounts_dmds.txt",row.names = T,col.names = T,quote=F,sep="\t")

cal_con_div <- function(counts,st,tissue1,tissue2,file_name){
  #browser()
  st_ym = subset(st,merged_tissue==tissue1|merged_tissue==tissue2)
  st_ym$species = factor(st_ym$species)
  st_ym$species = relevel(st_ym$species,ref="dsi")
  group=paste(st_ym$species,st_ym$merged_tissue,sep=".")
  y <- DGEList(counts=counts[,st_ym[,3]],group = group)
  keep <- filterByExpr(y)
  table(keep)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  
  #write.table(rownames(y),file_name,sep="\t",row.names = F, col.names = F,quote = F)
  
  # tissue evolving: Dyak vs Dmel/Dsim fold change
  expr=cpm(y)
  st12 = subset(st_ym, merged_tissue==tissue1&species=="dsi")
  st13 = subset(st_ym, merged_tissue==tissue2&species=="dsi")
  st2 = subset(st_ym, merged_tissue==tissue1&species!="dsi")
  st3 = subset(st_ym, merged_tissue==tissue2&species!="dsi")
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
  fc[g_div,"type"]=fc[g_div,"type"]+2
  print(table(fc$type))
  fc[fc$type==3,"type"]=0
  fc$type = factor(fc$type,levels=c(0,1,2))
  
  fc$direct=0
  fc$direct=ifelse(fc$type==1,ifelse(fc$fc1>-fc$fc2,"1","2"),fc$direct)
  fc$direct=ifelse(fc$type==2,ifelse(fc$fc1>fc$fc2,"3","4"),fc$direct)
  print(table(fc[,"direct"]))
  
  write.table(fc,file_name,row.names = T,col.names = T,sep = "\t",quote=F)
  return(fc)
}

fc = cal_con_div(adjusted,st,"AG","OV","AGvsOV.ConsDiv.txt")
p<-ggplot()+ geom_abline(intercept = 0, slope = 1, color="grey",linetype="dashed")+ geom_abline(intercept = 0, slope = -1, color="grey",linetype="dashed") + geom_vline(xintercept = 0, color="grey",linetype="dashed")+ geom_hline(yintercept = 0, color="grey",linetype="dashed") + geom_point(data=subset(fc,type==0), aes(x=fc1,y=fc2,color=type),size=0.5)+ geom_point(data=subset(fc,type!=0), aes(x=fc1,y=fc2,color=type),size=0.5)  + xlim(-15,15)+ylim(-15,15) + theme_bw(12) + scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) + xlab("Accessory gland evolving log2(dme/dsi)") + ylab("Ovary evolving log2(dme/dsi)")
ggsave("FC_AGvsOV.pdf", width=4,height = 3)

fc = cal_con_div(adjusted,st,"TE","OV","TEvsOV.ConsDiv.txt")
p<-ggplot()+ geom_abline(intercept = 0, slope = 1, color="grey",linetype="dashed")+ geom_abline(intercept = 0, slope = -1, color="grey",linetype="dashed") + geom_vline(xintercept = 0, color="grey",linetype="dashed")+ geom_hline(yintercept = 0, color="grey",linetype="dashed") + geom_point(data=subset(fc,type==0), aes(x=fc1,y=fc2,color=type),size=0.5)+ geom_point(data=subset(fc,type!=0), aes(x=fc1,y=fc2,color=type),size=0.5)  + xlim(-15,15)+ylim(-15,15) + theme_bw(12) + scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) + xlab("Testis evolving log2(dme/dsi)") + ylab("Ovary evolving log2(dme/dsi)")
ggsave("FC_TEvsOV.pdf", width=4,height = 3)


fc = cal_con_div(adjusted,st,"AG","TE","AGvsTE.ConsDiv.txt")
p<-ggplot()+ geom_abline(intercept = 0, slope = 1, color="grey",linetype="dashed")+ geom_abline(intercept = 0, slope = -1, color="grey",linetype="dashed") + geom_vline(xintercept = 0, color="grey",linetype="dashed")+ geom_hline(yintercept = 0, color="grey",linetype="dashed") + geom_point(data=subset(fc,type==0), aes(x=fc1,y=fc2,color=type),size=0.5)+ geom_point(data=subset(fc,type!=0), aes(x=fc1,y=fc2,color=type),size=0.5)  + xlim(-15,15)+ylim(-15,15) + theme_bw(12) + scale_color_manual(values=c('#999999','#E69F00', '#56B4E9')) + xlab("Accessory gland evolving log2(dme/dsi)") + ylab("Testis evolving log2(dme/dsi)")
ggsave("FC_AGvsTE.pdf", width=4,height = 3)


