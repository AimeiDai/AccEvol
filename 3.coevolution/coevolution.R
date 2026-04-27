library(ggplot2)
library(reshape)
library(ComplexHeatmap)
options(stringsAsFactors=F)
require(RColorBrewer)
library(circlize)
col_fun = colorRamp2(c(0, 1.5, 3), c("blue", "white", "red"))
set.seed(1)

#setwd("C:/Users/Administrator/OneDrive/10.shangrui/remove_batch_effect")
setwd("/Users/aimeidai/OneDrive/10.shangrui/remove_batch_effect")
sample_table = read.table("sample_table_rc.txt",sep="\t",header=T)
rownames(sample_table) = sample_table[,3]
cpm_adj = read.table("CPM_adjusted_5species.txt",header=T)
colnames(cpm_adj)=sample_table[,3]

st=sample_table
tpm_s = cpm_adj[,st[,3]]
tpm_m = do.call("cbind",tapply(st[,3],st[,c(4,13)],function(x){
  #browser()
  if(length(x)>1){
    y = as.matrix(apply(tpm_s[,x],1,mean))
  }else{
    y = as.matrix(tpm_s[,x])
  }
  colnames(y)=paste(unique(st[st[,3]%in%x,c(4,13)]),collapse =".");
  return(y)
}))


setwd("../coevolution")

permutation_test <- function(m,n,M,N,real,times=1000){
  x=rep(0,times)
  for(i in 1:times){
    x1 = sample(M, m)
    x2 = sample(N, n)
    x[i] = length(intersect(x1,x2))
  }
  #browser()
  return(length(x[x>=real])/times)
}

########## clustering along species
########## AG
mat_ag = tpm_m[,c("dme.AG","dsi.AG","dya.AG","dan.AG","dps.AG")]
colnames(mat_ag) = substr(colnames(mat_ag),1,3)
pdf("AG.Heatmap.pdf",height = 5,width=3)
ht_ag=draw(Heatmap(log10(mat_ag+1),name = "log10(TPM)", km =6, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=F,col = col_fun))
dev.off()
#, row_title_gp = gpar(fill = cols)
rcl.ag <- row_order(ht_ag)

###########OV
mat_ov = tpm_m[,c("dme.OV","dsi.OV","dya.OV","dan.OV","dps.OV")]
colnames(mat_ov) = substr(colnames(mat_ov),1,3)
pdf("OV.Heatmap.pdf",height = 5,width=3)
ht_ov=draw(Heatmap(log10(mat_ov+1),name = "log10(TPM)", km =6, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=F,col = col_fun))
dev.off()
rcl.ov <- row_order(ht_ov)
#### TE
mat_te = tpm_m[,c("dme.TE","dsi.TE","dya.TE","dan.TE","dps.TE")]
colnames(mat_te) = substr(colnames(mat_te),1,3)
pdf("TE.Heatmap.pdf",height = 5,width=3)
ht_te=draw(Heatmap(log10(mat_te+1),name = "log10(TPM)", km =6, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=F,col = col_fun))
dev.off()
rcl.te <- row_order(ht_te)
################ write clusters into files
te_cl = melt(lapply(rcl.te,function(x){
  rownames(mat_te)[x]}))
write.table(te_cl,"TE_cluster.txt",row.names = F,col.names = F,sep="\t",quote=F)
ov_cl = melt(lapply(rcl.ov,function(x){
  rownames(mat_ov)[x]}))
write.table(ov_cl,"OV_cluster.txt",row.names = F,col.names = F,sep="\t",quote=F)
ag_cl = melt(lapply(rcl.ag,function(x){
  rownames(mat_ag)[x]}))
write.table(ag_cl,"AG_cluster.txt",row.names = F,col.names = F,sep="\t",quote=F)

########## foldchanges, dme/dps
fc_ag = log2(tpm_m[,"dme.AG"]/tpm_m[,"dps.AG"])
fc_te = log2(tpm_m[,"dme.TE"]/tpm_m[,"dps.TE"])
fc_ov = log2(tpm_m[,"dme.OV"]/tpm_m[,"dps.OV"])
fc = cbind(AG=fc_ag,TE=fc_te,OV=fc_ov)
write.table(fc,"logFC_dmVSdp.txt",row.names = T,col.names = T,quote=F,sep="\t")
pdf_name="clusterFCdmVSdp.pdf"

#########################
ag_gg = melt(lapply(rcl.ag,function(x){
  fc_ag[rownames(mat_ag)[x]]}))
ag_gg=rbind(ag_gg,data.frame(value=fc_ag,L1="0"))

ov_gg = melt(lapply(rcl.ov,function(x){
  fc_ov[rownames(mat_ov)[x]]}))
ov_gg=rbind(ov_gg,data.frame(value=fc_ov,L1="0"))

te_gg = melt(lapply(rcl.te,function(x){
  fc_te[rownames(mat_te)[x]]}))
te_gg=rbind(te_gg,data.frame(value=fc_te,L1="0"))

fc_all = rbind(data.frame(ag_gg,tissue="AG"),data.frame(te_gg,tissue="TE"))
fc_all = rbind(fc_all,data.frame(ov_gg,tissue="OV"))
 
###############
fc_all = subset(fc_all,!is.na(value)&abs(value)!=Inf)
########################

sig=by(fc_all[,1:2],fc_all[,3],function(x){
  #browser()
  sig=data.frame()
  x0 = subset(x,L1=="0")[,1]
  for(i in 1:6){
    x1 = subset(x,L1==i)[,1]
    p = wilcox.test(x0,x1)$p.value
    s=ifelse(p<=0.001,"***",ifelse(p<=0.01,"**",ifelse(p<0.05,"*","ns")))
    d=median(x1,na.rm=T)-median(x0,na.rm=T)
    sig=rbind(sig,c(s,d))
  }
  colnames(sig)=c("signif","direct")
  return(sig)
})
sig_gg = rbind(data.frame(sig$TE,tissue="TE",L1=rownames(sig$TE)),data.frame(sig$OV,tissue="OV",L1=rownames(sig$OV)))
sig_gg = rbind(sig_gg,data.frame(sig$AG,tissue="AG",L1=rownames(sig$AG)))

########################
y = data.frame(tissue="TE",melt(lapply(rcl.te, length)))
y = rbind(y,data.frame(tissue="OV",melt(lapply(rcl.ov, length))))
y = rbind(y,data.frame(tissue="AG",melt(lapply(rcl.ag, length))))
y$clusterSize = y$value/max(y$value)*6

y_sig = merge(y,sig_gg,by.x=c(1,3),by.y=c(3,4))
y_sig$weight = ifelse(y_sig$signif=="ns",0,1)
y_sig$color = ifelse(y_sig$signif == "ns","0",ifelse(y_sig$direct>0,"1","-1"))

fc_all_color = merge(fc_all,y_sig[,c(1,2,8)],by.x=c(2,3),by.y=c(2,1),all.x=T)
fc_all_color[is.na(fc_all_color$color),"color"]="background"
cols <- c(brewer.pal(9, "Paired")[c(2,6)],"darkgrey","black")
names(cols)=c("-1","1","0","background")

data_summary <- function(x) {
  m <- boxplot(x,plot=F)$stats
  return(c(y=m[3],ymin=m[2],ymax=m[4]))
}
ds_median = subset(melt(tapply(fc_all_color[,3],fc_all_color[,1:2],function(x){boxplot(x,plot=F)$stats[3]})),L1==0)
dp_fc = melt(tapply(fc_all_color[,3],fc_all_color[,1:2],function(x){boxplot(x,plot=F)$stats[3]}))
cols_other = brewer.pal(9, "Greys")[c(3,5,7,9)]
p <- ggplot() + geom_hline(data=ds_median,aes(yintercept = value),color="darkgrey",size=0.2) + geom_line(data=dp_fc,aes(x=as.factor(L1),y=value),group = 1,color=cols_other[4],linetype="dashed")+ stat_summary(data=fc_all_color,aes(x=as.factor(L1), y= value,color = color), fun.data=data_summary,size=1.2) + facet_wrap(tissue~.,scales="free") + theme_bw(12) + theme(axis.text.x=element_text(angle=0,hjust=0.5),strip.background = element_blank(), axis.title.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + xlab("Clusters") + ylab("log2FC (dme/outgroup)")+ coord_cartesian(ylim=c(-4.5,3)) +scale_color_manual(values=cols)
#ggsave(pdf_name,width = 7.2, height = 2.6)

######## add lines of FCs using other outgroups
######### foldchanges, dme/dsi
fc_ag = log2(tpm_m[,"dme.AG"]/tpm_m[,"dsi.AG"])
fc_te = log2(tpm_m[,"dme.TE"]/tpm_m[,"dsi.TE"])
fc_ov = log2(tpm_m[,"dme.OV"]/tpm_m[,"dsi.OV"])

ag_gg = melt(lapply(rcl.ag,function(x){
  fc_ag[rownames(mat_ag)[x]]}))
ag_gg=rbind(ag_gg,data.frame(value=fc_ag,L1="0"))

ov_gg = melt(lapply(rcl.ov,function(x){
  fc_ov[rownames(mat_ov)[x]]}))
ov_gg=rbind(ov_gg,data.frame(value=fc_ov,L1="0"))

te_gg = melt(lapply(rcl.te,function(x){
  fc_te[rownames(mat_te)[x]]}))
te_gg=rbind(te_gg,data.frame(value=fc_te,L1="0"))

fc_all = rbind(data.frame(ag_gg,tissue="AG"),data.frame(te_gg,tissue="TE"),data.frame(ov_gg,tissue="OV"))
ds_fc = melt(tapply(fc_all[,1],fc_all[,2:3],function(x){boxplot(x,plot=F)$stats[3]}))

p_ds <- p+geom_point(data=ds_fc,aes(x=as.factor(L1),y=value),color=cols_other[1],size=0.8) + facet_wrap(tissue~.) + geom_line(data=ds_fc,aes(x=as.factor(L1),y=value),group = 1,color=cols_other[1],linetype="dashed")
######### foldchanges, dme/dsi
fc_ag = log2(tpm_m[,"dme.AG"]/tpm_m[,"dya.AG"])
fc_te = log2(tpm_m[,"dme.TE"]/tpm_m[,"dya.TE"])
fc_ov = log2(tpm_m[,"dme.OV"]/tpm_m[,"dya.OV"])

ag_gg = melt(lapply(rcl.ag,function(x){
  fc_ag[rownames(mat_ag)[x]]}))
ag_gg=rbind(ag_gg,data.frame(value=fc_ag,L1="0"))

ov_gg = melt(lapply(rcl.ov,function(x){
  fc_ov[rownames(mat_ov)[x]]}))
ov_gg=rbind(ov_gg,data.frame(value=fc_ov,L1="0"))

te_gg = melt(lapply(rcl.te,function(x){
  fc_te[rownames(mat_te)[x]]}))
te_gg=rbind(te_gg,data.frame(value=fc_te,L1="0"))

fc_all = rbind(data.frame(ag_gg,tissue="AG"),data.frame(te_gg,tissue="TE"),data.frame(ov_gg,tissue="OV"))
dy_fc = melt(tapply(fc_all[,1],fc_all[,2:3],function(x){boxplot(x,plot=F)$stats[3]}))
p_dy <- p_ds+geom_point(data=dy_fc,aes(x=as.factor(L1),y=value),color=cols_other[2],size=0.8) + facet_wrap(tissue~.) + geom_line(data=dy_fc,aes(x=as.factor(L1),y=value),group = 1,color=cols_other[2],linetype="dashed")

######### foldchanges, dme/dan
fc_ag = log2(tpm_m[,"dme.AG"]/tpm_m[,"dan.AG"])
fc_te = log2(tpm_m[,"dme.TE"]/tpm_m[,"dan.TE"])
fc_ov = log2(tpm_m[,"dme.OV"]/tpm_m[,"dan.OV"])

ag_gg = melt(lapply(rcl.ag,function(x){
  fc_ag[rownames(mat_ag)[x]]}))
ag_gg=rbind(ag_gg,data.frame(value=fc_ag,L1="0"))

ov_gg = melt(lapply(rcl.ov,function(x){
  fc_ov[rownames(mat_ov)[x]]}))
ov_gg=rbind(ov_gg,data.frame(value=fc_ov,L1="0"))

te_gg = melt(lapply(rcl.te,function(x){
  fc_te[rownames(mat_te)[x]]}))
te_gg=rbind(te_gg,data.frame(value=fc_te,L1="0"))

fc_all = rbind(data.frame(ag_gg,tissue="AG"),data.frame(te_gg,tissue="TE"),data.frame(ov_gg,tissue="OV"))
da_fc = melt(tapply(fc_all[,1],fc_all[,2:3],function(x){boxplot(x,plot=F)$stats[3]}))
p_da <- p_dy+geom_point(data=da_fc,aes(x=as.factor(L1),y=value),color=cols_other[3],size=0.8) + facet_wrap(tissue~.) + geom_line(data=da_fc,aes(x=as.factor(L1),y=value),group = 1,color=cols_other[3],linetype="dashed")
ggsave("clusterFC_all2.pdf",width = 8.2, height = 2.8)


##########################number of DEGs Dmel/Dpse
library(edgeR)
cal_con_div <- function(counts,st,species1,species2,tissue1){
  #browser()
  st_ym = subset(st,merged_tissue==tissue1&(species==species1|species==species2))
  st_ym$species = factor(st_ym$species)
  st_ym$species = relevel(st_ym$species,ref=species2)
  group=paste(st_ym$species,st_ym$merged_tissue,sep=".")
  y <- DGEList(counts=counts[,st_ym[,3]],group = group)
  keep <- filterByExpr(y)
  table(keep)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  expr=cpm(y)
  #concordant evolution
  design = model.matrix(~0+group)
  y <- estimateGLMCommonDisp(y, design)
  fit <- glmQLFit(y, design)
  colnames(fit)
  qlf <- glmQLFTest(fit,contrast = c(1,-1))
  con = topTags(qlf,nrow(y))$table
  con$ifSig=ifelse(con$FDR>0.05,"N.S.",ifelse(con$logFC>0,"Up","Down"))
  table(con$ifSig)
  g_con = rownames(con)[con$ifSig]
  return(con)
}

sample_table = read.table("../remove_batch_effect/sample_table_rc.txt",sep="\t",header=T)
rownames(sample_table) = sample_table[,3]

count_matrix = as.matrix(read.table("../remove_batch_effect/adjusted_raw_count_5species.txt",header=T, sep="\t"))
colnames(count_matrix) = sample_table[,3]

deg_ag=cal_con_div(count_matrix,sample_table,"dme","dps","AG")
write.table(deg_ag,"DEG_AG_dmdp.txt",row.names = T, col.names = F,sep="\t",quote=F)
ag_cl$ifSig = deg_ag[ag_cl[,1],"ifSig"]
ag_cl[is.na(ag_cl$ifSig),"ifSig"] = "Filtered"
x = table(ag_cl[,c(2,3)])
x = rbind("0"=colSums(x),x)
for(i in 2:nrow(x)){
  m = rbind(x[1,c(1,3)],x[i,c(1,3)])
  p=fisher.test(m)$p.value
  print(ifelse(p<0.005,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","ns"))))
}
x_gg = data.frame(melt(x/rowSums(x)),tissue="AG")

deg_ov=cal_con_div(count_matrix,sample_table,"dme","dps","OV")
write.table(deg_ov,"DEG_OV_dmdp.txt",row.names = T, col.names = F,sep="\t",quote=F)

ov_cl$ifSig = deg_ov[ov_cl[,1],"ifSig"]
ov_cl[is.na(ov_cl$ifSig),"ifSig"] = "Filtered"
x = table(ov_cl[,c(2,3)])
x = rbind("0"=colSums(x),x)
for(i in 2:nrow(x)){
  m = rbind(x[1,c(1,4)],x[i,c(1,4)])
  p=fisher.test(m)$p.value
  print(ifelse(p<0.005,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","ns"))))
}
x_gg = rbind(x_gg,data.frame(melt(x/rowSums(x)),tissue="OV"))

deg_te=cal_con_div(count_matrix,sample_table,"dme","dps","TE")
write.table(deg_te,"DEG_TE_dmdp.txt",row.names = T, col.names = F,sep="\t",quote=F)

te_cl$ifSig = deg_te[te_cl[,1],"ifSig"]
te_cl[is.na(te_cl$ifSig),"ifSig"] = "Filtered"
x = table(te_cl[,c(2,3)])
x = rbind("0"=colSums(x),x)
for(i in 2:nrow(x)){
  m = rbind(x[1,c(1,4)],x[i,c(1,4)])
  p=fisher.test(m)$p.value
  print(ifelse(p<0.005,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","ns"))))
}
x_gg = rbind(x_gg,data.frame(melt(x/rowSums(x)),tissue="TE"))
x_gg$X2 = factor(x_gg$X2,levels=c("Down","N.S.","Filtered","Up"))
cols <- c(brewer.pal(9, "Paired")[c(1,5)],"lightgrey",brewer.pal(9, "Pastel1")[9])
names(cols)=c("Down","Up","N.S.","Filtered")
p <- ggplot() + geom_bar(data=x_gg,aes(x=as.character(X1),y=value,fill = X2),stat="identity",color="black",linewidth=0.01)+ facet_wrap(tissue~.,scales="free") + theme_bw(12) + theme(axis.text.x=element_text(angle=0,hjust=0.5)) + xlab("Clusters") + ylab("Fraction of DEGs (dme/dps)")+scale_fill_manual(values=cols)
ggsave("Proportion_of_DEGs_by_Cluster_dmdp.pdf",width = 7, height = 2.6)

######## dme vs dsi
deg_ag=cal_con_div(count_matrix,sample_table,"dme","dsi","AG")
ag_cl$ifSig = deg_ag[ag_cl[,1],"ifSig"]
ag_cl[is.na(ag_cl$ifSig),"ifSig"] = "Filtered"
x = table(ag_cl[,c(2,3)])
x = rbind("0"=colSums(x),x)
for(i in 2:nrow(x)){
  m = rbind(x[1,c(1,3)],x[i,c(1,3)])
  p=fisher.test(m)$p.value
  print(ifelse(p<0.005,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","ns"))))
}
x_gg = data.frame(melt(x/rowSums(x)),tissue="AG")

deg_ov=cal_con_div(count_matrix,sample_table,"dme","dsi","OV")
ov_cl$ifSig = deg_ov[ov_cl[,1],"ifSig"]
ov_cl[is.na(ov_cl$ifSig),"ifSig"] = "Filtered"
x = table(ov_cl[,c(2,3)])
x = rbind("0"=colSums(x),x)
for(i in 2:nrow(x)){
  m = rbind(x[1,c(1,4)],x[i,c(1,4)])
  p=fisher.test(m)$p.value
  print(ifelse(p<0.005,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","ns"))))
}
x_gg = rbind(x_gg,data.frame(melt(x/rowSums(x)),tissue="OV"))

deg_te=cal_con_div(count_matrix,sample_table,"dme","dsi","TE")
te_cl$ifSig = deg_te[te_cl[,1],"ifSig"]
te_cl[is.na(te_cl$ifSig),"ifSig"] = "Filtered"
x = table(te_cl[,c(2,3)])
x = rbind("0"=colSums(x),x)
for(i in 2:nrow(x)){
  m = rbind(x[1,c(1,4)],x[i,c(1,4)])
  p=fisher.test(m)$p.value
  print(ifelse(p<0.005,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","ns"))))
}
x_gg = rbind(x_gg,data.frame(melt(x/rowSums(x)),tissue="TE"))
x_gg$X2 = factor(x_gg$X2,levels=c("Down","N.S.","Filtered","Up"))
cols <- c(brewer.pal(9, "Paired")[c(1,5)],"lightgrey",brewer.pal(9, "Pastel1")[9])
names(cols)=c("Down","Up","N.S.","Filtered")
p <- ggplot() + geom_bar(data=x_gg,aes(x=as.character(X1),y=value,fill = X2),stat="identity",color="black",linewidth=0.01)+ facet_wrap(tissue~.,scales="free") + theme_bw(12) + theme(axis.text.x=element_text(angle=0,hjust=0.5)) + xlab("Clusters") + ylab("Fraction of DEGs (dme/dsi)")+scale_fill_manual(values=cols)
ggsave("Proportion_of_DEGs_by_Cluster_dmds.pdf",width = 7, height = 3)



#####
ov_ag_c=matrix(nrow=6,ncol=6)
ov_ag_p=matrix(nrow=6,ncol=6)
for (i in 1:6){
  for (j in 1:6){
    gs1 = rownames(mat_ov)[rcl.ov[[i]]]
    gs2 = rownames(mat_ag)[rcl.ag[[j]]]
    overlap = length(intersect(gs1,gs2))
    ov_ag_c[i,j]=overlap/min(length(gs1),length(gs2))
    ov_ag_p[i,j]=permutation_test(length(gs1),length(gs2),nrow(tpm_m),nrow(tpm_m),overlap)
  }
}

#####
te_ag_c=matrix(nrow=6,ncol=6)
te_ag_p=matrix(nrow=6,ncol=6)
for (i in 1:6){
  for (j in 1:6){
    gs1 = rcl.te[[i]]
    gs2 = rcl.ag[[j]]
    overlap = length(intersect(gs1,gs2))
    te_ag_c[i,j]=overlap/min(length(gs1),length(gs2))
    te_ag_p[i,j]=permutation_test(length(gs1),length(gs2),nrow(tpm_m),nrow(tpm_m),overlap)
  }
}

#####
ov_te_c=matrix(nrow=6,ncol=6)
ov_te_p=matrix(nrow=6,ncol=6)
for (i in 1:6){
  for (j in 1:6){
    gs1 = rcl.ov[[i]]
    gs2 = rcl.te[[j]]
    overlap = length(intersect(gs1,gs2))
    ov_te_c[i,j]=overlap/min(length(gs1),length(gs2))
    ov_te_p[i,j]=permutation_test(length(gs1),length(gs2),nrow(tpm_m),nrow(tpm_m),overlap)
  }
}

#####################
x = data.frame(sec1="TE",sec2="AG",merge(melt(te_ag_c),melt(te_ag_p),by.x=c(1,2),by.y=c(1,2)))
x = rbind(x,data.frame(sec1="OV",sec2="AG",merge(melt(ov_ag_c),melt(ov_ag_p),by.x=c(1,2),by.y=c(1,2))))
x = rbind(x,data.frame(sec1="OV",sec2="TE",merge(melt(ov_te_c),melt(ov_te_p),by.x=c(1,2),by.y=c(1,2))))
# x$direct=apply(x,1,function(a){
#   #browser()
#   a1=y_sig[y_sig[,1]==a[1]&y_sig[,2]==a[3],]
#   a2=y_sig[y_sig[,1]==a[2]&y_sig[,2]==a[4],]
#   if (as.numeric(a["value.y"]) > 0.05){return(0)}
#   if (as.numeric(a["value.y"]) <= 0.05){
#     return(as.numeric(a1["direct"])*as.numeric(a2["direct"]))
#   }
# })
x$direct=apply(x,1,function(a){
  #browser()
  a1=y_sig[y_sig[,1]==a[1]&y_sig[,2]==a[3],]
  a2=y_sig[y_sig[,1]==a[2]&y_sig[,2]==a[4],]
  if (as.numeric(a["value.y"]) > 0.05){return("grey")}
  if (as.numeric(a["value.y"]) <= 0.05){
    
    if (a1["signif"]!="ns" & a2["signif"] != "ns"){
      if (as.numeric(a1["direct"])*as.numeric(a2["direct"])<0){
        return("blue")
      }else{return("red")}
    }
    if(a1["signif"]=="ns" & a2["signif"] == "ns"){return("darkgrey")}
    if (a1["signif"]!="ns" | a2["signif"] != "ns"){
      if (as.numeric(a1["direct"])*as.numeric(a2["direct"])<0){
        return("lightblue")
      }else{return("pink")}
    }
  }
})
# x_s = subset(x, col!="grey")
###############
# 
# # circus plot
# pdf("circos.plot_dmdp.pdf",width = 4.5,height = 4.5)
# #pdf("circos.plot_dmds.pdf",width = 4.5,height = 4.5)
# par(xpd=T)
# circos.par("track.height" = 0.3,gap.degree = 45)
# circos.initialize(y_sig$tissue, x=y_sig$L1,xlim = c(0.5, 6.5))
# circos.track(ylim = c(0, 1))
# for(sector.index in c("TE","AG","OV")) {
#   #browser()
#   y = subset(y_sig,tissue==sector.index)
#   circos.points(1:6, 0.5,sector.index,cex=y$size,col=y$color,pch=19)
#   circos.axis(sector.index=sector.index,major.at = 1:6,labels=1:6,labels.cex = 0.8,h="bottom",minor.ticks = 0,direction = "inside")
#   circos.text(sector.index=sector.index,
#     CELL_META$xcenter,
#     CELL_META$cell.ylim[2] + convert_y(2, "mm"),
#     labels = sector.index,
#     facing = "bending.inside", cex = 1.2,
#     adj = c(0.5, 0), niceFacing = TRUE
#   )
# }
# for(i in 1:nrow(x_s)){
#   circos.link(x_s[i,1],x_s[i,3],x_s[i,2],x_s[i,4],col=x_s[i,"col"],rou1 = 0.6,rou2 = 0.6,lwd=x_s[i,"value.x"]*5)
# }
# circos.clear()
# dev.off()

#############################
#cols=c(red="red",blue="blue",grey="grey")
cols=c(brewer.pal(9, "Paired")[c(1,2,5,6)],"grey")
names(cols)=c("lightblue","blue","pink","red","grey")
p<- ggplot(subset(x,sec1=="TE"&sec2=="AG")) + geom_point(aes(x=as.factor(X1),y=as.factor(X2),size=value.x*100,color=direct)) + scale_color_manual(values=cols)+theme_bw(10) + xlab("Testis clusters") + ylab("Accessory gland clusters")+scale_size_area(max_size = 10)
ggsave("TEvsAG_cluster_correlation.pdf",width=3.5,height = 2.5)
p<- ggplot(subset(x,sec1=="OV"&sec2=="AG")) + geom_point(aes(x=as.factor(X2),y=as.factor(X1),size=value.x*100,color=direct)) + scale_color_manual(values=cols)+theme_bw(10)+ xlab("Accessory gland clusters") + ylab("Ovary clusters") +scale_size_area(max_size = 10)
ggsave("AGvsOV_cluster_correlation.pdf",width=3.5,height = 2.5)
p<- ggplot(subset(x,sec1=="OV"&sec2=="TE")) + geom_point(aes(x=as.factor(X2),y=as.factor(X1),size=value.x*100,color=direct)) + scale_color_manual(values=cols)+theme_bw(10) + xlab("Testis clusters") + ylab("Ovary clusters")+scale_size_area(max_size = 10)
ggsave("TEvsOV_cluster_correlation.pdf",width=3.5,height = 2.5)
################################
##############################













ix=which(te_ag_p<=0.05,arr.ind=T)
a=data.frame(sec1="TE",ix1=ix[,1],sec2="AG",ix2=ix[,2],ol=te_ag_c[ix])

ix=which(ov_ag_p<=0.05,arr.ind=T)
a=rbind(a,data.frame(sec1="OV",ix1=ix[,1],sec2="AG",ix2=ix[,2],ol=ov_ag_c[ix]))

ix=which(ov_te_p<=0.05,arr.ind=T)
a=rbind(a,data.frame(sec1="OV",ix1=ix[,1],sec2="TE",ix2=ix[,2],ol=ov_te_c[ix]))

a$ifsig1=apply(a,1,function(x){
  #browser()
  x0 = as.numeric(subset(fc_all,L1==0&tissue==x[1])[,1])
  x1 = as.numeric(subset(fc_all,L1==x[2]&tissue==x[1])[,1])
  p = wilcox.test(x0,x1)$p.value
  return(ifelse(p<=0.05,1,-1))
})
a$ifsig2=apply(a,1,function(x){
  #browser()
  x0 = as.numeric(subset(fc_all,L1==0&tissue==x[3])[,1])
  x1 = as.numeric(subset(fc_all,L1==x[4]&tissue==x[3])[,1])
  p = wilcox.test(x0,x1)$p.value
  return(ifelse(p<=0.05,1,-1))
})


a$direct=apply(a,1,function(x){
  #browser()
  x0 = as.numeric(subset(fc_all,L1==0&tissue==x[1])[,1])
  x1 = as.numeric(subset(fc_all,L1==x[2]&tissue==x[1])[,1])
  y0 = as.numeric(subset(fc_all,L1==0&tissue==x[3])[,1])
  y1 = as.numeric(subset(fc_all,L1==x[4]&tissue==x[3])[,1])
  d1 = ifelse(median(x0)>median(x1),-1,1)
  d2 = ifelse(median(y0)>median(y1),-1,1)
  return(d1*d2)
})
a$col = ifelse(a$ix1%in%c(1,2,3)&a$ix2%in%c(1,2,3),"grey",ifelse(a$ifsig1==-1|a$ifsig2==-1,"grey",ifelse(a$direct==-1,"blue","red")))


#write.table(a,"cluster_corresponding.txt",row.names = F, col.names = T,sep="\t",quote=F)

####### circos plot
#a=read.table("cluster_corresponding.txt",header=T)
pdf("circos.plot.pdf",width = 4.5,height = 4.5)
par(xpd=T)
circos.par("track.height" = 0.5,gap.degree = 60)
circos.initialize(fc_all$tissue, x=fc_all$L1,xlim = c(-0.5, 6.5))
circos.track(fc_all$tissue, x = fc_all$L1, y = fc_all$value,ylim=c(-10,10),
             panel.fun = function(x, y) {
               #browser()
               a = by(y,x,function(x)x)
               circos.boxplot(a,0:6,outline = F,box_width = 0.5,lwd=1.2)
               circos.axis(major.at = 0:6,labels=0:6,labels.cex = 0.8,h="bottom",minor.ticks = 0,direction = "inside")
               circos.yaxis(labels.cex = 0.8)
               circos.lines(x=-0.5:6.5,y=rep(0,8),col="red",lty = 2)
               circos.text(
                 CELL_META$xcenter,
                 CELL_META$cell.ylim[2] + convert_y(2, "mm"),
                 CELL_META$sector.index,
                 facing = "bending.inside", cex = 1.2, 
                 adj = c(0.5, 0), niceFacing = TRUE
               )
               circos.text(
                 CELL_META$cell.xlim[1] - convert_x(10, "mm"),
                 CELL_META$ycenter,
                 "log2FC (dme/dps)",
                 facing = "clockwise", cex = 0.65, 
                 adj = c(0.5, 0), niceFacing = TRUE
               )
             })
for(i in 1:nrow(a)){
  circos.link(a[i,1],a[i,2],a[i,3],a[i,4],col=a[i,9],rou1 = 0.42,rou2 = 0.42)
}
circos.clear()
dev.off()

