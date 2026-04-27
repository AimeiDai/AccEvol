setwd("/Users/aimeidai/Library/CloudStorage/OneDrive-个人/10.shangrui/remove_batch_effect")
library(sva)
library(pamr)
library(limma)
library(ggplot2)
library(reshape)
library(edgeR)
library(ape)
options(stringsAsFactors=F)
require(RColorBrewer)
cols1 <- c("black","grey",rev(brewer.pal(12, "Paired")[c(1:6,10,8)]))
names(cols1)=c("EA","WI","LB","AG","GO.f","GO.m","OV","TE","HD.f","HD.m")

cols2 <- c(rev(brewer.pal(12, "Paired")[c(1:6,10)]))
names(cols2)=c("AG","GO.f","GO.m","OV","TE","HD.f","HD.m")


sample_table = read.table("sample_table_rc.txt",sep="\t",header=T)
rownames(sample_table) = sample_table[,3]

count_matrix = as.matrix(read.table("raw_counts_5_species.txt",header=T, sep="\t"))
colnames(count_matrix) = sample_table[,3]

############# total samples
y <- DGEList(counts=count_matrix, group=sample_table$merged_SRR)
keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
a = plotMDS(y, plot=F)
summary(prcomp(a$distance.matrix))
a_gg = data.frame(a$x,a$y)
rownames(a_gg)=rownames(a$distance.matrix.squared)
a_gg$tissue = factor(sample_table[rownames(a_gg),"tissue"])
a_gg$species = factor(sample_table[rownames(a_gg),"species"],levels=c("dme","dsi","dya","dan","dps"))
p <- ggplot(a_gg, aes(x=a.x, y= a.y, color=tissue, shape=species)) + 
  geom_point(size=2,stroke=1.2) + scale_color_manual(values=cols1)  + theme_bw(10) +
  scale_shape_manual(values=0:9) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text=element_text(size=8),legend.title=element_blank(), legend.position="right") + 
  labs(title="",x="PCoA1: 49.70% variance", y="PCoA2: 32.16% variance")
ggsave("PCA_with_batch_8tissue.pdf",width=6,height=5.2)


##########
tmp_s = subset(sample_table,batch!=10)
tmp = count_matrix[,tmp_s[,3]]
########before remove batch effect
y <- DGEList(counts=tmp, group=tmp_s$merged_SRR)
keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
a = plotMDS(y, plot=F)
summary(prcomp(a$distance.matrix))
a_gg = data.frame(a$x,a$y)
rownames(a_gg)=rownames(a$distance.matrix.squared)
a_gg$tissue = factor(tmp_s[rownames(a_gg),"tissue"])
a_gg$species = factor(tmp_s[rownames(a_gg),"species"],levels=c("dme","dsi","dya","dan","dps"))
p <- ggplot(a_gg, aes(x=a.x, y= a.y, color=tissue, shape=species)) + 
  geom_point(size=2,stroke=1.2) + scale_color_manual(values=cols2)  + theme_bw(10) +
  scale_shape_manual(values=0:9) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text=element_text(size=8),legend.title=element_blank(), legend.position="right") + 
  labs(title="",x="PCoA1: 51.46% variance", y="PCoA2: 31.45% variance")
ggsave("PCA_with_batch_4tissue.pdf",width=6,height=5.2)

### remove batch effects
#covar_mat <- model.matrix(~tissue,data=tmp_s)
adjusted <- ComBat_seq(tmp, batch=tmp_s$batch, group=NULL, covar_mod = tmp_s[,c("merged_tissue2","species")])

y <- DGEList(counts=adjusted, group=tmp_s$merged_SRR)
keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
a = plotMDS(y, plot=F)
summary(prcomp(a$distance.matrix))
a_gg = data.frame(a$x,a$y)
rownames(a_gg)=rownames(a$distance.matrix.squared)
a_gg$tissue = factor(tmp_s[rownames(a_gg),"tissue"])
a_gg$species = factor(tmp_s[rownames(a_gg),"species"],levels=c("dme","dsi","dya","dan","dps"))
p <- ggplot(a_gg, aes(x=a.x, y= a.y, color=tissue, shape=species)) + 
  geom_point(size=2,stroke=1.2) + scale_color_manual(values=cols2)  + theme_bw(10) +
  scale_shape_manual(values=0:9) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text=element_text(size=8),legend.title=element_blank(), legend.position="right") + 
  labs(title="",x="PCoA1: 50.31% variance", y="PCoA2: 31.66% variance")
ggsave("PCA_remove_batch_4tissues.pdf",width=6,height=5.2)

#############################
x = subset(sample_table,batch==10)
adj_count=cbind(adjusted, count_matrix[,x[,3]])
write.table(adj_count[,sample_table[,3]],"adjusted_raw_count_5species.txt",row.names = T, col.names = T, sep="\t",quote=F)


y <- DGEList(counts=adj_count, group=sample_table$merged_SRR)
keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

cpm_adj = sweep(adj_count,2,colSums(adj_count),"/")*10^6
write.table(cpm_adj[,sample_table[,3]],"CPM_adjusted_5species.txt",row.names = T, col.names = T, sep="\t",quote=F)

##########################################
#cpm_adj = read.table("CPM_adjusted_5species.txt",header=T)
#colnames(cpm_adj)=sample_table[,3]
tpm = cpm_adj
#############################
# tmp = count_matrix[,x[,3]]
# y <- DGEList(counts=tmp, group=x$merged_SRR)
# keep <- filterByExpr(y)
# table(keep)
# y <- y[keep, , keep.lib.sizes=FALSE]
# y <- calcNormFactors(y)
# tpm = cpm(y$counts)
#################################
# M = 1-cor(tpm,method="spearman")
# tr <- nj(M)
# pdf("nj_expression.pdf",width=7,height=7)
# plot(tr,"u",TRUE, label.offset = 1,cex=0.5)
# add.scale.bar(cex = 0.7, font = 2)
# layout(1)
# dev.off()


##############################
### a less theoretical example
## bootstrapping
cal_tl_bt <- function(expr,times=1000){
  tl=vector(length=times)
  for(i in 1:times){
    ix = sample(1:nrow(expr),0.85*nrow(expr))
    M = 1-cor(expr[ix,],method="spearman")
    tr <- nj(M)
    tl[i]=sum(tr$edge.length)
  }
  return(tl)
}

divergence <- function(x){
  #return(sd(x)/mean(x))
  # y = apply(x, 2, function(a){
  #   m=mean(a)
  #   s=sd(a)
  #   return((a-m)/s)
  # })
  return(rowSums(abs(x[,"dps"]-x)))
}

#####
tpm_var=matrix(nrow=nrow(tpm),ncol=8,dimnames=list(rownames(tpm),c("WI","LB","EA","HD.f","HD.m","TE","AG","OV")))

pdf("nj_by_tissue.pdf",width=5.5,height=3.8)
layout(matrix(1:8,2,4,byrow=T))
par(mar=c(1.5, 1.5, 1.5, 1.5), xpd=TRUE)
#layout.show(n=8)
# wing imaginal disc
st = subset(sample_table,merged_tissue=="WI")
tpm_s = tpm[,st[,3]]

tpm_m = do.call("cbind",tapply(st[,3],st[,4],function(x){
  #browser()
  if(length(x)==1){
    return(tpm_s[,x])
  }
  y = apply(tpm_s[,x],1,mean)
  return(y)
}))
tpm_var[,"WI"] = divergence(tpm_m)

M = 1-cor(tpm_m, method="spearman")
tr <- nj(M)
#pdf("nj_WI.pdf",width=7,height=7)
tr=rotate(tr,c("dan","dps"))
plot(tr,"u",TRUE, xpd=T,cex=1, edge.width = 2, main="Wing imaginal disc")
add.scale.bar(cex = 1, lwd=2, length=0.02, font = 2)
#dev.off()
sum(tr$edge.length)
# 0.1569415
tot_tl = data.frame(totTreeLen=cal_tl_bt(tpm_m),tissue="WI")

#### LB
st = subset(sample_table,merged_tissue=="LB")
tpm_s = tpm[,st[,3]]
tpm_m = do.call("cbind",tapply(st[,3],st[,4],function(x){
  #browser()
  if(length(x)==1){
    return(tpm_s[,x])
  }
  y = apply(tpm_s[,x],1,mean)
  return(y)
}))
tpm_var[,"LB"] = divergence(tpm_m)

M = 1-cor(tpm_m, method="spearman")
tr <- nj(M)
#pdf("nj_LB.pdf",width=7,height=7)
tr=rotate(tr,c("dan","dps"))
#tr=reorder(tr,order = "pruningwise")
plot(tr,"u",TRUE, cex=1, edge.width = 2, main="Larval brain")
add.scale.bar(cex = 1, lwd=2, length=0.02, font = 2)
#dev.off()
sum(tr$edge.length)
#  0.1550548
tot_tl = rbind(tot_tl,data.frame(totTreeLen=cal_tl_bt(tpm_m),tissue="LB"))

####male head
st = subset(sample_table,merged_tissue=="HD.m")
tpm_s = tpm[,st[,3]]
tpm_m = do.call("cbind",tapply(st[,3],st[,4],function(x){
  #browser()
  y = apply(tpm_s[,x],1,mean)
  return(y)
}))
tpm_var[,"HD.m"] = divergence(tpm_m)

M = 1-cor(tpm_m, method="spearman")
tr <- nj(M)
#pdf("nj_HD.pdf",width=7,height=7)
tr=rotate(tr,c("dan","dps"))
plot(tr,"u",TRUE, cex=1, edge.width = 2, main="Male head")
add.scale.bar(cex = 1, lwd=2, length=0.02, font = 2)
#dev.off()
sum(tr$edge.length)
# 0.2397096
tot_tl = rbind(tot_tl,data.frame(totTreeLen=cal_tl_bt(tpm_m),tissue="HD.m"))

####female head
st = subset(sample_table,merged_tissue=="HD.f")
tpm_s = tpm[,st[,3]]
tpm_m = do.call("cbind",tapply(st[,3],st[,4],function(x){
  #browser()
  y = apply(tpm_s[,x],1,mean)
  return(y)
}))
#tpm_var[,"HD.f"] = apply(tpm_m,1,divergence)
tpm_var[,"HD.f"] = divergence(tpm_m)

M = 1-cor(tpm_m, method="spearman")
tr <- nj(M)
#pdf("nj_HD.pdf",width=7,height=7)
tr=rotate(tr,c("dan","dps"))
plot(tr,"u",TRUE, cex=1, edge.width = 2, main="Female head")
add.scale.bar(cex = 1, lwd=2, length=0.02, font = 2)
#dev.off()
sum(tr$edge.length)
# 0.1807609
tot_tl = rbind(tot_tl,data.frame(totTreeLen=cal_tl_bt(tpm_m),tissue="HD.f"))

#### EA
st = subset(sample_table,merged_tissue=="EA")
tpm_s = tpm[,st[,3]]
tpm_m = do.call("cbind",tapply(st[,3],st[,4],function(x){
  #browser()
  if(length(x)==1){
    return(tpm_s[,x])
  }
  y = apply(tpm_s[,x],1,mean)
  return(y)
}))
#tpm_var[,"EA"] = apply(tpm_m,1,divergence)
tpm_var[,"EA"] = divergence(tpm_m)

M = 1-cor(tpm_m, method="spearman")
tr <- nj(M)
#pdf("nj_EA.pdf",width=7,height=7)
tr=rotate(tr,c("dan","dps"))
plot(tr,"u",TRUE, cex=1, edge.width = 2, main="Eye-antennal imaginal disc")
add.scale.bar(cex = 1, lwd=2, length=0.02, font = 2)
#dev.off()
sum(tr$edge.length)
# 0.1730053
tot_tl = rbind(tot_tl,data.frame(totTreeLen=cal_tl_bt(tpm_m),tissue="EA"))

#### ovary
st = subset(sample_table,merged_SRR!="SRR330569"&(merged_tissue=="OV"))
tpm_s = tpm[,st[,3]]
tpm_m = do.call("cbind",tapply(st[,3],st[,4],function(x){
  #browser()
  y = apply(tpm_s[,x],1,mean)
  return(y)
}))
#tpm_var[,"OV"] = apply(tpm_m,1,divergence)
tpm_var[,"OV"] = divergence(tpm_m)

M = 1-cor(tpm_m, method="spearman")
tr <- nj(M)
#pdf("nj_OV.pdf",width=7,height=7)
tr=rotate(tr,c("dan","dps"))
plot(tr,"u",TRUE, cex=1, edge.width = 2, main="Ovary")
add.scale.bar(cex = 1, lwd=2, length=0.02, font = 2)
#dev.off()
sum(tr$edge.length)
# 0.145766
tot_tl = rbind(tot_tl,data.frame(totTreeLen=cal_tl_bt(tpm_m),tissue="OV"))

#### testis
st = subset(sample_table,merged_tissue=="TE")
tpm_s = tpm[,st[,3]]
tpm_m = do.call("cbind",tapply(st[,3],st[,4],function(x){
  #browser()
  y = apply(tpm_s[,x],1,mean)
  return(y)
}))
#tpm_var[,"TE"] = apply(tpm_m,1,divergence)
tpm_var[,"TE"] = divergence(tpm_m)

M = 1-cor(tpm_m, method="spearman")
tr <- nj(M)
#pdf("nj_TE.pdf",width=7,height=7)
tr=rotate(tr,c("dan","dps"))
plot(tr,"u",TRUE,cex=1, edge.width = 2, main="Testis")
add.scale.bar(cex = 1, lwd=2, length=0.02, font = 2)
#dev.off()
sum(tr$edge.length)
# 0.3593882
tot_tl = rbind(tot_tl,data.frame(totTreeLen=cal_tl_bt(tpm_m),tissue="TE"))

#### AG
st = subset(sample_table,merged_tissue=="AG")
tpm_s = tpm[,st[,3]]
tpm_m = do.call("cbind",tapply(st[,3],st[,4],function(x){
  #browser()
  y = apply(tpm_s[,x],1,mean)
  return(y)
}))
#tpm_var[,"AG"] = apply(tpm_m,1,divergence)
tpm_var[,"AG"] = divergence(tpm_m)

M = 1-cor(tpm_m, method="spearman")
tr <- nj(M)
#pdf("nj_AG.pdf",width=7,height=7)
tr=rotate(tr,c("dan","dps"))
plot(tr,"u",FALSE, cex=1, edge.width = 2, main="Accessory gland")
add.scale.bar(cex = 1, lwd=2, length=0.02, font = 2)
#dev.off()
sum(tr$edge.length)
# 0.358567
tot_tl = rbind(tot_tl,data.frame(totTreeLen=cal_tl_bt(tpm_m),tissue="AG"))

dev.off()

tot_tl$tissue = factor(tot_tl$tissue, levels=c("WI","LB","EA","HD.f","HD.m","OV","TE","AG"))
medians=data.frame(median=tapply(tot_tl[,1],tot_tl[,2],median),max=tapply(tot_tl[,1],tot_tl[,2],max))
medians$ix = rownames(medians)
p <- ggplot() + geom_boxplot(data=tot_tl, aes(x=tissue,y=totTreeLen)) + geom_text(data=medians,aes(x=ix,y=max+0.02,label=round(median,4)),size=3.5) + theme_minimal(12) +  ylab("Total tree length") + theme(axis.text.x=element_text(angle=45,hjust=1)) + ylim(c(0.1,0.4))
ggsave("totTreeLen_bootstrapping.pdf",width=3.8,height=4)

comparisons = combn(unique(tot_tl$tissue),2)
p_vals=rep(NA,ncol(comparisons))
for (i in 1:ncol(comparisons)) {
  x1 = subset(tot_tl,tissue==comparisons[1,i])$totTreeLen
  x2 = subset(tot_tl,tissue==comparisons[2,i])$totTreeLen
  p_vals[i]=wilcox.test(x1,x2)$p.value
}
p.adjust(p_vals)
tpm_var[which(is.nan(tpm_var),arr.ind=T)]=0
write.table(tpm_var,"divergenceByTissue.txt",row.names=T,col.names=T,sep="\t",quote=F)

library(ComplexHeatmap)
set.seed(1)
mat = as.matrix(read.table("divergenceByTissue.txt",header = T))
pdf("heatmapDivergenceByTissue.pdf",height = 5,width=3.5)
ht=draw(Heatmap(log10(mat+1),name = "log10Div", km =6, show_row_names = FALSE, show_column_names = TRUE))
dev.off()



