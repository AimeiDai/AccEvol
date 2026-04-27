setwd("/Users/aimeidai/Library/CloudStorage/OneDrive-个人/10.shangrui/CNE")
library("reshape")
library("ggplot2")
require(RColorBrewer)
library(gghalves)
options(stringsAsFactors = F)

agCNE = subset(read.table("Accgland.CNE.ol.bed"),V7>=-1500&V7<=100&V9!=-1)
teCNE = subset(read.table("Testis.CNE.ol.bed"),V7>=-1500&V7<=100&V9!=-1)
ovCNE = subset(read.table("Ovary.CNE.ol.bed"),V7>=-1500&V7<=100&V9!=-1)

pCs = read.table("merged.CNE.cons.phastCons.bed",stringsAsFactors = F)
pCs$avgPhastcons = pCs$V5/pCs$V6

agpCs = merge(agCNE[,6:10],pCs[,c(1:4,7)],by.x=3:5,by.y=1:3)
ovpCs = merge(ovCNE[,6:10],pCs[,c(1:4,7)],by.x=3:5,by.y=1:3)
tepCs = merge(teCNE[,6:10],pCs[,c(1:4,7)],by.x=3:5,by.y=1:3)

allpCs = rbind(data.frame(agpCs,tissue="AG",type="All"),data.frame(tepCs,tissue="TE",type="All"),data.frame(ovpCs,tissue="OV",type="All"))
allpCs$tissue=factor(allpCs$tissue, levels=c("AG","TE","OV"))
cols <- brewer.pal(9, "Set1")
p1 <- ggplot(allpCs) + geom_density(aes(x=avgPhastcons,color=tissue),linewidth=0.8)+facet_wrap(~type)+theme_classic(10)+scale_color_manual(values=cols)+xlab("Averaged phastCons")+ylab("Density")
p2 <- ggplot(allpCs) + stat_ecdf(aes(x=avgPhastcons,color=tissue),linewidth=0.8,geom="step")+facet_wrap(~type)+theme_classic(10)+scale_color_manual(values=cols)+xlab("Averaged phastCons")+ylab("Empirical cummulative density")
p <- p1+p2
ggsave("whole_avgPhastcons_density_plot.pdf",width = 5.5,height = 3)
ks.test(agpCs[,7],ovpCs[,7]) #D = 0.052238, p-value = 0.01295
ks.test(tepCs[,7],ovpCs[,7]) #D = 0.051615, p-value = 0.01501
ks.test(agpCs[,7],tepCs[,7]) #D = 0.014008, p-value = 0.989

########
ifag = apply(pCs[,1:3],1,function(x){paste(x,collapse = ":")}) %in% apply(agCNE[,8:10],1,function(x){paste(x,collapse = ":")})
ifte = apply(pCs[,1:3],1,function(x){paste(x,collapse = ":")}) %in% apply(teCNE[,8:10],1,function(x){paste(x,collapse = ":")})
ifov = apply(pCs[,1:3],1,function(x){paste(x,collapse = ":")}) %in% apply(ovCNE[,8:10],1,function(x){paste(x,collapse = ":")})

ix=apply(cbind(as.numeric(ifag),as.numeric(ifte),as.numeric(ifov)),1,function(x){paste(x,collapse = "")})

pCs$ix=factor(ix,levels=c("000","111","110","101","011","100","010","001"))

mds = melt(tapply(pCs$avgPhastcons,pCs$ix,median))
mds$ix=ifelse(mds$indices==0,"000",ifelse(mds$indices==11,"011",ifelse(mds$indices=="10","010",ifelse(mds$indices=="1","001",mds$indices))))

cols1=c("darkgrey",cols[1],rep(cols[2],3),rep(cols[3],3))
p <- ggplot() + geom_violin(data=pCs,aes(x=ix,y=avgPhastcons,color=ix),width=1.2)+geom_point(data=mds,aes(x=ix,y=value,colour = ix))+scale_color_manual(values=cols1) + theme_bw(10)+coord_cartesian(ylim=c(0,1.1)) + ylab("Averaged Phastcons")+ coord_flip()
ggsave("avgPhastcons_vs_CNE_distribution.pdf",width = 2.8,height = 3.5)

ixes=combn(c("001","010","100","011","101","110","111","000"),2)

pval=matrix(NA,nrow = 8,ncol=8,dimnames = list(c("001","010","100","011","101","110","111","000"),c("001","010","100","011","101","110","111","000")))
for(i in 1:ncol(ixes)){
  x1 = subset(pCs,ix==ixes[1,i])$avgPhastcons
  x2 = subset(pCs,ix==ixes[2,i])$avgPhastcons
  pval[ixes[1,i],ixes[2,i]]=wilcox.test(x1,x2)$p.value
}
ifsig = ifelse(pval<=0.001,"***",ifelse(pval<=0.01,"**",ifelse(pval<=0.05,"*","ns")))
cols2=c("grey",brewer.pal(9, "Reds")[c(1,3,5,7)])
names(cols2)=c("NA","ns","*","**","***")

library(ComplexHeatmap)
pdf("WilcoxTest.heatmap.pdf",width=3,height = 2)
h=draw(Heatmap(ifsig, name = "Significance", col = cols2,rect_gp = gpar(col = "white", lwd = 1), show_row_names = T, show_column_names = T,cluster_rows = F,cluster_columns=F,show_column_dend =F,show_row_dend =F))
dev.off()
mds
# p <- ggplot() + geom_half_violin(data=pCs,aes(x=ix,y=avgPhastcons,color=ix,fill = ix),width=1.8,alpha = 0.5,side = "r")+geom_point(data=mds,aes(x=ix,y=value,colour = ix))+scale_color_manual(values=cols1)+scale_fill_manual(values=cols1) + theme_bw(10)+coord_cartesian(ylim=c(0,1.1)) + ylab("Averaged Phastcons") + coord_flip()
# ggsave("avgPhastcons_vs_CNE_distribution2.pdf",width = 2.5,height = 3)
#+geom_boxplot(data=pCs,aes(x=ix,y=avgPhastcons,color=ix),width = 0.2, cex = 0.5, alpha = 0.5, position = position_nudge(x = 0.22, y = 0))
