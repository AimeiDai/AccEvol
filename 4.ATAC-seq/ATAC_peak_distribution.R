library(ggplot2)
library(reshape)
library(ComplexHeatmap)
library(edgeR)
options(stringsAsFactors = F)
#setwd("C:/Users/Administrator/OneDrive/10.shangrui/ATAC-seq")
setwd("/Users/aimeidai/Library/CloudStorage/OneDrive-个人/10.shangrui/ATAC-seq")
require(RColorBrewer)
library(circlize)
###########
tab = subset(read.table("readCounts.tab",header = T),chr %in% c("2L","2R","3L","3R","4","X"))
tab = tab[,c(1:9,12:14)]
y <- DGEList(counts=tab[,4:ncol(tab)])
keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
y_rpm = cpm(y)
corr = cor(y_rpm)
col_fun = colorRamp2(c(0.6, 0.8, 1), c("white","pink", "red"))
pdf("ATAC.cor.heatmap.pdf",width=3.5,height = 2.8)
h=draw(Heatmap(corr,name = "Pearson's r", show_row_names = T, show_column_names = TRUE,cluster_columns=F,cluster_rows=F,col = col_fun))
dev.off()

corr_gg = melt(data.frame(AGvsTE=c(corr[1:3,7:9]),AGvsOV=c(corr[1:3,4:6]),TEvsOV=c(corr[4:6,7:9])))
p <- ggplot() + geom_boxplot(data=corr_gg,aes(y=value, x=variable)) + theme_bw(7) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab("Pearson's r")+ylim(0.58,1.05)+coord_flip()
ggsave("correlation_boxplot.pdf",width=1.5,height=0.8)

x1=subset(corr_gg,variable=="AGvsTE")[,2]
x2=subset(corr_gg,variable=="AGvsOV")[,2]
x3=subset(corr_gg,variable=="TEvsOV")[,2]
wilcox.test(x1,x2)$p.value
wilcox.test(x1,x3)$p.value
wilcox.test(x2,x3)$p.value

a = plotMDS(y, plot=F)
summary(prcomp(a$distance.matrix))
a_gg = data.frame(a$x,a$y)
rownames(a_gg)=rownames(a$distance.matrix.squared)
a_gg$tissue = c(rep("AG",3),rep("OV",3),rep("TE",3))
a_gg$replicate = rep(1:3,times=3)
p <- ggplot(a_gg, aes(x=a.x, y= a.y, color=tissue ))+ theme_bw(10) + 
  geom_text(aes(label=as.character(replicate)),size=6) + scale_color_brewer(palette="Set1")   + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text=element_text(size=6),legend.title=element_blank(), legend.position="right") + ylim(-2,2) +  labs(title="",x="PCoA1: 84.47% variance", y="PCoA2: 10.14% variance")
ggsave("PCAamongSamples.pdf",width=3.3,height=2.6)
############## 
## loading packages
library(ChIPseeker)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
library(clusterProfiler)
#########
tePeak = readPeakFile("Testis-overlapping.2rep")
ovPeak = readPeakFile("Ovary-overlapping.2rep")
agPeak = readPeakFile("Accgland-overlapping.2rep")
###### peak annotation
tePeakAnno <- annotatePeak(tePeak, tssRegion=c(-3000, 100),
                         TxDb=txdb, annoDb="org.Dm.eg.db")
ovPeakAnno <- annotatePeak(ovPeak, tssRegion=c(-3000, 100),
                           TxDb=txdb, annoDb="org.Dm.eg.db")
agPeakAnno <- annotatePeak(agPeak, tssRegion=c(-3000, 100),
                           TxDb=txdb, annoDb="org.Dm.eg.db")
sum(tePeakAnno@annoStat$Frequency[1:3])
sum(ovPeakAnno@annoStat$Frequency[1:3])
sum(agPeakAnno@annoStat$Frequency[1:3])
plotAnnoPie(tePeakAnno)
plotAnnoPie(ovPeakAnno)
plotAnnoPie(agPeakAnno)


# pdf("PeaksTSSDistribution.pdf",width = 6,height = 1)
# plotDistToTSS(tePeakAnno, title="Distribution of testis ATAC-seq peaks relative to TSS")
# plotDistToTSS(ovPeakAnno, title="Distribution of ovary ATAC-seq peaks relative to TSS")
# plotDistToTSS(agPeakAnno, title="Distribution of accessory gland ATAC-seq peaks relative to TSS")
# dev.off()

files=list(AG="Accgland-overlapping.2rep",TE="Testis-overlapping.2rep",OV="Ovary-overlapping.2rep")
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
pdf("PeaksGenomicDistribution.pdf",width = 5,height = 1.8)
plotAnnoBar(peakAnnoList, title="Genomic distribution of ATAC-seq peaks")
plotDistToTSS(peakAnnoList,title="Distribution of ATAC-seq peaks relative to TSS")
dev.off()

#################
# peak distribution
agPeak = unique(subset(read.table("../CNE/Accgland.CNE.ol.bed"),V7>=-1500&V7<=100)[,c(1:7)])
tePeak = unique(subset(read.table("../CNE/Testis.CNE.ol.bed"),V7>=-1500&V7<=100)[,c(1:7)])
ovPeak = unique(subset(read.table("../CNE/Ovary.CNE.ol.bed"),V7>=-1500&V7<=100)[,c(1:7)])

######## ATACpeaks assigned to genes by promoter region
################ peak number distribution
# agn = tapply(agPeak[,6],agPeak[,14],length)
# ten = tapply(tePeak[,6],tePeak[,14],length)
# ovn = tapply(ovPeak[,6],ovPeak[,14],length)
##################
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
###############
setwd("/Users/aimeidai/Library/CloudStorage/OneDrive-个人/10.shangrui/ATAC-seq")

cal_expr_divergence <- function(expr,genelist,times=1000){
  #browser()
  rhos=rep(NA,times)
  for(i in 1:times){
    x = sample(genelist,length(genelist)*0.8)
    rhos[i]= 1-cor(expr[x,1],expr[x,2],method = "spearman")
  }
  return(rhos)
}
expr_div = rbind(data.frame(value=cal_expr_divergence(tpm_m[,c("dme.AG","dps.AG")],intersect(unique(agPeak[,"V6"]),rownames(tpm_m))),type="+",tissue="AG"),data.frame(value=cal_expr_divergence(tpm_m[,c("dme.AG","dps.AG")],setdiff(rownames(tpm_m),unique(agPeak[,"V6"]))),type="-",tissue="AG"),
                 data.frame(value=cal_expr_divergence(tpm_m[,c("dme.TE","dps.TE")],intersect(unique(tePeak[,"V6"]),rownames(tpm_m))),type="+",tissue="TE"),data.frame(value=cal_expr_divergence(tpm_m[,c("dme.TE","dps.TE")],setdiff(rownames(tpm_m),unique(tePeak[,"V6"]))),type="-",tissue="TE"),
                 data.frame(value=cal_expr_divergence(tpm_m[,c("dme.OV","dps.OV")],intersect(unique(ovPeak[,"V6"]),rownames(tpm_m))),type="+",tissue="OV"),data.frame(value=cal_expr_divergence(tpm_m[,c("dme.OV","dps.OV")],setdiff(rownames(tpm_m),unique(ovPeak[,"V6"]))),type="-",tissue="OV"))

expr_div$tissue = factor(expr_div$tissue,levels=c("AG","TE","OV"))
p<-ggplot(expr_div) + geom_boxplot(aes(x=tissue,y=value,colour =interaction(type,tissue)),outlier.size = 0.5)+ theme_bw(12) + ylab(expression(paste("Expression divergence (1-",rho,")",sep="")))+ theme(axis.text.x=element_text(angle=45,hjust=1)) + scale_color_manual(values=brewer.pal(9, "Paired")[c(5,6,3,4,1,2)])
ggsave("peakEffect_in_evolutiondivergence.pdf",width = 4.5,height = 3.5)
by(expr_div[,1:2],expr_div[,3],function(x){
  wilcox.test(x[,1]~x[,2])
})
by(expr_div[,1],expr_div[,2:3],function(x){
  median(x)
})
####################
fc_all = read.table("../coevolution/logFC_dmVSdp.txt")
fc_all$ifAG = ifelse(rownames(fc_all) %in% agPeak[,6],"+","-")
fc_all$ifTE = ifelse(rownames(fc_all) %in% tePeak[,6],"+","-")
fc_all$ifOV = ifelse(rownames(fc_all) %in% ovPeak[,6],"+","-")

fc_gg = rbind(data.frame(log2FC = fc_all[,1],ifPeak=fc_all[,4],Organ="AG"),data.frame(log2FC = fc_all[,2],ifPeak=fc_all[,5],Organ="TE"),data.frame(log2FC = fc_all[,3],ifPeak=fc_all[,6],Organ="OV"))
fc_gg$Organ = factor(fc_gg$Organ,levels=c("AG","TE","OV"))
p<-ggplot(fc_gg) + geom_boxplot(aes(x=Organ,y=log2FC,colour =interaction(ifPeak,Organ)),outlier.size = 0.5)+ theme_bw(12) + ylab("log2FC (dme/dps)")+ theme(axis.text.x=element_text(angle=45,hjust=1)) + scale_color_manual(values=brewer.pal(9, "Paired")[c(5,6,3,4,1,2)])+coord_cartesian(ylim=c(-10,20))
ggsave("peakEffect_in_evolution.pdf",width = 5,height = 3.5)
tapply(fc_gg[,1],fc_gg[,2:3],function(x){median(x,na.rm=T)})
# Organ
# ifPeak         AG        TE           OV
# - -0.4353580 0.2218161 -0.002730526
# + -0.2600045 0.3554009  0.074995883
tapply(fc_gg[,1:2],fc_gg[,3],function(x){wilcox.test(x[,1]~x[,2])$p.value})
# AG           TE           OV 
# 1.797402e-08 1.617351e-03 6.776337e-05 
#####################
peak_frac = data.frame()
#TE
tec1 = read.table("../coevolution/TE_cluster.txt",sep="\t")
tec2 = read.table("../coevolution/DEG_TE_dmdp.txt",sep="\t")
tec = merge(tec2,tec1,by.x=1,by.y=1)
tec$ifPeak = ifelse(tec[,1]%in%tePeak[,"V6"],"Yes","No")
x=table(tec[,c(8,9)])
x = rbind(x,"0"=colSums(x))
x0=x["0",]
for(i in 1:6){
  xi = x[as.character(i),]
  m = rbind(x0-xi,xi)
  cat(fisher.test(m)$p.value," ")
}
#1.605526e-21  0.06177427  1.536403e-10  8.894508e-10  0.1415541  1.939804e-08 
peak_frac=data.frame(proportion=x[,2]/(x[,1]+x[,2]),organ="TE",Cluster=rownames(x))

#OV
ovc1 = read.table("../coevolution/OV_cluster.txt",sep="\t")
ovc2 = read.table("../coevolution/DEG_OV_dmdp.txt",sep="\t")
ovc = merge(ovc2,ovc1,by.x=1,by.y=1)
ovc$ifPeak = ifelse(ovc[,1]%in%ovPeak[,"V6"],"Yes","No")
#table(ovc[,c(8,7,9)])
x=table(ovc[,c(8,9)])
x = rbind(x,"0"=colSums(x))
x0=x["0",]
for(i in 1:6){
  xi = x[as.character(i),]
  m = rbind(x0-xi,xi)
  cat(fisher.test(m)$p.value," ")
}
#4.441535e-14  0.03646776  5.530267e-07  0.005200785  0.2366445  0.6006852 
peak_frac=rbind(peak_frac,data.frame(proportion=x[,2]/(x[,1]+x[,2]),organ="OV",Cluster=rownames(x)))

#AG
agc1 = read.table("../coevolution/AG_cluster.txt",sep="\t")
agc2 = read.table("../coevolution/DEG_AG_dmdp.txt",sep="\t")
agc = merge(agc2,agc1,by.x=1,by.y=1)
agc$ifPeak = ifelse(agc[,1]%in%agPeak[,"V6"],"Yes","No")
#table(agc[,c(8,7,9)])
x=table(agc[,c(8,9)])
x = rbind(x,"0"=colSums(x))
x0=x["0",]
for(i in 1:6){
  xi = x[as.character(i),]
  m = rbind(x0-xi,xi)
  cat(fisher.test(m)$p.value," ")
}
#1.031532e-38  4.891524e-40  0.004649907  4.579298e-20  4.486828e-06  0.0001222739  
peak_frac=rbind(peak_frac,data.frame(proportion=x[,2]/(x[,1]+x[,2]),organ="AG",Cluster=rownames(x)))

########
p <- ggplot(data=peak_frac) + geom_line(aes(x=Cluster,y=proportion,group=organ),color="darkgrey",linewidth=0.4) + geom_point(aes(x=Cluster,y=proportion,color=organ),size=4)+ theme_bw(12) + ylab("Proportion of genes\nwith ATAC-seq peaks")+xlab("Coevolution clusters")+scale_color_manual(values = brewer.pal(3, "Set1")) + ylim(c(0,0.5))
ggsave("peakDistributionAlongWithClusters.pdf",width = 4.5,height = 3)






##############
rcl.s = read.table("../ConsDivEvolution/coevolve_state_by_species.txt",header=F)
rownames(rcl.s)=rcl.s[,2]
table(rcl.s$V1)

all_p =  read.table("all_organ.narrowPeak.genomic.region",sep="\t",quote="",header=T)
all_p2 = subset(all_p, Distance.to.TSS>=-1500&Distance.to.TSS<=100)

all_p2$rcl.s = rcl.s[all_p2$Entrez.ID,1]
all_p2 = subset(all_p2,!is.na(rcl.s))
table(all_p2$Peak.Score)
n=table(all_p2[,c("rcl.s","Peak.Score")])
p = sweep(n,2,colSums(n),"/")
######################
rcl.s$ifAGPeak = ifelse(rcl.s[,2] %in% agPeak[,12],"+","-")
rcl.s$ifTEPeak = ifelse(rcl.s[,2] %in% tePeak[,12],"+","-")
rcl.s$ifOVPeak = ifelse(rcl.s[,2] %in% ovPeak[,12],"+","-")

rcl.s$type=apply(rcl.s[,c(3:5)],1,function(x){paste("AG",x[1],"TE",x[2],"OV",x[3],sep = "")})
x=rbind(table(rcl.s[,c(1,6)]),"Whole"=table(rcl.s[,6]))
r=x/rowSums(x)
x_gg = melt(r)

cols <- rev(brewer.pal(9, "Set1"))[1:8]
names(cols)=c("AG-TE-OV-","AG-TE-OV+","AG-TE+OV-","AG-TE+OV+","AG+TE-OV-","AG+TE-OV+","AG+TE+OV-","AG+TE+OV+")
x_gg$X1=factor(x_gg$X1,levels=c("Whole",1:6))
p <- ggplot(x_gg) + geom_bar(aes(x=X1,y=value,fill=X2),stat="identity",color="white")+ theme_bw(12) + theme(axis.text.x=element_text(angle=0,hjust=0.5,vjust = 0.5), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="right", legend.title=element_blank()) + labs(title="",y="Proportion of peaks", x="Coevolution cluster") +scale_fill_manual(values=cols)
ggsave("peakInCoevolveCluster.distribution.pdf",width=4.2,height=3.4)

s0 = x["Whole",]
p_vals=matrix("n.s.",nrow=ncol(x),ncol=nrow(x),dimnames=list(colnames(x),c("Whole",1:6)))
for (i in c("Whole",1:6)) {
  if (i == "Whole"){next}
  s1 = x[i,]
  #print(i)
  #print(chisq.test(cbind(s0,s1))$p.value)
  for(j in 1:ncol(x)){
    #  cat(colnames(x)[j],": ")
    m=cbind(c(s1[j],sum(s1[-j])),c(s0[j],sum(s0[-j])))
    p=fisher.test(m)$p.value
    p_vals[j,i]=ifelse(p<=0.001,"***",ifelse(p<=0.01,"**",ifelse(p<=0.05,"*","n.s.")))
  }
}

cols=c("white",brewer.pal(9, "Reds")[c(3,5,7)])
names(cols)=c("n.s.","*","**","***")
r=r[c("Whole",1:6),]
pdf("FisherTest.heatmap.pdf",width=4,height = 2.5)
h=draw(Heatmap(p_vals, name = "Proportion", col = cols, show_row_names = T, show_column_names = T,cluster_rows = F,cluster_columns=F,show_column_dend =F,show_row_dend =F,
               cell_fun = function(j,i, x, y, width, height, fill) {
                 grid.text(sprintf("%.4f", r[j,i]), x,y, gp = gpar(fontsize = 6))}))
dev.off()
############# k2p values
setwd("../Ks")
ovK2p = read.table("Ovary.narrowPeak.2rep.sorted.k2p",header = T)
ov_cl = merge(ovPeak[,c(1,12)],ovK2p[,c(1,5)],by.x=1,by.y=1)

agK2p = read.table("Accgland.narrowPeak.2rep.sorted.k2p",header = T)
ag_cl = merge(agPeak[,c(1,12)],agK2p[,c(1,5)],by.x=1,by.y=1)

teK2p = read.table("Testis.narrowPeak.2rep.sorted.k2p",header = T)
te_cl = merge(tePeak[,c(1,12)],teK2p[,c(1,5)],by.x=1,by.y=1)

all_cl = rbind(data.frame(ov_cl[,2:3],Organ="OV"),data.frame(ag_cl[,2:3],Organ="AG"),data.frame(te_cl[,2:3],Organ="TE"))
all_cl$type = rcl.s[all_cl[,1],1]
#all_cl$type = factor(all_cl$type,levels=c(1:6))
all_cl = rbind(all_cl,data.frame(all_cl[,1:3],type="0"))
all_cl$Organ = factor(all_cl$Organ,levels=c("AG","TE","OV"))
cols <- c(rev(brewer.pal(9, "Set1"))[c(1,8,9)],brewer.pal(8, "Dark2")[1:8])
#cols <- brewer.pal(8, "Dark2")[1:8]
p <- ggplot(subset(all_cl,!is.na(type))) + geom_boxplot(aes(x=as.character(type),y=K2P,color=as.character(type)),outlier.size = 0.5)+ theme_classic(12)+facet_grid(~Organ) + ylab("K2P value")+ theme(axis.text.x=element_text(angle=0,hjust=0.5)) + scale_color_manual(values=cols)+ylim(c(0,1)) + xlab("Evolution cluster")
ggsave("Ks_in_evolution.pdf",width = 7.5,height = 4)

data=subset(all_cl,!is.na(type))
for (j in c("AG","TE","OV")){
  cat(j,"\n")
  x0 = subset(data,type=="0"&Organ==j)
  for(i in 0:6){
    if (i == 0){next}
    x1 = subset(data,type==i&Organ==j)
    p = wilcox.test(x0$K2P,x1$K2P)$p.value
    cat(i,":",ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","NS")))," ",sep="")
  }
  cat("\n")
}

##################
all = read.table("../ConsDivEvolution/all_coevolve_scores.txt")
cols1 <- brewer.pal(9, "Set1")[c(9,1:3)]
cols2=c("white",brewer.pal(9, "Reds")[c(3,5,7)])
names(cols2)=c("n.s.","*","**","***")

####### AG versus TE
agVste = apply(all[,c(1,4)],1,function(x){
  if(x[2]=="ambiguous"){return("ambiguous")}
  else if(as.numeric(x[1])>0){return("concordant")}
  else if(as.numeric(x[1])<0){return("divergent")}
  else if(as.numeric(x[1])==0){return(x[2])}
  return("opps")
})

agVste =  data.frame(type=agVste,ifAGPeak = ifelse(names(agVste)%in%agPeak[,12],"+","-"),ifTEPeak = ifelse(names(agVste)%in%tePeak[,12],"+","-"))
agVste$ifPeak = apply(agVste[,2:3],1,function(x){paste("AG",x[1],"TE",x[2],sep="")})

tab = table(agVste[,c(1,4)])
r = tab/rowSums(tab)
r_gg = melt(r)
r_gg$type = factor(r_gg$type,levels = c("ns","concordant","divergent","ambiguous"))
names(cols1)=c("AG-TE-","AG-TE+","AG+TE-","AG+TE+")
p <- ggplot(r_gg) + geom_bar(aes(x=type,y=value,fill=ifPeak),stat="identity",color="white")+ theme_bw(12) + theme(axis.text.x=element_text(angle=45,hjust=1), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="right", legend.title=element_blank()) + labs(title="",y="Proportion of peaks", x="Coevolution state") +scale_fill_manual(values=cols1)
ggsave("AGvsTE.peak.distribution.pdf",width = 3.1,height = 3.5)


s0 = tab["ns",]
p_vals=matrix("n.s.",nrow=ncol(tab),ncol=nrow(tab),dimnames=list(colnames(tab),c("ns","concordant","divergent","ambiguous")))
for (i in rownames(tab)) {
  if (i == "ns"){next}
  s1 = tab[i,]
  #print(i)
  #print(chisq.test(cbind(s0,s1))$p.value)
  for(j in 1:ncol(tab)){
    #  cat(colnames(x)[j],": ")
    m=cbind(c(s1[j],sum(s1[-j])),c(s0[j],sum(s0[-j])))
    p=fisher.test(m)$p.value
    p_vals[j,i]=ifelse(p<=0.001,"***",ifelse(p<=0.01,"**",ifelse(p<=0.05,"*","n.s.")))
  }
}

r=r[c("ns","concordant","divergent","ambiguous"),]
pdf("FisherTest.heatmap.AGvsTE.pdf",width=2.8,height = 2)
h=draw(Heatmap(p_vals, name = "Proportion", col = cols2, show_row_names = T, show_column_names = T,cluster_rows = F,cluster_columns=F,show_column_dend =F,show_row_dend =F,
               cell_fun = function(j,i, x, y, width, height, fill) {
                 grid.text(sprintf("%.4f", r[j,i]), x,y, gp = gpar(fontsize = 6))}))
dev.off()

####### AG versus OV
agVsov = apply(all[,c(2,5)],1,function(x){
  if(x[2]=="ambiguous"){return("ambiguous")}
  else if(as.numeric(x[1])>0){return("concordant")}
  else if(as.numeric(x[1])<0){return("divergent")}
  else if(as.numeric(x[1])==0){return(x[2])}
  return("opps")
})

agVsov =  data.frame(type=agVsov,ifAGPeak = ifelse(names(agVsov)%in%agPeak[,12],"+","-"),ifOVPeak = ifelse(names(agVsov)%in%ovPeak[,12],"+","-"))
agVsov$ifPeak = apply(agVsov[,2:3],1,function(x){paste("AG",x[1],"OV",x[2],sep="")})

tab = table(agVsov[,c(1,4)])
r = tab/rowSums(tab)
r_gg = melt(r)
r_gg$type = factor(r_gg$type,levels = c("ns","concordant","divergent","ambiguous"))
names(cols1)=c("AG-OV-","AG-OV+","AG+OV-","AG+OV+")
p <- ggplot(r_gg) + geom_bar(aes(x=type,y=value,fill=ifPeak),stat="identity",color="white")+ theme_bw(12) + theme(axis.text.x=element_text(angle=45,hjust=1), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="right", legend.title=element_blank()) + labs(title="",y="Proportion of peaks", x="Coevolution state") +scale_fill_manual(values=cols1)
ggsave("AGvsOV.peak.distribution.pdf",width = 3.1,height = 3.5)


s0 = tab["ns",]
p_vals=matrix("n.s.",nrow=ncol(tab),ncol=nrow(tab),dimnames=list(colnames(tab),c("ns","concordant","divergent","ambiguous")))
for (i in rownames(tab)) {
  if (i == "ns"){next}
  s1 = tab[i,]
  #print(i)
  #print(chisq.test(cbind(s0,s1))$p.value)
  for(j in 1:ncol(tab)){
    #  cat(colnames(x)[j],": ")
    m=cbind(c(s1[j],sum(s1[-j])),c(s0[j],sum(s0[-j])))
    p=fisher.test(m)$p.value
    p_vals[j,i]=ifelse(p<=0.001,"***",ifelse(p<=0.01,"**",ifelse(p<=0.05,"*","n.s.")))
  }
}

r=r[c("ns","concordant","divergent","ambiguous"),]
pdf("FisherTest.heatmap.AGvsOV.pdf",width=2.8,height = 2)
h=draw(Heatmap(p_vals, name = "Proportion", col = cols2, show_row_names = T, show_column_names = T,cluster_rows = F,cluster_columns=F,show_column_dend =F,show_row_dend =F,
               cell_fun = function(j,i, x, y, width, height, fill) {
                 grid.text(sprintf("%.4f", r[j,i]), x,y, gp = gpar(fontsize = 6))}))
dev.off()

####### TE versus OV
teVsov = apply(all[,c(3,6)],1,function(x){
  if(x[2]=="ambiguous"){return("ambiguous")}
  else if(as.numeric(x[1])>0){return("concordant")}
  else if(as.numeric(x[1])<0){return("divergent")}
  else if(as.numeric(x[1])==0){return(x[2])}
  return("opps")
})

teVsov =  data.frame(type=teVsov,ifTEPeak = ifelse(names(teVsov)%in%tePeak[,12],"+","-"),ifOVPeak = ifelse(names(teVsov)%in%ovPeak[,12],"+","-"))
teVsov$ifPeak = apply(teVsov[,2:3],1,function(x){paste("TE",x[1],"OV",x[2],sep="")})

tab = table(teVsov[,c(1,4)])
r = tab/rowSums(tab)
r_gg = melt(r)
r_gg$type = factor(r_gg$type,levels = c("ns","concordant","divergent","ambiguous"))
names(cols1)=c("TE-OV-","TE-OV+","TE+OV-","TE+OV+")
p <- ggplot(r_gg) + geom_bar(aes(x=type,y=value,fill=ifPeak),stat="identity",color="white")+ theme_bw(12) + theme(axis.text.x=element_text(angle=45,hjust=1), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="right", legend.title=element_blank()) + labs(title="",y="Proportion of peaks", x="Coevolution state") +scale_fill_manual(values=cols1)
ggsave("TEvsOV.peak.distribution.pdf",width = 3.1,height = 3.5)


s0 = tab["ns",]
p_vals=matrix("n.s.",nrow=ncol(tab),ncol=nrow(tab),dimnames=list(colnames(tab),c("ns","concordant","divergent","ambiguous")))
for (i in rownames(tab)) {
  if (i == "ns"){next}
  s1 = tab[i,]
  #print(i)
  #print(chisq.test(cbind(s0,s1))$p.value)
  for(j in 1:ncol(tab)){
    #  cat(colnames(x)[j],": ")
    m=cbind(c(s1[j],sum(s1[-j])),c(s0[j],sum(s0[-j])))
    p=fisher.test(m)$p.value
    p_vals[j,i]=ifelse(p<=0.001,"***",ifelse(p<=0.01,"**",ifelse(p<=0.05,"*","n.s.")))
  }
}

r=r[c("ns","concordant","divergent","ambiguous"),]
pdf("FisherTest.heatmap.TEvsOV.pdf",width=2.8,height = 2)
h=draw(Heatmap(p_vals, name = "Proportion", col = cols2, show_row_names = T, show_column_names = T,cluster_rows = F,cluster_columns=F,show_column_dend =F,show_row_dend =F,
               cell_fun = function(j,i, x, y, width, height, fill) {
                 grid.text(sprintf("%.4f", r[j,i]), x,y, gp = gpar(fontsize = 6))}))
dev.off()

####### K2P
setwd("../Ks")
agKs = read.table("Accgland.narrowPeak.2rep.sorted.k2p",header = T)
teKs = read.table("Testis.narrowPeak.2rep.sorted.k2p",header = T)
ovKs = read.table("Ovary.narrowPeak.2rep.sorted.k2p",header = T)

ag_merged = merge(agPeak[,c(1,12)],agKs[,c(1,5)],by.x=1,by.y=1)
te_merged = merge(tePeak[,c(1,12)],teKs[,c(1,5)],by.x=1,by.y=1)
ov_merged = merge(ovPeak[,c(1,12)],ovKs[,c(1,5)],by.x=1,by.y=1)

x = rbind(data.frame(state=agVste[ag_merged[,2],1],K2P=ag_merged[,3],peak="AG",compare="AGvsTE"),
          data.frame(state=agVste[te_merged[,2],1],K2P=te_merged[,3],peak="TE",compare="AGvsTE"),
          data.frame(state=agVsov[ag_merged[,2],1],K2P=ag_merged[,3],peak="AG",compare="AGvsOV"),
          data.frame(state=agVsov[ov_merged[,2],1],K2P=ov_merged[,3],peak="OV",compare="AGvsOV"),
          data.frame(state=teVsov[te_merged[,2],1],K2P=te_merged[,3],peak="TE",compare="TEvsOV"),
          data.frame(state=teVsov[ov_merged[,2],1],K2P=ov_merged[,3],peak="OV",compare="TEvsOV"))
data = subset(x,!is.na(state))
data$state=factor(data$state,levels=c("ns","concordant","divergent","ambiguous"))
data$compare=factor(data$compare,levels=c("AGvsTE","AGvsOV","TEvsOV"))
p <- ggplot(data) + geom_boxplot(aes(x=state,y=K2P,color=state)) + facet_grid(peak~compare) + theme_bw(12) + theme(axis.text.x=element_text(angle=45,hjust=1), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="right", legend.title=element_blank()) + labs(title="",y="Proportion of peaks", x="Coevolution state") +scale_color_manual(values=c("grey",brewer.pal(9, "Set1")[c(1:4)]))
ggsave("K2P.distribution.pdf",width = 6.5,height = 3.5)

by(data[,1:3],data[,4],function(x){
  by(x[,c(1,2)],x[,3],function(y){
    #browser()
    #wilcox.test(y[,1]~y[,2])
    t0 = subset(y,state=="ns")[,2]
    t1 = subset(y,state=="concordant")[,2]
    t2 = subset(y, state=="divergent")[,2]
    t3 = subset(y, state=="ambiguous")[,2]
    p=c(wilcox.test(t0,t1)$p.value,wilcox.test(t0,t2)$p.value,wilcox.test(t0,t3)$p.value)
    return(p)
  })
})



#######peak distribution
peaks = subset(read.table("../coevolution/all_organ.narrowPeak.genomic.region",sep = "\t",quote = "",header = T,colClasses = "character"),as.numeric(Distance.to.TSS)>=-3000&as.numeric(Distance.to.TSS)<=100)
ix=tapply(peaks[,c(1,10)],peaks[,12],function(x){
  x[which.min(abs(as.numeric(x[,2]))),1]
})
peaks_nearest=peaks[peaks[,1]%in%ix,]
rownames(peaks_nearest)=peaks_nearest$Entrez.ID
table(peaks_nearest[,6])
##############################

#####################
fbgn2sym = read.table("../ConsDivEvolution/gene_id_2_gene_symbol.tab",row.names = 2)
rcl.s = read.table("../ConsDivEvolution/coevolve_state_by_species.txt",header=F)
rcl.s$symbol=fbgn2sym[rcl.s[,1],1]
rownames(rcl.s)=rcl.s[,1]

rcl.s$Peak.Score = peaks_nearest[rcl.s[,1],"Peak.Score"]
rcl.s[is.na(rcl.s$Peak.Score),"Peak.Score"]="000"

rcl.s$type = ifelse(rcl.s$Peak.Score=="111","AG+TE+OV+",
                      ifelse(rcl.s$Peak.Score=="110","AG+TE+OV-",
                             ifelse(rcl.s$Peak.Score=="101","AG+TE-OV+",
                                    ifelse(rcl.s$Peak.Score=="011","AG-TE+OV+",
                                           ifelse(rcl.s$Peak.Score=="100","AG+TE-OV-",
                                                  ifelse(rcl.s$Peak.Score=="010","AG-TE+OV-",
                                                         ifelse(rcl.s$Peak.Score=="001","AG-TE-OV+","AG-TE-OV-"
                                                         )))))))
states=names(sort(table(rcl.s$V2),decreasing = T))[1:11]
rcl.s$state = ifelse(rcl.s$V2%in%states,rcl.s$V2,"Others")
x=table(rcl.s[,c(6,5)])
r=x/rowSums(x)
x_gg = melt(r)
x_gg$state = factor(x_gg$state,levels=c(states,"Others"))
cols <- rev(brewer.pal(9, "Set1"))[1:8]
names(cols)=c("AG-TE-OV-","AG-TE-OV+","AG-TE+OV-","AG-TE+OV+","AG+TE-OV-","AG+TE-OV+","AG+TE+OV-","AG+TE+OV+")
p <- ggplot(x_gg) + geom_bar(aes(x=state,y=value,fill=as.character(type)),stat="identity",color="white")+ theme_bw(12) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="right", legend.title=element_blank()) + labs(title="",y="Proportion of peaks", x="Coevolution state") +scale_fill_manual(values=cols)
ggsave("peakInCoevolveCluster.distribution.pdf",width=6,height=3.4)


s0 = x["000",]
p_vals=matrix("n.s.",nrow=ncol(x),ncol=nrow(x),dimnames=list(colnames(x),c(states,"Others")))
for (i in c(states,"Others")) {
  if (i == "000"){next}
  s1 = x[i,]
  #print(i)
  #print(chisq.test(cbind(s0,s1))$p.value)
  for(j in 1:ncol(x)){
  #  cat(colnames(x)[j],": ")
    m=cbind(c(s1[j],sum(s1[-j])),c(s0[j],sum(s0[-j])))
    p=fisher.test(m)$p.value
    p_vals[j,i]=ifelse(p<=0.001,"***",ifelse(p<=0.01,"**",ifelse(p<=0.05,"*","n.s.")))
  }
}

cols=c("white",brewer.pal(9, "Reds")[c(3,5,7)])
names(cols)=c("n.s.","*","**","***")
r=r[c(states,"Others"),]
pdf("FisherTest.heatmap.pdf",width=5.8,height = 2.8)
h=draw(Heatmap(p_vals, name = "Proportion", col = cols, show_row_names = T, show_column_names = T,cluster_rows = F,cluster_columns=F,show_column_dend =F,show_row_dend =F,
               cell_fun = function(j,i, x, y, width, height, fill) {
                 grid.text(sprintf("%.4f", r[j,i]), x,y, gp = gpar(fontsize = 6.5))}))
dev.off()
#################################################################
## coevolution clusters within organs
k2p = read.table("../Ks/all.narrowPeak.merged.sorted.k2p",header=T,row.names = 1)
peaks_nearest$K2P = k2p[peaks_nearest[,1],4]




#################
te_cl = read.table("../coevolution/TE_cluster.txt",header = F)
te_peak = merge(te_cl,peaks_nearest[,c(6,12,20)],by.x=1,by.y=2,all.x=T)
te_peak[is.na(te_peak[,3]),3]="000"
# te_peak$type = ifelse(te_peak$Peak.Score=="111","AG+TE+OV+",
#                       ifelse(te_peak$Peak.Score=="110","AG+TE+OV-",
#                              ifelse(te_peak$Peak.Score=="101","AG+TE-OV+",
#                                   ifelse(te_peak$Peak.Score=="011","AG-TE+OV+",
#                                          ifelse(te_peak$Peak.Score=="100","AG+TE-OV-",
#                                                 ifelse(te_peak$Peak.Score=="010","AG-TE+OV-",
#                                                        ifelse(te_peak$Peak.Score=="001","AG-TE-OV+","AG-TE-OV-"
#                              )))))))

te_peak$type = ifelse(te_peak$Peak.Score%in%c("111","110","011","010"),"+","-")

ag_cl = read.table("../coevolution/AG_cluster.txt",header = F)
ag_peak = merge(ag_cl,peaks_nearest[,c(6,12,20)],by.x=1,by.y=2,all.x=T)
ag_peak[is.na(ag_peak[,3]),3]="000"
ag_peak$type = ifelse(ag_peak$Peak.Score%in%c("111","110","100","101"),"+","-")

ov_cl = read.table("../coevolution/OV_cluster.txt",header = F)
ov_peak = merge(ov_cl,peaks_nearest[,c(6,12,20)],by.x=1,by.y=2,all.x=T)
ov_peak[is.na(ov_peak[,3]),3]="000"
ov_peak$type = ifelse(ov_peak$Peak.Score%in%c("111","001","011","101"),"+","-")


all_peaks = rbind(data.frame(te_peak,organ="TE"),data.frame(ag_peak,organ="AG"),data.frame(ov_peak,organ="OV"))

p <- ggplot() + geom_boxplot(data= subset(all_peaks,as.character(Peak.Score)!="000"&K2P!=0&type!="-"),aes(x=as.character(V2),y=K2P,color=organ),position = position_dodge(0.9)) + theme_classic(12)  + ylab("K2P value")+scale_color_brewer(palette="Dark2")
ggsave("K2P_within_organs.pdf",width=8,height = 3)
data= subset(all_peaks,as.character(Peak.Score)!="000"&K2P!=0&type!="-")
for(i in c("AG","OV","TE")){
  x1 = 
}



#####################################################

te_dis = table(te_peak[,c(2,4)])
te_dis = rbind(te_dis,Whole=colSums(te_dis[,1:2]))
te_dis = te_dis/rowSums(te_dis)
te_dis_gg = melt(te_dis)
te_dis_gg$X1=factor(te_dis_gg$X1,levels=c("Whole",1:6))
p <- ggplot(te_dis_gg) + geom_bar(aes(x=X1,y=value,fill=as.character(X2)),stat="identity",color="black")+ theme_bw(10) + theme(axis.text.x=element_text(angle=0, hjust = 0.5), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="right", legend.title=element_blank()) + labs(title="",y="Proportion of peaks", x="Coevolution clusters in testis") +scale_fill_manual(values=c("white","black"))
ggsave("peakInTECluster.distribution.pdf",width=3.5,height=3)


te_cl = read.table("../coevolution/OV_cluster.txt",header = F)
te_peak = merge(te_cl,peaks_nearest[,c(6,12)],by.x=1,by.y=2,all.x=T)
te_peak[is.na(te_peak[,3]),3]="000"
# te_peak$type = ifelse(te_peak$Peak.Score=="111","AG+TE+OV+",
#                       ifelse(te_peak$Peak.Score=="110","AG+TE+OV-",
#                              ifelse(te_peak$Peak.Score=="101","AG+TE-OV+",
#                                     ifelse(te_peak$Peak.Score=="011","AG-TE+OV+",
#                                            ifelse(te_peak$Peak.Score=="100","AG+TE-OV-",
#                                                   ifelse(te_peak$Peak.Score=="010","AG-TE+OV-",
#                                                          ifelse(te_peak$Peak.Score=="001","AG-TE-OV+","AG-TE-OV-"
#                                                          )))))))
te_peak$type = ifelse(te_peak$Peak.Score%in%c("111","001","011","101"),"OV+","OV-")

te_dis = table(te_peak[,c(2,4)])
te_dis = rbind(te_dis,Whole=colSums(te_dis[,1:2]))
te_dis = te_dis/rowSums(te_dis)
te_dis_gg = melt(te_dis)
te_dis_gg$X1=factor(te_dis_gg$X1,levels=c("Whole",1:6))
p <- ggplot(te_dis_gg) + geom_bar(aes(x=X1,y=value,fill=X2),stat="identity",color="black")+ theme_bw(10) + theme(axis.text.x=element_text(angle=0, hjust = 0.5), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="right", legend.title=element_blank()) + labs(title="",y="Proportion of peaks", x="Coevolution clusters in ovary") +scale_fill_manual(values=c("white","black"))
ggsave("peakInOVCluster.distribution.pdf",width=3.5,height=3)

te_cl = read.table("../coevolution/AG_cluster.txt",header = F)
te_peak = merge(te_cl,peaks_nearest[,c(6,12)],by.x=1,by.y=2,all.x=T)
te_peak[is.na(te_peak[,3]),3]="000"
# te_peak$type = ifelse(te_peak$Peak.Score=="111","AG+TE+OV+",
#                       ifelse(te_peak$Peak.Score=="110","AG+TE+OV-",
#                              ifelse(te_peak$Peak.Score=="101","AG+TE-OV+",
#                                     ifelse(te_peak$Peak.Score=="011","AG-TE+OV+",
#                                            ifelse(te_peak$Peak.Score=="100","AG+TE-OV-",
#                                                   ifelse(te_peak$Peak.Score=="010","AG-TE+OV-",
#                                                          ifelse(te_peak$Peak.Score=="001","AG-TE-OV+","AG-TE-OV-"
#                                                          )))))))
te_peak$type = ifelse(te_peak$Peak.Score%in%c("111","110","100","101"),"AG+","AG-")

te_dis = table(te_peak[,c(2,4)])
te_dis = rbind(te_dis,Whole=colSums(te_dis[,1:2]))
te_dis = te_dis/rowSums(te_dis)
te_dis_gg = melt(te_dis)
p <- ggplot(te_dis_gg) + geom_bar(aes(x=X1,y=value,fill=as.character(X2)),stat="identity",color="black")+ theme_bw(10) + theme(axis.text.x=element_text(angle=0, hjust = 0.5), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="right", legend.title=element_blank()) + labs(title="",y="Proportion of peaks", x="Coevolution clusters in Accessory gland") +scale_fill_manual(values=c("white","black"))
ggsave("peakInAGCluster.distribution.pdf",width=3.5,height=3)


########################################



