setwd("/Users/aimeidai/Library/CloudStorage/OneDrive-个人/10.shangrui/CNE")
library("reshape")
library("ggplot2")
library("ggtree")
library(txdbmaker)
require(RColorBrewer)
options(stringsAsFactors = F)
t = read.tree("dm6.27way.commonNames.nh")
st = read.tree("dm6.5way.commonNames.nwk")
p <- ggtree(t,size=0.7)+ geom_tiplab(aes(label=gsub("rosophila_",". ",label)),offset = 0.05, font="Italic")+ theme_tree2() + ggplot2::xlim(0, 42)
ggsave("5way_phylogenetic_tree.pdf",width = 3.4,height = 3)

######## CNEs
# genomic distribution
#txdb=makeTxDbFromGFF("canonical_isoforms2.gtf")
library(ChIPseeker)
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
files=list(G2="G2",G3="G3",G4="G4",Merged="whole.CNE.bedops.merge.bed")
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
pdf("CNEGenomicDistribution.pdf",width = 6,height = 2)
plotAnnoBar(peakAnnoList, title="Genomic distribution of CNEs")
plotDistToTSS(peakAnnoList,title="Distribution of CNEs relative to TSS")
dev.off()

##############
#G1 = read.table("G1")
G2 = read.table("G2")
G3 = read.table("G3")
G4 = read.table("G4")
#G5 = read.table("G5")
all = read.table("whole.CNE.bedops.merge.bed")

#all[all[,1]%in%G1[,1]&all[,2]%in%G1[,2]&all[,3]%in%G1[,3],4]=1
all[all[,1]%in%G2[,1]&all[,2]%in%G2[,2]&all[,3]%in%G2[,3],4]=1
all[all[,1]%in%G3[,1]&all[,2]%in%G3[,2]&all[,3]%in%G3[,3],5]=1
all[all[,1]%in%G4[,1]&all[,2]%in%G4[,2]&all[,3]%in%G4[,3],6]=1
#all[all[,1]%in%G5[,1]&all[,2]%in%G5[,2]&all[,3]%in%G5[,3],8]=1
all[which(is.na(all),arr.ind = T)]=0
all$cons = apply(all[,4:6],1,function(x){
  paste(x,collapse ="")})
sort(table(all$cons))

all$consLevel = ifelse(all$cons == "111","High",ifelse(all$cons %in% c("011","110","101"),"Median","Low"))
write.table(all[,c(1:3,7:8)],"merged.CNE.cons.bed",sep="\t",quote = F,row.names=F,col.names = F)

# #### length distribution
# all$length = all[,3]-all[,2]+1
# all$consLevel=factor(all$consLevel,levels=c("High","Median","Low"))
# p <- ggplot(all,aes(x=length)) + geom_histogram(aes(y=..density..,color=consLevel),fill="white",position="identity", alpha=0.25,bins = 800)+ scale_color_brewer(palette="Dark2") + theme_classic(10) + coord_cartesian(xlim=c(0,1000))+theme(legend.position=c(0.75,0.7))
# ggsave("CNE_length_distribution.pdf",width = 3.5,height = 3)

##### ATAC-seq peaks overlapping CNEs
# count = rbind(AG=c(13300, 10686, 10074, 7561),TE=c(10858, 9356, 8854, 6737),OV=c(12139, 10343, 9783, 7546))
# pp = count[,2:4]/count[,1]
# colnames(pp)=paste("G",seq(2,4),sep="")
# gg = melt(pp)
# 
# #cols <- brewer.pal(11, "grey")[c(8,9,11)]
# cols=c("#D0E1F9", "#74A9CF", "#0570B0")
# p <- ggplot(data=gg) + geom_bar(aes(x=X1,y=value,fill = X2),color="white",lwd = 0.98,width = 0.8,position=position_dodge(0.8),stat="identity")+theme_bw(10)+scale_fill_manual(values=cols) + ylab("Proportion of peaks\noverlapping CNE") + xlab("")
# ggsave("Proportion_peaks_CNE.pdf",width = 4,height = 2)

##############CNEs that overlap with ATAC-seq peaks
# ag = read.table("CNE.Accgland.ol.bed")
# te = read.table("CNE.Testis.ol.bed")
# ov = read.table("CNE.Ovary.ol.bed")
# ############### CNE usage
# #c = c(AG=nrow(ag),TE=nrow(te),OV=nrow(ov))/335623
# 
# #################
# tot = unlist(table(all$consLevel))
# c = cbind(AG=table(ag[,4]),TE=table(te[,4]),OV=table(ov[,4]),total=tot)
# #prop = melt(sweep(c,2,colSums(c),"/"))
# prop = melt(c[,1:3]/257639)
# prop$X1=factor(prop$X1,levels = c("High","Median","Low"))
# p <- ggplot(data=prop) + geom_bar(aes(x=X2,y=value,fill = X1),color="white",lwd = 0.98,width = 0.8,stat="identity")+theme_bw(10)+scale_fill_manual(values=rev(cols)) + ylab("Proportion of CNEs\noverlapping ATAC peaks") + xlab("")
# ggsave("Proportion_ages_CNE.pdf",width = 2.5,height = 2)

########### assign CNEs to genes
# length distribution
agCNE = subset(read.table("Accgland.CNE.ol.bed"),V7>=-2000&V7<=2000)
teCNE = subset(read.table("Testis.CNE.ol.bed"),V7>=-2000&V7<=2000)
ovCNE = subset(read.table("Ovary.CNE.ol.bed"),V7>=-2000&V7<=2000)

lengths = rbind(data.frame(value=unique(agCNE[,c(1:4,7)])[,"V7"],tissue="AG"),data.frame(value=unique(teCNE[,c(1:4,7)])[,"V7"],tissue="TE"),data.frame(value=unique(ovCNE[,c(1:4,7)])[,"V7"],tissue="OV"))

p <- ggplot(lengths)+geom_histogram(aes(x=value,color = tissue,y=..density..),fill="white",bins = 200,alpha=0.4)+ scale_color_brewer(palette="Dark2") + theme_classic(10) + coord_cartesian(xlim=c(-2000,2000))+theme(legend.position=c(0.75,0.7)) + xlab("Distance to nearest TSS")
ggsave("peak_distance_distribution.pdf",width = 3.5,height = 3)

####### CNE regulations in genes
all = read.table("merged.CNE.cons.bed")
# CNE densities
agCNE = subset(read.table("Accgland.CNE.ol.bed"),V9!=-1)
teCNE = subset(read.table("Testis.CNE.ol.bed"),V9!=-1)
ovCNE = subset(read.table("Ovary.CNE.ol.bed"),V9!=-1)
##### conservation spectrum
cs=rbind(data.frame(melt(table(all[,"V4"])),tissue="All"),data.frame(melt(table(agCNE[,"V11"])),tissue="AG"),data.frame(melt(table(teCNE[,"V11"])),tissue="TE"),data.frame(melt(table(ovCNE[,"V11"])),tissue="OV"))
tot = c(All=nrow(all),AG=nrow(agCNE),TE=nrow(teCNE),OV=nrow(ovCNE))
cs$Fraction = 0
cs[cs[,3]=="All","Fraction"]=cs[cs[,3]=="All","value"]/tot["All"]
cs[cs[,3]=="AG","Fraction"]=cs[cs[,3]=="AG","value"]/tot["AG"]
cs[cs[,3]=="TE","Fraction"]=cs[cs[,3]=="TE","value"]/tot["TE"]
cs[cs[,3]=="OV","Fraction"]=cs[cs[,3]=="OV","value"]/tot["OV"]
# cs$Fraction=0
# cs[1:3,"Fraction"]=cs[1:3,2]/sum(cs[1:3,2])
# cs[4:6,"Fraction"]=cs[4:6,2]/sum(cs[4:6,2])
# cs[7:9,"Fraction"]=cs[7:9,2]/sum(cs[7:9,2])
cs$tissue = factor(cs$tissue, levels=c("All","AG","TE","OV"))
#cs$Var.1 = factor(cs$Var.1, levels=c("High","Median","Low"))
cs$Var = ifelse(cs$Var.1==1,"001",ifelse(cs$Var.1==10,"010",ifelse(cs$Var.1==11,"011",cs$Var.1)))
cols=brewer.pal(9, "Dark2")
p <- ggplot() + geom_bar(data=cs,aes(x=tissue,y=Fraction,fill=as.vector(Var)),stat="identity",color="white")+ scale_fill_manual(values=rev(cols)) + theme_classic(12)
ggsave("CNEConservationSpectrum.pdf",width = 3.5,height = 3.5)

for (i in unique(cs[,5])){
  cat(i,": \n")
  x01 = subset(cs,tissue=="All"&Var==i)[,2]
  x02 = sum(subset(cs,tissue=="All"&Var!=i)[,2])
  for (j in c("AG","TE","OV")) {
    cat(j,": ")
    xi1 = subset(cs,tissue==j&Var==i)[,2]
    xi2 = sum(subset(cs,tissue==j&Var!=i)[,2])
    m=cbind(c(xi2,xi1),c(x02,x01))
    #print(m)
    cat(wilcox.test(m)$p.value," ")
  }
  cat("\n")
}

# CNE densities
agCNE = subset(read.table("Accgland.CNE.ol.bed"),V7>=-1500&V7<=100&V9!=-1)
teCNE = subset(read.table("Testis.CNE.ol.bed"),V7>=-1500&V7<=100&V9!=-1)
ovCNE = subset(read.table("Ovary.CNE.ol.bed"),V7>=-1500&V7<=100&V9!=-1)
### CNE effects on expression divergence
# fc_all = read.table("../coevolution/logFC_dmVSdp.txt")
# 
# fc_all$ifAG = ifelse(rownames(fc_all) %in% agCNE[,6],"+","-")
# fc_all$ifTE = ifelse(rownames(fc_all) %in% teCNE[,6],"+","-")
# fc_all$ifOV = ifelse(rownames(fc_all) %in% ovCNE[,6],"+","-")
# 
# fc_gg = rbind(data.frame(log2FC = fc_all[,1],ifPeak=fc_all[,4],Organ="AG"),data.frame(log2FC = fc_all[,2],ifPeak=fc_all[,5],Organ="TE"),data.frame(log2FC = fc_all[,3],ifPeak=fc_all[,6],Organ="OV"))
# fc_gg$Organ = factor(fc_gg$Organ,levels=c("AG","TE","OV"))
# p<-ggplot(fc_gg) + geom_boxplot(aes(x=Organ,y=log2FC,colour =interaction(ifPeak,Organ)),outlier.size = 0.5)+ theme_bw(12) + ylab("log2FC (dme/dps)")+ theme(axis.text.x=element_text(angle=45,hjust=1)) + scale_color_manual(values=brewer.pal(9, "Paired")[c(5,6,3,4,1,2)])+coord_cartesian(ylim = c(-10,20))
# 
# ggsave("CNEEffect_in_evolution.pdf",width = 5,height = 3.5)
# tapply(fc_gg[,1],fc_gg[,2:3],function(x){median(x,na.rm=T)})
# # Organ
# # ifPeak         AG        TE            OV
# # - -0.4284080 0.2219565 -0.0007592768
# # + -0.2612377 0.3631382  0.0786466684
# tapply(fc_gg[,1:2],fc_gg[,3],function(x){wilcox.test(x[,1]~x[,2])$p.value})
# AG           TE           OV 
# 4.228862e-07 1.133707e-03 9.464313e-05 
agCount = tapply(agCNE[,8:9],agCNE[,4],nrow)
teCount = tapply(teCNE[,8:9],teCNE[,4],nrow)
ovCount = tapply(ovCNE[,8:9],ovCNE[,4],nrow)

agtot = unique(agCNE[,1:7])
agtot$Count = agCount[agtot[,4]]
tetot = unique(teCNE[,1:7])
tetot$Count = teCount[tetot[,4]]
ovtot = unique(ovCNE[,1:7])
ovtot$Count = ovCount[ovtot[,4]]

totCount = rbind(data.frame(agtot,tissue="AG"),data.frame(tetot,tissue="TE"),data.frame(ovtot,tissue="OV"))
totCount$peakLen = totCount[,3]-totCount[,2]+1

dens = tapply(totCount[,c("Count","peakLen")],totCount[,c("V6","tissue")],function(x){
  #browser()
  if(nrow(x)==0){return(NA)}
  return(sum(x[,1])/sum(x[,2])*1000)
})
##### CNE density on expression divergence
gs = intersect(rownames(dens),rownames(fc_all))
dens_gg = rbind(data.frame(CNEdensity=dens[gs,"AG"],log2FC=fc_all[gs,"AG"],tissue="AG"),data.frame(CNEdensity=dens[gs,"OV"],log2FC=fc_all[gs,"OV"],tissue="OV"),data.frame(CNEdensity=dens[gs,"TE"],log2FC=fc_all[gs,"TE"],tissue="TE"))
# p <- ggplot(dens_gg,aes(x=CNEdensity,y=log2FC)) + geom_point(size=0.5)+facet_wrap(~tissue) + theme_bw(10)+ geom_smooth(method=lm)
# ggsave("CNEdensityvsExpressionDivergence.pdf",width = 8,height = 3.5)
# 
# by(dens_gg[,1:2],dens_gg[,3],function(x){
#   #browser()
#   x=subset(x,!is.na(CNEdensity)&!is.infinite(log2FC))
#   summary(lm(x[,2]~x[,1]))
# })
# ##### CNE counts on expression evolution
# agNum = by(agCNE[,8:9],agCNE[,6],nrow)
# teNum = by(teCNE[,8:9],teCNE[,6],nrow)
# ovNum = by(ovCNE[,8:9],ovCNE[,6],nrow)
# 
# num_gg = rbind(data.frame(count=agNum,log2FC=fc_all[names(agNum),"AG"], tissue = "AG"),data.frame(count=teNum,log2FC=fc_all[names(teNum),"TE"], tissue = "TE"),data.frame(count=ovNum,log2FC=fc_all[names(ovNum),"OV"], tissue = "OV"))
# p <- ggplot(num_gg,aes(x=count,y=log2FC)) + geom_point(size=0.5)+facet_wrap(~tissue) + theme_bw(10)+ geom_smooth(method=lm)
# ggsave("CNECountvsExpressionDivergence.pdf",width = 8,height = 3.5)

######## CNE density along with coevolution clusters
cne_dens = data.frame()
#TE
tec1 = read.table("../coevolution/TE_cluster.txt",sep="\t")
tec2 = read.table("../coevolution/DEG_TE_dmdp.txt",sep="\t")
tec = merge(tec2,tec1,by.x=1,by.y=1)
rownames(tec)=tec[,1]
gs = intersect(tec[,1],rownames(dens))
cne_dens = rbind(cne_dens,data.frame(tec[gs,c(1,2,8)],CNEdensity=dens[gs,1],tissue="TE"))
cne_dens = rbind(cne_dens,data.frame(tec[gs,c(1,2)],V2.y=0,CNEdensity=dens[gs,2],tissue="TE"))

#OV
ovc1 = read.table("../coevolution/OV_cluster.txt",sep="\t")
ovc2 = read.table("../coevolution/DEG_OV_dmdp.txt",sep="\t")
ovc = merge(ovc2,ovc1,by.x=1,by.y=1)
rownames(ovc)=ovc[,1]
gs = intersect(ovc[,1],rownames(dens))
cne_dens = rbind(cne_dens,data.frame(ovc[gs,c(1,2,8)],CNEdensity=dens[gs,2],tissue="OV"))
cne_dens = rbind(cne_dens,data.frame(ovc[gs,c(1,2)],V2.y=0,CNEdensity=dens[gs,2],tissue="OV"))

#AG
agc1 = read.table("../coevolution/AG_cluster.txt",sep="\t")
agc2 = read.table("../coevolution/DEG_AG_dmdp.txt",sep="\t")
agc = merge(agc2,agc1,by.x=1,by.y=1)
rownames(agc)=agc[,1]
gs = intersect(agc[,1],rownames(dens))
cne_dens = rbind(cne_dens,data.frame(agc[gs,c(1,2,8)],CNEdensity=dens[gs,2],tissue="AG"))
cne_dens = rbind(cne_dens,data.frame(agc[gs,c(1,2)],V2.y=0,CNEdensity=dens[gs,2],tissue="AG"))
cne_dens$tissue = factor(cne_dens$tissue,levels=c("AG","TE","OV"))
cne_mean = melt(tapply(cne_dens[,4],cne_dens[,c(3,5)],function(x){median(x,na.rm=T)}))
cne_sr = melt(tapply(cne_dens[,4],cne_dens[,c(3,5)],function(x){sd(x,na.rm=T)}))
cne_ms = merge(cne_mean,cne_sr,by.x=c(1,2),by.y=c(1,2))
cne_ms$tissue = factor(cne_ms$tissue,levels=c("AG","TE","OV"))
p <- ggplot()+ geom_boxplot(data=cne_dens,aes(x=as.character(V2.y),y=CNEdensity,color=tissue),position = position_dodge(0.55),width=0.6,outlier.shape = NA) + geom_line(data=cne_ms,aes(x=as.character(V2.y),y=value.x,group=tissue,color=tissue),position = position_dodge(0.55),linetype="dashed") + geom_point(data = cne_ms,aes(x=as.character(V2.y),y=value.x,color=tissue),position = position_dodge(0.55))+ theme_bw(12) + ylab("CNE densities")+xlab("Coevolution clusters")+scale_color_manual(values = brewer.pal(3, "Set1"))+ coord_cartesian(ylim=c(0,10.5))
#+ geom_pointrange(aes(x=as.character(V2.y),y=value.x,ymin = value.x-value.y,ymax=value.x+value.y,color=tissue),position = position_dodge(0.55)) 
ggsave("CNEdensityAlongWithClusters.pdf",width = 4.5,height = 2.8)
do.call("rbind",by(cne_dens[,c(3:4)],cne_dens[,5],function(x){
  by(x[,2],x[,1],function(y){median(y,na.rm=T)})
}))
for(i in c("AG","TE","OV")){
  x0 = subset(cne_dens,tissue==i&V2.y==0)$CNEdensity
  print(i)
  for(j in seq(1,6)){
    xj = subset(cne_dens,tissue==i&V2.y==j)$CNEdensity
    print(wilcox.test(x0,xj)$p.value)
  }
}

for(i in seq(0,6)){
  cat(i,": ")
  xag = subset(cne_dens,tissue=="AG"&V2.y==i)$CNEdensity
  xte = subset(cne_dens,tissue=="TE"&V2.y==i)$CNEdensity
  xov = subset(cne_dens,tissue=="OV"&V2.y==i)$CNEdensity
  cat(wilcox.test(xag,xte)$p.value,wilcox.test(xag,xov)$p.value,wilcox.test(xov,xte)$p.value,"\n")
}
# 0 : 0.73584 0.8590473 0.6074316 
# 1 : 0.7738474 0.8725058 0.6148703 
# 2 : 0.4537304 0.9849126 0.3713921 
# 3 : 5.534691e-05 0.8411319 2.200142e-05 
# 4 : 6.311698e-05 0.5524237 0.0007225867 
# 5 : 7.118124e-05 0.4356416 4.947996e-06 
# 6 : 0.157584 0.3685042 0.7769277 

##### phastCons score on coevolutionay cluster
pCs = read.table("merged.CNE.cons.phastCons.bed",stringsAsFactors = F)
pCs$avgPhastcons = pCs$V5/pCs$V6
agpCs = merge(agCNE[,6:10],pCs[,c(1:4,7)],by.x=3:5,by.y=1:3)
ovpCs = merge(ovCNE[,6:10],pCs[,c(1:4,7)],by.x=3:5,by.y=1:3)
tepCs = merge(teCNE[,6:10],pCs[,c(1:4,7)],by.x=3:5,by.y=1:3)

###########
# agpCs$cluster = agc[agpCs$V6,"V2.y"]
# ovpCs$cluster = ovc[ovpCs$V6,"V2.y"]
# tepCs$cluster = tec[tepCs$V6,"V2.y"]
# 
# allpCs = rbind(data.frame(agpCs[,c("V6","V4","cluster")],tissue="AG"),data.frame(tepCs[,c("V6","V4","cluster")],tissue="TE"),data.frame(ovpCs[,c("V6","V4","cluster")],tissue="OV"),data.frame(agpCs[,c("V6","V4")],cluster=0,tissue="AG"),data.frame(tepCs[,c("V6","V4")],cluster=0,tissue="TE"),data.frame(ovpCs[,c("V6","V4")],cluster=0,tissue="OV"))
# allpCs_gg=melt(table(allpCs[,c("V4","cluster","tissue")]))
# allpCs_gg$Fraction = NA
# allpCs_gg[allpCs_gg[,3]=="AG","Fraction"] = allpCs_gg[allpCs_gg[,3]=="AG","value"]/nrow(subset(agpCs,!is.na(cluster)))
# allpCs_gg[allpCs_gg[,3]=="TE","Fraction"] = allpCs_gg[allpCs_gg[,3]=="TE","value"]/nrow(subset(tepCs,!is.na(cluster)))
# allpCs_gg[allpCs_gg[,3]=="OV","Fraction"] = allpCs_gg[allpCs_gg[,3]=="OV","value"]/nrow(subset(ovpCs,!is.na(cluster)))
# p <- ggplot() + geom_bar(data=subset(allpCs_gg,!is.na(cluster)&cluster!=0),aes(x=as.character(cluster),y=Fraction,color=V4),fill="white",stat = "identity",position = position_dodge()) + facet_grid(tissue~.) + theme_classic(12)
# ggsave("conservationSpectrumVSCluster.pdf",width = 4.5,height = 3.5)
###############
cal_bootstrapping <- function(x,times=1000){
  cvs=rep(NA,times)
  for(i in 1:times){
    xsub = sample(x,length(x)*0.8)
    cvs[i] = sd(xsub,na.rm=T)/mean(xsub,na.rm=T)
  }
  return(cvs)
}

CVs=data.frame()
for(i in unique(agc1[,2])){
  gl = subset(agc1,V2==i)[,1]
  values = agpCs[agpCs[,4]%in%gl,7]
  CVs = rbind(CVs,data.frame(CV=cal_bootstrapping(values),cluster=i,tissue="AG"))
}
for(i in unique(tec1[,2])){
  gl = subset(tec1,V2==i)[,1]
  values = tepCs[tepCs[,4]%in%gl,7]
  CVs = rbind(CVs,data.frame(CV=cal_bootstrapping(values),cluster=i,tissue="TE"))
}
for(i in unique(ovc1[,2])){
  gl = subset(ovc1,V2==i)[,1]
  values = ovpCs[ovpCs[,4]%in%gl,7]
  CVs = rbind(CVs,data.frame(CV=cal_bootstrapping(values),cluster=i,tissue="OV"))
}
gl = agc1[,1]
values = agpCs[agpCs[,4]%in%gl,7]
CVs = rbind(CVs,data.frame(CV=cal_bootstrapping(values),cluster=0,tissue="AG"))
gl = tec1[,1]
values = tepCs[tepCs[,4]%in%gl,7]
CVs = rbind(CVs,data.frame(CV=cal_bootstrapping(values),cluster=0,tissue="TE"))
gl = ovc1[,1]
values = ovpCs[ovpCs[,4]%in%gl,7]
CVs = rbind(CVs,data.frame(CV=cal_bootstrapping(values),cluster=0,tissue="OV"))
CVs$tissue = factor(CVs$tissue,levels=c("AG","TE","OV"))
cv_ms=melt(do.call("rbind",tapply(CVs[,1:2],CVs[,3],function(x){
  tapply(x[,1],x[,2],function(y){return(median(y))})
})))
colnames(cv_ms)=c("tissue","cluster","median")
p <- ggplot() + geom_line(data=cv_ms,aes(x=as.character(cluster),y=median,color=tissue,group = tissue),linewidth=0.5,linetype="dashed",position=position_dodge(0.5)) + geom_boxplot(data=CVs,aes(x=as.character(cluster),y=CV,color=tissue),width=0.6,position = position_dodge(0.5),outlier.shape = NA) + geom_point(data=cv_ms,aes(x=as.character(cluster),y=median,color=tissue),position = position_dodge(0.5)) + theme_bw(12) + ylab("Coefficient of variation\n(CNE phastCons)")+xlab("Coevolution clusters")+scale_color_manual(values = brewer.pal(3, "Set1"))
ggsave("cvphastConsAlongWithClusters.pdf",width = 4.5,height = 2.8)

do.call("rbind",by(CVs[,1:2],CVs[,3],function(x){
  by(x[,1],x[,2],median)
}))

# 0         1         2         3         4         5         6
# AG 0.6971407 0.5468067 0.8080219 0.7103633 0.7035636 0.6775834 0.7699147
# TE 0.7024809 0.5368545 0.6372400 0.7099841 0.7457639 0.7063721 0.7153385
# OV 0.6643147 0.4692403 0.5929064 0.6511566 0.6787369 0.7151097 0.6671280

for(i in c("AG","TE","OV")){
  print(i)
  x0 = subset(CVs,cluster==0&tissue==i)[,1]
  for(j in 1:6){
    xi = subset(CVs,cluster==j&tissue==i)[,2]
    cat(j,":",wilcox.test(x0,xi)$p.value,"\t")
  }
  cat("\n")
}

# x=tapply(agpCs[,7],agpCs[,4],function(x){
#   #if(length(x)==1){return(0)}
#   return(sd(x,na.rm=T)/mean(x,na.rm=T))})
# agc1$pC_cv = x[agc1[,1]]
# 
# x=tapply(tepCs[,7],tepCs[,4],function(x){
#   #if(length(x)==1){return(0)}
#   return(sd(x,na.rm=T)/mean(x,na.rm=T))})
# tec1$pC_cv = x[tec1[,1]]
# 
# x=tapply(ovpCs[,7],ovpCs[,4],function(x){
#   #if(length(x)==1){return(0)}
#   return(sd(x,na.rm=T)/mean(x,na.rm=T))})
# ovc1$pC_cv = x[ovc1[,1]]
# 
# pCs_cv=rbind(data.frame(agc1,tissue="AG"),
#   data.frame(tec1,tissue="TE"),
#   data.frame(ovc1,tissue="OV"),data.frame(agc1[,c(1,3)],V2=0,tissue="AG"),data.frame(tec1[,c(1,3)],V2=0,tissue="TE"),data.frame(ovc1[,c(1,3)],V2=0,tissue="OV"))
# 
# 
# cv_mean = merge(melt(tapply(pCs_cv[,3],pCs_cv[,c(2,4)],function(x){mean(x,na.rm=T)})),melt(tapply(pCs_cv[,3],pCs_cv[,c(2,4)],function(x){sd(x,na.rm=T)})),by.x=c(1,2),by.y=c(1,2))
# 
# cv_mean$tissue = factor(cv_mean$tissue,levels=c("AG","TE","OV"))
# p <- ggplot(data=cv_mean) + geom_line(aes(x=as.character(V2),y=value.x,group=tissue),color="grey",linewidth=0.4,position = position_dodge(0.55)) + geom_pointrange(aes(x=as.character(V2),y=value.x,ymin = value.x-value.y,ymax=value.x+value.y,color=tissue),position = position_dodge(0.55)) + theme_bw(12) + ylab("Coefficient of variation\n(CNE phastCons)")+xlab("Coevolution clusters")+scale_color_manual(values = brewer.pal(3, "Set1")) + coord_cartesian(ylim=c(0,1.5))
# ggsave("cvphastConsAlongWithClusters.pdf",width = 4.5,height = 2.8)

for(i in seq(0,6)){
  cat(i,": ")
  xag = subset(pCs_cv,tissue=="AG"&V2==i)[,3]
  xov = subset(pCs_cv,tissue=="OV"&V2==i)[,3]
  xte = subset(pCs_cv, tissue=="TE"&V2==i)[,3]
  cat(wilcox.test(xag,xte)$p.value,wilcox.test(xag,xov)$p.value,wilcox.test(xov,xte)$p.value,"\n")
}
# 0 : 0.99077 0.4865127 0.4945701 
# 1 : 0.6174142 0.290689 0.4061528 
# 2 : 0.6857328 0.8228454 0.8906401 
# 3 : 0.8200191 0.585136 0.7106219 
# 4 : 0.353863 0.9838108 0.3633914 
# 5 : 0.7616644 0.6885012 0.6115786 
# 6 : 0.4361657 0.1575458 0.5105043 

# for(i in c("AG","TE","OV")){
#   cat(i,": ")
#   x0 = subset(pCs_cv, tissue==i&V2==0)
#   for (j in seq(1,6)){
#     xj = subset(pCs_cv, tissue==i&V2==j)
#     cat(wilcox.test(x0[,3],xj[,3])$p.value," ")
#   }
#   cat("\n")
# }
# AG : 0.05008722  0.6783017  0.6101588  0.4351318  0.4829192  0.04298121  
# TE : 0.1602471  0.7894156  0.2420592  0.7226197  0.3442203  0.337125  
# OV : 0.7792986  0.9591621  0.2803865  0.8452779  0.4520362  0.8825559  
########## 
# Sfps, TSGs, OSGs
sfp = read.table("../Sfp_SBG/wigby2020.txt")
tsg = read.table("../Sfp_SBG/testis_specific_genes.txt")
osg = read.table("../Sfp_SBG/ovary_specific_genes.txt")

dens_gg = subset(melt(dens),!is.na(value))
dens_gg=rbind(data.frame(dens_gg,type="All"),data.frame(dens_gg[dens_gg[,1]%in%sfp[,3]&dens_gg[,2]=="AG",],type="TSG"),data.frame(dens_gg[dens_gg[,1]%in%tsg[,1]&dens_gg[,2]=="TE",],type="TSG"),data.frame(dens_gg[dens_gg[,1]%in%osg[,1]&dens_gg[,2]=="OV",],type="TSG"))
dens_gg$tissue = factor(dens_gg$tissue,levels=c("AG","TE","OV"))
cols <- brewer.pal(12,"Set1")[c(1,2,3)]
p <- ggplot(dens_gg) + geom_boxplot(aes(x=type,y=value,fill = tissue)) + theme_classic(10) + scale_fill_manual(values = cols) + ylab("CNE density")
ggsave("CNE_density.pdf",width = 3,height = 3)

for(i in c("All","TSG")){
  cat(i, ": ")
  xag = subset(dens_gg,type==i&tissue=="AG")$value
  xov = subset(dens_gg,type==i&tissue=="OV")$value
  xte = subset(dens_gg,type==i&tissue=="TE")$value
  cat(median(xag,na.rm=T),median(xov,na.rm=T),median(xte,na.rm=T),"\n",wilcox.test(xag,xte)$p.value,wilcox.test(xag,xov)$p.value,wilcox.test(xte,xov)$p.value,"\n")
}

for (i in c("AG","TE","OV")){
  cat(i, ": ")
  xall = subset(dens_gg,type=="All"&tissue==i)$value
  xtsg = subset(dens_gg, type="TSG"&tissue==i)$value
  cat(wilcox.test(xall,xtsg)$p.value,"\n")
}

## heatmap of CNE densities
dens[which(is.na(dens),arr.ind=T)]=0
col_fun = colorRamp2(c(0, 2, 5), c("white", "pink", "red"))
pdf("whole_genome_CNE_densities.pdf",width = 3,height = 5)
ht_ag=draw(Heatmap(as.matrix(dens),name = "CNE densities", km =4, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=T,cluster_rows = T,show_column_dend = F,col = col_fun))
dev.off()
dens=dens[,c(1,3,2)]
ix = ifelse(dens>=4,1,0)
types=melt.list(apply(ix,1,function(x){paste(x,collapse = "")}))
types$ifTSG = ifelse(types$L1%in%sfp[,3],"Sfps","Other")
types$ifTSG = ifelse(types$L1%in%tsg[,1],"TSGs",types$ifTSG)
types$ifTSG = ifelse(types$L1%in%osg[,1],"OSGs",types$ifTSG)

prop = data.frame(sweep(table(types[,c(1,3)]),1,table(types[,1]),"/"))
prop$value = factor(prop$value,levels=c("000","100","010","001","110","101","011","111"))
prop$ifTSG = factor(prop$ifTSG, levels=c("Sfps","TSGs","OSGs"))
cols <- brewer.pal(9, "Set1")
p <- ggplot(subset(prop,ifTSG!="Other")) + geom_bar(aes(x=value,y=Freq,fill = ifTSG),stat = "identity")+ theme_classic(10) + scale_fill_manual(values=cols)+xlab("High CNE density (AG|TE|OV)") + ylab("Proportion of Tissue specific genes")
ggsave("distribution_TSG_in_CNE_density.pdf",width = 3.5,height = 4)
##### conservation spectrum
cs=rbind(data.frame(numbers=table(agCNE$V11),tissue="AG",type="All"),data.frame(numbers=table(ovCNE$V11),tissue="OV",type="All"),data.frame(numbers=table(teCNE$V11),tissue="TE",type="All"),data.frame(numbers=table(agCNE[agCNE$V6%in%sfp[,3],"V11"]),tissue="AG",type="TSG"),data.frame(numbers=table(teCNE[teCNE$V6%in%tsg[,1],"V11"]),tissue="TE",type="TSG"),data.frame(numbers=table(ovCNE[ovCNE$V6%in%osg[,1],"V11"]),tissue="OV",type="TSG"))
cs_p = tapply(cs[,1:2],cs[,3:4],function(x){
  #browser()
  p = x[,2]/sum(x[,2])
  names(p) = x[,1]
  return(p)
})
cs_p=rbind(data.frame(melt(do.call("rbind",cs_p[,"All"])),type="All"),data.frame(melt(do.call("rbind",cs_p[,"TSG"])),type="TSG"))
cs_p$X2 = factor(cs_p$X2, levels=rev(c("High","Median","Low")))
cols <- brewer.pal(9, "Set1")[c(9,1)]
p <- ggplot(cs_p) + geom_bar(aes(x=X2,y=value,fill=type),color="white",stat="identity",position = position_dodge()) + facet_wrap(~X1) + theme_classic(10) + scale_fill_manual(values=cols)+xlab("Conservation level") + ylab("Proportion of CNEs")
ggsave("conservation_spectrum.pdf",width=5,height = 3.5)
#######
pCs = read.table("merged.CNE.cons.phastCons.bed",stringsAsFactors = F)
pCs$avgPhastcons = pCs$V5/pCs$V6
agpCs = merge(agCNE[,6:10],pCs[,c(1:4,7)],by.x=3:5,by.y=1:3)
ovpCs = merge(ovCNE[,6:10],pCs[,c(1:4,7)],by.x=3:5,by.y=1:3)
tepCs = merge(teCNE[,6:10],pCs[,c(1:4,7)],by.x=3:5,by.y=1:3)

allpCs = rbind(data.frame(agpCs,tissue="AG",type="All"),data.frame(tepCs,tissue="TE",type="All"),data.frame(ovpCs,tissue="OV",type="All"),data.frame(agpCs[agpCs[,4]%in%sfp[,3],],tissue="AG",type="TSG"),data.frame(tepCs[tepCs[,4]%in%tsg[,1],],tissue="TE",type="TSG"),data.frame(ovpCs[ovpCs[,4]%in%osg[,1],],tissue="OV",type="TSG"))
allpCs$tissue=factor(allpCs$tissue, levels=c("AG","TE","OV"))
cols <- brewer.pal(9, "Set1")
# p <- ggplot(allpCs) + geom_violin(aes(x=type,y=avgPhastcons,color=tissue))+theme_classic(10) + xlab("") + ylab("Average phastCons score")+scale_color_manual(values=cols)
# ggsave("avgPhastcons_boxplot.pdf",width=4.5,height = 3)
p <- ggplot(allpCs) + geom_density(aes(x=avgPhastcons,color=tissue),linewidth=0.8)+facet_wrap(~type)+theme_classic(10)+scale_color_manual(values=cols)+xlab("Averaged phastCons")
ggsave("avgPhastcons_density_plot.pdf",width = 4,height = 3)

#####
ranges = melt(tapply(allpCs$avgPhastcons,allpCs$V6,function(x){
  #browser()
  return(sd(x)/mean(x))
}))
ranges$ifTSG = ifelse(ranges[,1]%in%sfp[,3],"Sfps","Other")
ranges$ifTSG = ifelse(ranges[,1]%in%tsg[,1],"TSGs",types$ifTSG)
ranges$ifTSG = ifelse(ranges[,1]%in%osg[,1],"OSGs",types$ifTSG)
ranges$ifTSG = factor(ranges$ifTSG,levels=c("Sfps","TSGs","OSGs","Other"))
p <- ggplot(subset(ranges,value!=0)) + geom_boxplot(aes(x=ifTSG,y=value,color=ifTSG))+theme_classic(10) + xlab("") + ylab("Coefficient of variance of phastCons")+scale_color_manual(values=c(cols[1:3],"grey"))
ggsave("CVofPhastcons_boxplot.pdf",width=2.8,height = 3)

x0 = subset(ranges,ifTSG=="Other")[,2]
for(i in c("Sfps","TSGs","OSGs")){
  xi = subset(ranges,ifTSG==i)[,2]
  cat(wilcox.test(x0,xi)$p.value," ")
}
