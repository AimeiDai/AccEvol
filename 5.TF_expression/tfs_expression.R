setwd("/Users/aimeidai/Library/CloudStorage/OneDrive-个人/10.shangrui/remove_batch_effect")
library(sva)
library(pamr)
library(limma)
library(ggplot2)
library(reshape)
library(edgeR)
library(ape)
library(circlize)
options(stringsAsFactors=F)
require(RColorBrewer)
library(ComplexHeatmap)
#############
##############
sample_table = read.table("sample_table_rc.txt",sep="\t",header=T)
rownames(sample_table) = sample_table[,3]

count_matrix = as.matrix(read.table("raw_counts_5_species.txt",header=T, sep="\t"))
colnames(count_matrix) = sample_table[,3]


 setwd("../TF_expression")
count_tfs = as.matrix(read.table("raw_counts_tfs_5_species.txt",header=T, sep="\t"))
colnames(count_tfs) = sample_table[,3]
count_matrix = rbind(count_matrix,count_tfs)

tmp_s = subset(sample_table,batch!=10)
tmp = count_matrix[,tmp_s[,3]]
adjusted <- ComBat_seq(tmp, batch=tmp_s$batch, group=NULL, covar_mod = tmp_s[,c("merged_tissue2","species")])

x = subset(sample_table,batch==10)
adj_count=cbind(adjusted, count_matrix[,x[,3]])

y <- DGEList(counts=adj_count, group=sample_table$merged_SRR)
keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
cpm_adj = cpm(y)

fbgn2sym = read.table("../Sfp_SBG/gene_id_2_gene_symbol.tab")
rownames(fbgn2sym)=fbgn2sym[,2]
gene_len = read.table("../Sfp_SBG/canonical_exon.length")
merged = merge(fbgn2sym,gene_len,by.x= 1,by.y=1)
rownames(merged)=merged[,2]

x = intersect(rownames(adj_count),merged[,2])
rpk = adj_count[x,]/merged[x,3]*1000
tpm = sweep(rpk,2,colSums(rpk),"/")*10^6
rownames(tpm) = merged[rownames(tpm),1]

st = subset(sample_table,merged_tissue%in%c("AG","TE","OV"))
tpm_m = do.call("cbind",tapply(st[,3],st[,c(4,13)],function(x){
  #browser()
  if(length(x)>1){
    y = as.matrix(apply(tpm[,x],1,mean))
  }else{
    y = as.matrix(tpm[,x])
  }
  colnames(y)=paste(unique(st[st[,3]%in%x,c(4,13)]),collapse =".");
  return(y)
}))
#############

plot_heatmap <- function(tpm_sfp,ix,pdf_file="SFP.Heatmap.pdf",k=5,w=4,h=4){
  x1 = tpm_sfp[,c("dme.AG","dsi.AG","dya.AG","dan.AG","dps.AG")]
  x2 = tpm_sfp[,c("dme.TE","dsi.TE","dya.TE","dan.TE","dps.TE")]
  x3 = tpm_sfp[,c("dme.OV","dsi.OV","dya.OV","dan.OV","dps.OV")]
  colnames(x1)=c("dme","dsi","dya","dan","dps")
  colnames(x2)=c("dme","dsi","dya","dan","dps")
  colnames(x3)=c("dme","dsi","dya","dan","dps")
  s1 = Heatmap(log10(x1+1),name = "log10(TPM)", km =k, show_row_names = T, show_column_names = TRUE,cluster_columns=F,col = col_fun,use_raster=F,column_names_gp = grid::gpar(fontsize = 8),row_names_gp = grid::gpar(fontsize = 6))
  s2 = Heatmap(log10(x2+1),name = "log10(TPM)", km =k, show_row_names = T, show_column_names = TRUE,cluster_columns=F,col = col_fun,use_raster=F,column_names_gp = grid::gpar(fontsize = 8),row_names_gp = grid::gpar(fontsize = 6))
  s3 = Heatmap(log10(x3+1),name = "log10(TPM)", km =k, show_row_names = T, show_column_names = TRUE,cluster_columns=F,col = col_fun,use_raster=F,column_names_gp = grid::gpar(fontsize = 8),row_names_gp = grid::gpar(fontsize = 6))
  ht_ix = Heatmap(ix,name = "Expression", km =k, show_row_names = T, show_column_names = TRUE,cluster_columns=F,col = col2,use_raster=F,column_names_gp = grid::gpar(fontsize = 8),row_names_gp = grid::gpar(fontsize = 6))
  pdf(pdf_file,height = h,width=w)
  ht_ag=draw(s1+s2+s3+ht_ix)
  dev.off()
}

#tfs = read.table("TFs_for_TFBS.tab")
allmotifs = read.table("whole_TFs_for_TFBS.tab",header = F)
tfs=allmotifs[,2]
tpm_tfs = tpm_m[intersect(rownames(tpm_m),tfs),]
x=setdiff(tfs,rownames(tpm_m))
fbgn2sym[fbgn2sym[,1]%in%x,]
paste(fbgn2sym[fbgn2sym[,1]%in%x,2],collapse = "|")
col_fun = colorRamp2(c(0, 1, 2), c("white","pink", "red"))


# tpm_cut = cbind(rowMeans(tpm_tfs[,1:5]),rowMeans(tpm_tfs[,6:10]),rowMeans(tpm_tfs[,11:15]))
# colnames(tpm_cut) = c("AG","OV","TE")
tpm_cut=tpm_tfs[,c(2,12,7)]
ix1=apply(ifelse(tpm_cut>=5,1,0),1,function(x){paste(x,collapse="")})
ix=ifelse(ix1=="001","OVuniq",ifelse(ix1=="010","TEuniq",ifelse(ix1=="100","AGuniq",ifelse(ix1=="110","AG-TE",ifelse(ix1=="101","AG-OV",ifelse(ix1=="011","TE-OV",ifelse(ix1=="111","AG-TE-OV","Others")))))))
#ix=cbind(ifelse(tpm_cut[,1]>=5,ifelse(tpm_cut[,1]>=20,"high","median"),"low"),ifelse(tpm_cut[,3]>=5,ifelse(tpm_cut[,3]>=20,"high","median"),"low"),ifelse(tpm_cut[,2]>=5,ifelse(tpm_cut[,2]>=20,"high","median"),"low"))
#colnames(ix) = c("AG","TE","OV")
col2 = c(brewer.pal(9, "Paired")[1:8])
names(col2)=c("Others","AGuniq","TEuniq","OVuniq","AG-TE","AG-OV","TE-OV","AG-TE-OV")
plot_heatmap(tpm_tfs,ix,"TF.heatmap.pdf",k=5,w=4.5,h=8)

#apply(ix,2,function(x){sum(x!="low")})
#########
agtfs = rownames(tpm_cut[tpm_cut[,1]>=5,])
tetfs = rownames(tpm_cut[tpm_cut[,2]>=5,])
ovtfs = rownames(tpm_cut[tpm_cut[,3]>=5,])
tf_ol=venn(list(AG=agtfs,TE=tetfs,OV=ovtfs))
######### regulatory complexity
# TFBS density
# whole genome
agcountp = read.table("AG.TFBS.counts",sep="\t",header = T)
agcountp[,1]=paste("chr",agcountp[,1],sep = "")
tecountp = read.table("TE.TFBS.counts",sep="\t",header = T)
tecountp[,1]=paste("chr",tecountp[,1],sep = "")
ovcountp = read.table("OV.TFBS.counts",sep="\t",header = T)
ovcountp[,1]=paste("chr",ovcountp[,1],sep = "")
### generate motifs and TF symbols, modify TF symbols by manual
# allmotifs = unique(c(colnames(agcount)[4:ncol(agcount)],colnames(tecount)[4:ncol(tecount)],colnames(ovcount)[4:ncol(ovcount)]))
# tfs = gsub("MA\\d+\\.\\d\\.","",motifs)
# write.table(cbind(allmotifs,tfs),"motifs_2_tfs.txt",row.names = F, col.names = T,sep="\t",quote=F)

#allmotifs = read.table("motifs_2_tfs.txt",header = T)
allmotifs = allmotifs[allmotifs[,2]%in%names(ix[ix!="Others"]),]
agmotifs = intersect(colnames(agcountp),allmotifs[allmotifs[,2]%in%agtfs,3])#
agcount = agcountp[,c("CHR","START","END",agmotifs)]
temotifs = intersect(colnames(tecountp),allmotifs[allmotifs[,2]%in%tetfs,3])#
tecount = tecountp[,c("CHR","START","END",temotifs)]
ovmotifs = intersect(colnames(ovcountp),allmotifs[allmotifs[,2]%in%ovtfs,3])#
ovcount = ovcountp[,c("CHR","START","END",ovmotifs)]

# col_fun = colorRamp2(c(0, 1, 3), c("white", "pink", "red"))
# pdf("whole_genome_TFBS_counts.pdf",width = 5,height = 4)
# ht_ag=draw(Heatmap(as.matrix(agcount[,4:ncol(agcount)]),name = "TFBS counts", km =4, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=T,cluster_rows = T,show_column_dend = F,col = col_fun))
# ht_te=draw(Heatmap(as.matrix(tecount[,4:ncol(tecount)]),name = "TFBS counts", km =4, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=T,cluster_rows = T,show_column_dend = F,col = col_fun))
# ht_ov=draw(Heatmap(as.matrix(ovcount[,4:ncol(ovcount)]),name = "TFBS counts", km =4, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=T,cluster_rows = T,show_column_dend = F,col = col_fun))
# dev.off()

agPeak = subset(read.table("../ATAC-seq/Accgland.ATACpeaks.distance"),V7>=-1500&V7<=100)
tePeak = subset(read.table("../ATAC-seq/Testis.ATACpeaks.distance"),V7>=-1500&V7<=100)
ovPeak = subset(read.table("../ATAC-seq/Ovary.ATACpeaks.distance"),V7>=-1500&V7<=100)
agmotif = merge(agcount,agPeak[,c(1:4,6)],by.x=1:3,by.y=1:3,all.y=T)
temotif = merge(tecount,tePeak[,c(1:4,6)],by.x=1:3,by.y=1:3,all.y=T)
ovmotif = merge(ovcount,ovPeak[,c(1:4,6)],by.x=1:3,by.y=1:3,all.y=T)
genes = unique(c(agmotif$V6,temotif$V6,ovmotif$V6))

#############
agCbyG = matrix(NA,length(genes),length(agmotifs),dimnames = list(genes,agmotifs))
tmp=do.call("rbind",tapply(agmotif[,agmotifs],agmotif$V6,function(x){
  if(nrow(x)==1){return(x)}
  else{
    # browser()
    return(colSums(x,na.rm = T))
  }
}))
agCbyG[rownames(tmp),colnames(tmp)]=as.matrix(tmp)
agCbyG[which(is.na(agCbyG),arr.ind = T)]=0


teCbyG = matrix(NA,length(genes),length(temotifs),dimnames = list(genes,temotifs))
tmp=do.call("rbind",tapply(temotif[,temotifs],temotif$V6,function(x){
  if(nrow(x)==1){return(x)}
  else{
    # browser()
    return(colSums(x,na.rm = T))
  }
}))
teCbyG[rownames(tmp),colnames(tmp)]=as.matrix(tmp)
teCbyG[which(is.na(teCbyG),arr.ind = T)]=0

ovCbyG = matrix(NA,length(genes),length(ovmotifs),dimnames = list(genes,ovmotifs))
tmp=do.call("rbind",tapply(ovmotif[,ovmotifs],ovmotif$V6,function(x){
  if(nrow(x)==1){return(x)}
  else{
    # browser()
    return(colSums(x,na.rm = T))
  }
}))
ovCbyG[rownames(tmp),colnames(tmp)]=as.matrix(tmp)
ovCbyG[which(is.na(ovCbyG),arr.ind = T)]=0

library(gplots)
v=venn(list(AG=agmotifs,TE=temotifs,OV=ovmotifs))

commonmotifs = attributes(v)$intersections$`AG:TE:OV`
#commonmotifs = c("MA1700.1.Clamp","MA0015.1.Cf2","MA0529.2.BEAF.32","MA1456.1.Dref","MA1842.1.msl.1","MA0535.1.Mad","MA1459.1.M1BP","MA0533.1.su.Hw.","MA1839.1.GATAd","MA1460.1.pho","MA0531.1.CTCF")
agtemotifs = attributes(v)$intersections$`AG:TE`
agovmotifs = attributes(v)$intersections$`AG:OV`
teovmotifs = attributes(v)$intersections$`TE:OV`
aguniq = attributes(v)$intersections$`AG`
teuniq = attributes(v)$intersections$`TE`
ovuniq = attributes(v)$intersections$`OV`

#####
col_fun = colorRamp2(c(0, 2, 4), c("white", "pink", "red"))
heatmap_evolType <- function(genelist=genes,w=6,h=4,kms=4,filename="promoter_TFBS_counts.pdf"){
    c_ag=Heatmap(as.matrix(agCbyG[genelist,commonmotifs]),name = "TFBS counts", km =kms, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=F,cluster_rows = T,show_column_dend = F,show_row_dend = F,column_names_side = "top",col = col_fun,column_names_gp = grid::gpar(fontsize = 5))
    c_te=Heatmap(as.matrix(teCbyG[genelist,commonmotifs]),name = "TFBS counts", km =kms, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=F,column_names_side = "top",cluster_rows = F,show_column_dend = F,col = col_fun,column_names_gp = grid::gpar(fontsize = 5))
    c_ov=Heatmap(as.matrix(ovCbyG[genelist,commonmotifs]),name = "TFBS counts", km =kms, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=F,column_names_side = "top",cluster_rows = F,show_column_dend = F,col = col_fun,column_names_gp = grid::gpar(fontsize = 5))
    
    u_agov1=Heatmap(as.matrix(agCbyG[genelist,agovmotifs]),name = "TFBS counts", km =kms, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=F,column_names_side = "top",cluster_rows = F,show_column_dend = F,col = col_fun,column_names_gp = grid::gpar(fontsize = 5))
    u_agov2=Heatmap(as.matrix(ovCbyG[genelist,agovmotifs]),name = "TFBS counts", km =kms, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=F,column_names_side = "top",cluster_rows = F,show_column_dend = F,col = col_fun,column_names_gp = grid::gpar(fontsize = 5))
    #u_agte1=Heatmap(as.matrix(agCbyG[genelist,agtemotifs]),name = "TFBS counts", km =kms, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=F,column_names_side = "top",cluster_rows = F,show_column_dend = F,col = col_fun,column_names_gp = grid::gpar(fontsize = 5))
    #u_agte2=Heatmap(as.matrix(teCbyG[genelist,agtemotifs]),name = "TFBS counts", km =kms, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=F,column_names_side = "top",cluster_rows = F,show_column_dend = F,col = col_fun,column_names_gp = grid::gpar(fontsize = 5))
    u_teov1=Heatmap(as.matrix(teCbyG[genelist,teovmotifs]),name = "TFBS counts", km =kms, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=F,column_names_side = "top",cluster_rows = F,show_column_dend = F,col = col_fun,column_names_gp = grid::gpar(fontsize = 5))
    u_teov2=Heatmap(as.matrix(ovCbyG[genelist,teovmotifs]),name = "TFBS counts", km =kms, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=F,column_names_side = "top",cluster_rows = F,show_column_dend = F,col = col_fun,column_names_gp = grid::gpar(fontsize = 5))
    
    u_ag=Heatmap(as.matrix(agCbyG[genelist,aguniq]),name = "TFBS counts", km =kms, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=T,column_names_side = "top",cluster_rows = F,show_column_dend = F,col = col_fun,column_names_gp = grid::gpar(fontsize = 5))
    u_te=Heatmap(as.matrix(teCbyG[genelist,teuniq]),name = "TFBS counts", km =kms, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=T,column_names_side = "top",cluster_rows = F,show_column_dend = F,col = col_fun,column_names_gp = grid::gpar(fontsize = 5))
    u_ov=Heatmap(as.matrix(ovCbyG[genelist,ovuniq]),name = "TFBS counts", km =kms, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=T,column_names_side = "top",cluster_rows = F,show_column_dend = F,col = col_fun,column_names_gp = grid::gpar(fontsize = 5))
    pdf(filename,width = w,height = h)
    h=draw(c_ag+c_te+c_ov+u_agov1+u_agov2+u_teov1+u_teov2+u_ag+u_te+u_ov)#+u_agte1+u_agte2
    dev.off()
    return(c_ag)
}

#
###########
c_ag=heatmap_evolType(w=7.2,h=3,kms=4)
rcl.ag <- melt(row_order(c_ag))

#########
evolType = read.table("../ConsDivEvolution/coevolve_state_by_species.txt")
evolType$type =  ifelse(evolType$V2%in%c("100","010","001","110","101","011","111"),"Concordant",ifelse(evolType$V2%in%c("-100","0-10","00-1","-1-10","-10-1","0-1-1"),"Divergent",ifelse(evolType$V2%in%c("000"),"ambiguous","Other")))
rownames(evolType) = evolType$V1

agmotif$type = evolType[agmotif$V6,3]
agmotif$ifCM=ifelse(rowSums(agmotif[,commonmotifs])==0,"No","Yes")
table(agmotif[,c("type","ifCM")])
write.table(agmotif[,c("CHR","START","END","type","ifCM")],"AG_ACR_evolType_CommonMotif.bed",sep="\t",row.names = F, col.names = F,quote=F)

temotif$type = evolType[temotif$V6,3]
temotif$ifCM=ifelse(rowSums(temotif[,commonmotifs])==0,"No","Yes")
table(temotif[,c("type","ifCM")])
write.table(temotif[,c("CHR","START","END","type","ifCM")],"TE_ACR_evolType_CommonMotif.bed",sep="\t",row.names = F, col.names = F,quote=F)

ovmotif$type = evolType[ovmotif$V6,3]
ovmotif$ifCM=ifelse(rowSums(ovmotif[,commonmotifs])==0,"No","Yes")
table(ovmotif[,c("type","ifCM")])
write.table(ovmotif[,c("CHR","START","END","type","ifCM")],"OV_ACR_evolType_CommonMotif.bed",sep="\t",row.names = F, col.names = F,quote=F)

#######
agphast = read.table("AG_ACR_evolType_CommonMotif.bed-phastCons")
tephast = read.table("TE_ACR_evolType_CommonMotif.bed-phastCons")
ovphast = read.table("OV_ACR_evolType_CommonMotif.bed-phastCons")
acrphast = rbind(data.frame(agphast,organ="AG"),data.frame(tephast,organ="TE"),data.frame(ovphast,organ="OV"))
acrphast$avgPhast = acrphast$V6/acrphast$V7
acrphast$organ=factor(acrphast$organ,levels=c("AG","TE","OV"))
acrphast$type=ifelse(acrphast$V5=="No","withoutCommonMotifs",acrphast$V4)
p<-ggplot(subset(acrphast,!is.na(type)&type!="Other"&type!="withoutCommonMotifs")) + geom_boxplot(aes(x=type,y=avgPhast,color=organ)) + facet_grid(~organ)+theme_bw(10)+ theme(axis.text.x=element_text(angle=45,hjust=1),strip.background = element_blank(), axis.title.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+ scale_color_brewer(palette="Set1") + ylim(c(0,1.5))+ylab("Averaged phastCons")
ggsave("ACR_phastCons_withCommonMotifs.pdf",width = 3.2,height = 3)
by(acrphast$avgPhast,acrphast[,c("type","organ")],median)
paired = combn(c("ambiguous","Concordant","Divergent","withoutCommonMotifs"),2)
for(i in c("AG","TE","OV")){
  p=NULL
  for (j in 1:ncol(paired)){
    x1 = subset(acrphast,organ==i&type==paired[1,j])
    x2 = subset(acrphast,organ==i&type==paired[2,j])
   p=c(p,wilcox.test(x1$avgPhast,x2$avgPhast)$p.value)
  }
  #p = p.adjust(p)
  cat(ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","ns"))))
  cat("\n")
}
paired
# xag=do.call("cbind",by(agmotif[,commonmotifs],agmotif$type,function(x){
#   return(colSums(x>=1,na.rm=T))
# }))
# xag=rbind(xag,"total"=table(evolType$type))
# #sweep(xag,1,tot,"/")
# tot = table(evolType$type)
# rag=melt(xag[,1:3]/rowSums(xag[,1:3]))
# p <- ggplot(rag)+geom_point(aes(x=X1,y=value,color=X2),stat = "identity")+ theme_bw(10)+ scale_color_manual(values=c('#56B4E9', '#E69F00', "grey")) + theme(axis.text.x=element_text(angle=45,hjust=1),strip.background = element_blank(), axis.title.x = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) 


################ regulatory complexity
agmotif$width=agmotif[,3]-agmotif[,2]+1
ovmotif$width=ovmotif[,3]-ovmotif[,2]+1
temotif$width=temotif[,3]-temotif[,2]+1
agden=tapply(agmotif[,c(agmotifs,"width")],agmotif$V6,function(x){
  #browser()
  d = sum(x[,agmotifs],na.rm=T)/sum(x[,"width"])*1000
  return(d)
})
teden = tapply(temotif[,c(temotifs,"width")],temotif$V6,function(x){
  #browser()
  d = sum(x[,temotifs],na.rm=T)/sum(x[,"width"])*1000
  return(d)
})
ovden = tapply(ovmotif[,c(ovmotifs,"width")],ovmotif$V6,function(x){
  #browser()
  d = sum(x[,ovmotifs],na.rm=T)/sum(x[,"width"])*1000
  return(d)
})
# ############# expression divergence vs regulatory complexity
# fc_all = read.table("../coevolution/logFC_dmVSdp.txt")
# aggene = intersect(rownames(fc_all),names(agden))
# tegene = intersect(rownames(fc_all),names(teden))
# ovgene = intersect(rownames(fc_all),names(ovden))
# fc_den = rbind(data.frame(FC=fc_all[aggene,"AG"],TFBSdensity=agden[aggene],tissue="AG"), data.frame(FC=fc_all[ovgene,"OV"],TFBSdensity=ovden[ovgene], tissue = "OV"), data.frame(FC=fc_all[tegene,"TE"],TFBSdensity=teden[tegene],tissue="TE"))
# 
# p <- ggplot(fc_den) + geom_point(aes(x=FC,y=log2(TFBSdensity)),size=0.5) + facet_grid(~tissue) + theme_bw(10)
# ggsave("FCvsTFBSdensity.pdf",width = 9,height = 3)
# by(fc_den[,1:2],fc_den[,3],function(x){
#   #browser()
#   x=subset(x,abs(FC)!=Inf)
#   summary(lm(x[,1]~x[,2]))
# })
############## coevolution clusters versus regulatory complexity
# clusters = c(names(sort(table(evolType$V2),decreasing = T))[1:11],"Others")
# 
# gc = intersect(names(agden),evolType[,1])
# evolType[gc,"AGden"] = agden[gc]
# gc = intersect(names(teden),evolType[,1])
# evolType[gc,"TEden"] = teden[gc]
# gc = intersect(names(ovden),evolType[,1])
# evolType[gc,"OVden"] = unlist(ovden[gc])
# 
# evolType_gg = melt(evolType,id=c("V1","V2"))
# evolType_gg$type = factor(ifelse(evolType_gg$V2 %in% clusters,evolType_gg$V2,"Others"),levels = clusters)
# data=evolType_gg[evolType_gg[,4]!="NA",]
# p <- ggplot(evolType_gg) + geom_boxplot(aes(x=type,y=value,color=variable),outlier.shape = NA,position = position_dodge()) + theme_bw(10)+ theme(axis.text.x=element_text(angle=45,hjust=1)) +scale_color_manual(values = brewer.pal(3, "Set1")) + coord_cartesian(ylim = c(0,30))
# ggsave("EvolType_TFBS_density.pdf",width = 7.5,height = 3)

# tec1 = read.table("../coevolution/TE_cluster.txt",sep="\t")
# rownames(tec1)=tec1[,1]
# regu_complexity = data.frame(Complexity=teden,tec1[names(teden),],tissue="TE")
# regu_complexity = rbind(regu_complexity,data.frame(Complexity=teden,V1=tec1[names(teden),1],V2=0,tissue="TE"))
# 
# ovc1 = read.table("../coevolution/OV_cluster.txt",sep="\t")
# 
# rownames(ovc1)=ovc1[,1]
# regu_complexity = rbind(regu_complexity,data.frame(Complexity=ovden,ovc1[names(ovden),],tissue="OV"))
# regu_complexity = rbind(regu_complexity,data.frame(Complexity=ovden,V1=ovc1[names(ovden),1],V2=0,tissue="OV"))
# 
# agc1 = read.table("../coevolution/AG_cluster.txt",sep="\t")
# rownames(agc1)=agc1[,1]
# regu_complexity = rbind(regu_complexity,data.frame(Complexity=agden,agc1[names(agden),],tissue="AG"))
# regu_complexity = rbind(regu_complexity,data.frame(Complexity=agden,V1=agc1[names(agden),1],V2=0,tissue="AG"))
# 
# data=subset(regu_complexity,!is.na(V1))
# 
# data_summary <- function(x) {
#   y = boxplot(x,plot=F)
#   m <- y$stats[3]
#   ymin <- y$stats[2]
#   ymax <- y$stats[4]
#   return(c(y=m,ymin=ymin,ymax=ymax))
# }
# x=tapply(data$Complexity,data[,c("V2","tissue")],data_summary)
# x_gg = data.frame()
# for(i in rownames(x)){
#   y=do.call("rbind",x[i,])
#   x_gg = rbind(x_gg,data.frame(y,tissue=rownames(y),cluster=i))
# }
# p <- ggplot(x_gg,aes(x=cluster,color=tissue)) + geom_hline(yintercept =x_gg[1:3,1],color="darkgrey",linetype="dashed")+ geom_line(aes(y=y,group=tissue),position = position_dodge(0.5))+ geom_point(aes(y=y),position = position_dodge(0.5),size=3) + geom_errorbar(aes(ymin = ymin,ymax = ymax),position = position_dodge(0.5)) + theme_bw(10) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab("TFBS density")+xlab("Coevolution clusters")+scale_color_manual(values = brewer.pal(3, "Set1")) + coord_cartesian(ylim = c(2,16))
# #+ scale_y_log10()
# ggsave("TFBSDensityAlongWithClusters.pdf",width = 4,height = 2.5)
# 
# by(data[,c(1,4)],data[,3],function(x){
#   x1 = subset(x,tissue=="AG")[,1]
#   x2 = subset(x, tissue == "OV")[,1]
#   x3 = subset(x, tissue == "TE")[,1]
#   cat(median(x1),median(x2),median(x3),"\n")
#   cat("AGvsOV:",wilcox.test(x1,x2)$p.val, "AGvsTE:",wilcox.test(x1,x3)$p.val, "OVvsTE:",wilcox.test(x2,x3)$p.val,"\n")
# })
# 6.122449 5.988024 6.30063 
# AGvsOV: 0.00834233 AGvsTE: 0.5767507 OVvsTE: 0.001330382 
# 5.56381 5.555556 5.681818 
# AGvsOV: 0.8722405 AGvsTE: 0.8667625 OVvsTE: 0.9553244 
# 4.914934 4.938272 4.885395 
# AGvsOV: 0.6682532 AGvsTE: 0.9405838 OVvsTE: 0.5974104 
# 5.111603 5.295675 5.651056 
# AGvsOV: 0.6847052 AGvsTE: 0.06562203 OVvsTE: 0.01685811 
# 6.068882 5.714286 6.842626 
# AGvsOV: 0.01032333 AGvsTE: 0.03620578 OVvsTE: 1.661667e-06 
# 6.706908 6.750965 6.530651 
# AGvsOV: 0.5026274 AGvsTE: 0.5537736 OVvsTE: 0.9955336 
# 6.802721 6.802721 7.034007 
# AGvsOV: 0.5381259 AGvsTE: 0.3142687 OVvsTE: 0.7335277 
# 
# by(data[,c(1,3)],data[,4],function(x){
#   x0 = subset(x,V2=="0")[,1]
#   for(i in 1:6){
#     xi = subset(x,V2==i)[,1]
#     cat(i,wilcox.test(x0,xi)$p.val,"\n")
#   }
# })
#AG
# 1 0.1589891 
# 2 0.02760903 
# 3 0.0003140284 
# 4 0.9060753 
# 5 0.001214554 
# 6 0.4172453 
#OV
# 1 0.4824121 
# 2 0.001147545 
# 3 0.001891589 
# 4 0.2085254 
# 5 1.095722e-05 
# 6 0.009689389 
#TE
# 1 0.1556551 
# 2 6.120842e-05 
# 3 0.04965069 
# 4 0.02344919 
# 5 0.09256584 
# 6 0.06678759 
######## regulatory complexity of Acps, TSG and OSG
acps = read.table("../Sfp_SBG/wigby2020.txt")
tsg = read.table("../Sfp_SBG/testis_specific_genes.txt")
osg = read.table("../Sfp_SBG/ovary_specific_genes.txt")

###### overlaps with cluster1 in TFBS counts
cls = matrix(NA,ncol=4)
cls[1,1] = length(intersect(acps[,3],genes))
cls[1,2] = length(intersect(tsg[,1],genes))
cls[1,3] = length(intersect(osg[,1],genes))
cls[1,4] = length(genes)
ols = matrix(NA,nrow=4,ncol=4,dimnames = list(1:4,c("Sfps","TSG","OSG","Background")))
for(i in 1:4){
  ci = subset(rcl.ag,L1==i)[,1]
  cat(length(ci),"\t")
  ols[as.numeric(i),1]=length(intersect(acps[,3],genes[ci]))
  ols[as.numeric(i),2]=length(intersect(tsg[,1],genes[ci]))
  ols[as.numeric(i),3]=length(intersect(osg[,1],genes[ci]))
  ols[as.numeric(i),4]=length(genes[ci])
}

p_vals = matrix(NA,ncol=3,nrow=4,dimnames=list(1:4,c("Sfps","TSG","OSG")))
for(j in c("Sfps","TSG","OSG")){
  for(i in 1:4){
    x10 = colSums(ols[-i,c("Background",j)])
    x11 = ols[i,c("Background",j)]
    p_vals[i,j]=fisher.test(rbind(x10,x11),alternative = "greater")$p.value
  }
}

signif = ifelse(p_vals<=0.001,"***",ifelse(p_vals<=0.01,"**",ifelse(p_vals<=0.05,"*","ns")))

#colnames(ols)=c("Sfps","TSG","OSG","Background")
r=sweep(ols,2,cls,"/")
r_gg = melt(r)
cols=brewer.pal(4,"Dark2")
p <- ggplot(r_gg) + geom_bar(aes(x=X2,y=value,fill = as.character(X1)),stat="identity",color="white")+ theme_bw(12) + theme(axis.text.x=element_text(angle=45,hjust=1,vjust = 1), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="right", legend.title=element_blank()) + labs(title="",y="Proportion of genes", x="") +scale_fill_manual(values=rev(cols[1:4]))
ggsave("distributionOfTSGsAlongTFBSClusters.pdf",width=2.6,height = 3.5)

cols=c("white",brewer.pal(9, "Reds")[c(3,5,7)])
names(cols)=c("ns","*","**","***")
#r=r[c("Background",1:3),]
r=t(r)
pdf("FisherTest.heatmap.pdf",width=2.3,height = 2)
h=draw(Heatmap(cbind(signif,"Background"=rep("ns",4)), name = "Proportion", col = cols, show_row_names = T, show_column_names = T,cluster_rows = F,cluster_columns=F,show_column_dend =F,show_row_dend =F,
               cell_fun = function(j,i, x, y, width, height, fill) {
                 grid.text(sprintf("%.4f", r[j,i]), x,y, gp = gpar(fontsize = 6))}))
dev.off()

#################
#acps[acps[,3]%in%tsg[,1],]
# a1 = intersect(acps[,3],names(agden))
# a2 = intersect(tsg[,1],names(teden))
# a3 = intersect(osg[,1],names(ovden))
# organ_specific=rbind(data.frame(TFBSdensity=agden[a1],tissue="AG",type="OSG"),data.frame(TFBSdensity=teden[a2],tissue="TE",type="OSG"),data.frame(TFBSdensity=ovden[a3],tissue="OV",type="OSG"),data.frame(TFBSdensity=agden[setdiff(names(agden),a1)],tissue="AG",type="all"),data.frame(TFBSdensity=teden[setdiff(names(teden),a2)],tissue="TE",type="all"),data.frame(TFBSdensity=ovden[setdiff(names(ovden),a3)],tissue="OV",type="all"))
# organ_specific$tissue = factor(organ_specific$tissue,levels=c("AG","TE","OV"))
# 
# p <- ggplot() + geom_boxplot(data=organ_specific,aes(x=tissue,y=TFBSdensity,colour = tissue),outlier.shape = NA)+facet_wrap(~type) + theme_classic(10) + ylab("TFBS density") + coord_cartesian(ylim = c(0,25))+ theme(axis.text.x=element_text(angle=45,hjust=1)) +scale_color_manual(values = brewer.pal(6, "Set1"))
# ggsave("organ_specific_vs_TFBS_density.pdf",width = 3,height = 3)
# 
# by(organ_specific[,c(1,2)],organ_specific[,3],function(x){
#   #browser()
#   x1 = subset(x,tissue=="AG")[,1]
#   x2 = subset(x,tissue=="TE")[,1]
#   x3 = subset(x,tissue=="OV")[,1]
#   p=c(wilcox.test(x1,x2)$p.value,
#   wilcox.test(x1,x3)$p.value,
#   wilcox.test(x2,x3)$p.value)
# })
# organ_specific[, 3]: all
# [1] 0.130907372 0.200809920 0.004478541
# ---------------------------------------------------------------------- 
#   organ_specific[, 3]: OSG
# [1] 0.0168465738 0.0007350254 0.2307820698
#######
heatmap_evolType <- function(genelist=genes,w=6,h=4,filename="promoter_TFBS_counts.pdf"){
  
  c_te=Heatmap(as.matrix(teCbyG[genelist,commonmotifs]),name = "TFBS counts", km =3, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=F,cluster_rows = F,show_column_dend = F,col = col_fun,column_names_gp = grid::gpar(fontsize = 5))
  c_ov=Heatmap(as.matrix(ovCbyG[genelist,commonmotifs]),name = "TFBS counts", km =3, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=F,cluster_rows = F,show_column_dend = F,col = col_fun,column_names_gp = grid::gpar(fontsize = 5))
  
  u_ag=Heatmap(as.matrix(agCbyG[genelist,aguniq]),name = "TFBS counts", km =3, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=T,cluster_rows = F,show_column_dend = F,col = col_fun,column_names_gp = grid::gpar(fontsize = 5))
  u_te=Heatmap(as.matrix(teCbyG[genelist,teuniq]),name = "TFBS counts", km =3, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=T,cluster_rows = F,show_column_dend = F,col = col_fun,column_names_gp = grid::gpar(fontsize = 5))
  u_ov=Heatmap(as.matrix(ovCbyG[genelist,ovuniq]),name = "TFBS counts", km =3, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=T,cluster_rows = F,show_column_dend = F,col = col_fun,column_names_gp = grid::gpar(fontsize = 5))
  pdf(filename,width = w,height = h)
  h=draw(c_ag+c_te+c_ov+u_ag+u_te+u_ov)
  dev.off()
  return(c_ag)
}

x = agCbyG[intersect(acps[,3],agmotif$V6),agmotifs]

agdeg = read.table("../coevolution/DEG_AG_dmdp.txt")
labels = fbgn2sym[intersect(rownames(x),subset(agdeg,V7!="N.S.")[,1]),1]

rownames(x)=fbgn2sym[rownames(x),1]

#ha <- rowAnnotation(link = anno_mark(at = match(labels,rownames(x)), labels = labels), width = unit(0.01, "cm") + max_text_width(labels))
ha = rowAnnotation(text = anno_mark( at = which(rownames(x)%in%labels), labels =  labels, labels_gp=gpar(fontsize=5)))
c_ag=Heatmap(as.matrix(x),name = "TFBS counts", km =2,show_row_names = F, show_column_names = TRUE,cluster_columns=T,cluster_rows = T,show_column_dend = F,col = col_fun,column_names_gp = gpar(fontsize = 6))
pdf("Acps_TFBS.pdf",width = 4,height=5.5)
draw(c_ag+ha)
dev.off()

                 