library(ggplot2)
library(reshape)
library(ComplexHeatmap)
library(colorRamp2)
options(stringsAsFactors=F)
require(RColorBrewer)
col_fun = colorRamp2(c(0, 1.5, 3), c("blue", "white", "red"))
set.seed(1)

cols <- brewer.pal(12, "Set1")[1:6]
setwd("/Users/aimeidai/OneDrive/10.shangrui/remove_batch_effect")
sample_table = read.table("sample_table_rc.txt",sep="\t",header=T)
rownames(sample_table) = sample_table[,3]
rc_adj = read.table("adjusted_raw_count_5species.txt",header=T)
colnames(rc_adj)=sample_table[,3]

setwd("../Sfp_SBG")
fbgn2sym = read.table("gene_id_2_gene_symbol.tab")
rownames(fbgn2sym)=fbgn2sym[,2]
gene_len = read.table("canonical_exon.length")
merged = merge(fbgn2sym,gene_len,by.x= 1,by.y=1)
rownames(merged)=merged[,2]

x = intersect(rownames(rc_adj),merged[,2])
rpk = rc_adj[x,]/merged[x,3]*1000
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
###### sfps
sfp = read.table("wigby2020.txt")
x = intersect(sfp[,1],rownames(tpm))
tpm_sfp = tpm_m[x,]

ppi = read.table("physical_interactions.tab")
partners = setdiff(unlist(ppi[ppi[,1]%in%sfp[,3]|ppi[,2]%in%sfp[,3],]),sfp[,3])
x = intersect(merged[partners,1],rownames(tpm))
tpm_partners = tpm_m[x,]
###########
plot_heatmap <- function(tpm_sfp,pdf_file="SFP.Heatmap.pdf",k=5,w=4,h=4){
  x1 = tpm_sfp[,c("dme.AG","dsi.AG","dya.AG","dan.AG","dps.AG")]
  x2 = tpm_sfp[,c("dme.TE","dsi.TE","dya.TE","dan.TE","dps.TE")]
  x3 = tpm_sfp[,c("dme.OV","dsi.OV","dya.OV","dan.OV","dps.OV")]
  colnames(x1)=c("dme","dsi","dya","dan","dps")
  colnames(x2)=c("dme","dsi","dya","dan","dps")
  colnames(x3)=c("dme","dsi","dya","dan","dps")
  s1 = Heatmap(log10(x1+1),name = "log10(TPM)", km =k, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=F,col = col_fun,use_raster=F)
  s2 = Heatmap(log10(x2+1),name = "log10(TPM)", km =k, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=F,col = col_fun,use_raster=F)
  s3 = Heatmap(log10(x3+1),name = "log10(TPM)", km =k, show_row_names = FALSE, show_column_names = TRUE,cluster_columns=F,col = col_fun,use_raster=F)
  pdf(pdf_file,height = h,width=w)
  ht_ag=draw(s1+s2+s3)
  dev.off()
}
###############
plot_heatmap(tpm_sfp)
plot_heatmap(tpm_partners,"SFP.partner.Heatmap.pdf",w=4,h=5.5)

########
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

fc_all = read.table("../coevolution/logFC_dmVSdp.txt")

others = setdiff(setdiff(rownames(fc_all),sfp[,3]),partners)
fc_sfp = fc_all[intersect(rownames(fc_all),sfp[,3]),]
fc_partners = fc_all[intersect(rownames(fc_all),partners),]
fc_others = fc_all[others,]
fc_gg = rbind(data.frame(melt(fc_sfp),type="Sfps"),data.frame(melt(fc_partners),type="PPI partners"),data.frame(melt(fc_others),type="Other genes"))
fc_gg$variable = factor(fc_gg$variable,levels=c("AG","TE","OV"))
fc_gg$type = factor(fc_gg$type,levels=c("Sfps","PPI partners","Other genes"))
cols <- c("#808A87",brewer.pal(12,"Paired")[c(5,6)])
  
# p<-ggplot(data=fc_gg,aes(x=variable,y=value,color=type)) + geom_violin(width=0.8,position = position_dodge(0.8),trim = T) + stat_summary(fun.data=mean_sdl, mult=1,geom="point",size=0.4,position = position_dodge(0.8)) + theme_bw(10) + ylab("log2FC (dme/dps)")+ theme(axis.text.x=element_text(angle=45,hjust=1)) + scale_color_manual(values=rev(cols))+ coord_cartesian(ylim = c(-10,20))
# ggsave("Sfp_ppi_evolution2.pdf",width = 3.6,height = 3.2)
fc_gg$type = factor(fc_gg$type, levels=c("Other genes","PPI partners","Sfps"))
fc_medians = melt(tapply(fc_gg[,2],fc_gg[,c(1,3)],function(x){median(x,na.rm=T)}))
fc_medians$type=factor(fc_medians$type,levels=c("Other genes","PPI partners","Sfps"))
p<-ggplot()+ geom_vline(data=fc_medians,aes(xintercept=value,color=type),linetype="dashed",linewidth=0.3) +stat_ecdf(data=fc_gg,aes(x=value,color=type),linewidth=0.5) + facet_grid(variable~.) + theme_bw(10) + ylab("Density")+ theme(axis.text.x=element_text(angle=45,hjust=1), panel.grid=element_blank()) + scale_color_manual(values=cols) + xlim(c(-5,5)) + xlab("log2 FC (dme/dps)")
ggsave("Sfp_ppi_evolution_ecdf.pdf",width = 3,height = 3)

# ggplot(data=subset(fc_gg,!is.infinite(value)), aes(x = value, y = variable, fill = type,colour = type)) +geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, gradient_lwd = 1.) + scale_color_manual(values = cols)+
#   scale_x_continuous(expand = c(0.01, 0),limits = c(-8,8)) +
#   scale_y_discrete(expand = c(0.05, 0)) + 
#   scale_fill_viridis(name = "type", option = "C",discrete = T,alpha = 0.06) +labs(title = 'Expression changes',x="log2 FC (D.mel/D.pse)") +
#   theme_ridges(font_size = 10, grid = TRUE) + theme(axis.title.y = element_blank())
# ggsave("Sfps_logFC_3d.pdf",width = 5,height = 3)

by(fc_gg[,2:3],fc_gg[,1],function(x){
  x1 = subset(x,type=="Sfps")[,1]
  x2 = subset(x,type=="PPI partners")[,1]
  x3 = subset(x,type=="Other genes")[,1]
  p = c(wilcox.test(x1,x3)$p.value,wilcox.test(x2,x3)$p.value,wilcox.test(x1,x2)$p.value)
  #return(list(c(median(x1,na.rm = T),median(x2,na.rm = T),median(x3,na.rm=T)),ifelse(p>0.05,"ns",ifelse(p<=0.001,"***",ifelse(p<=0.01,"**","*")))))
  return(list(c(median(x1,na.rm = T),median(x2,na.rm = T),median(x3,na.rm=T)),p))
})
# fc_gg[, 1]: AG
# [1]  0.3352476 -0.2091013 -0.3916919
# fc_gg[, 1]: TE
# [1] 1.0007348 0.3751906 0.2439568
# fc_gg[, 1]: OV
# [1] 0.206486982 0.205991432 0.008357767
######### testis biased genes
expr=read.table("190508_flyatlas_allgenes_alltissues_uniq_FPKM.txt",sep=" ")
expr$tau = apply(expr[,2:28],1,function(x){
  x= log2(x+1)
  return(sum(1-x/max(x))/(length(x)-1))
})

expr$tissue=apply(expr[,2:29],1,function(x){
  if(!is.nan(x[28])&x[28]>=0.8){
    #browser()
    return(names(which.max(x)))
  }else if(!is.nan(x[28])&x[28]<=0.2){return("Widely")}
  return("Other")
})
table(expr$tissue)
asg = subset(expr,tissue=="MaleAccessory")[,1]
tsg = subset(expr,tissue=="MaleTestis")[,1]
osg = subset(expr,tissue=="FemaleOvary")[,1]

write.table(tsg,"testis_specific_genes.txt",row.names = F,col.names = F,quote=F,sep = "\t")
write.table(osg,"ovary_specific_genes.txt",row.names = F,col.names = F,quote=F,sep = "\t")
ts_partners = setdiff(unlist(ppi[ppi[,1]%in%tsg|ppi[,2]%in%tsg,]),tsg)
os_partners = setdiff(unlist(ppi[ppi[,1]%in%osg|ppi[,2]%in%osg,]),osg)

tpm_ts = tpm_m[intersect(rownames(tpm_m),fbgn2sym[tsg,1]),]
tpm_os = tpm_m[intersect(rownames(tpm_m),fbgn2sym[osg,1]),]
tpm_tsp = tpm_m[intersect(rownames(tpm_m),fbgn2sym[ts_partners,1]),]
tpm_osp = tpm_m[intersect(rownames(tpm_m),fbgn2sym[os_partners,1]),]

##### heatmaps of tissue specific genes vs their partners
plot_heatmap(tpm_ts,"TSG.heatmap.pdf",k=4,h=10)
plot_heatmap(tpm_tsp,"TSG.partners.heatmap.pdf",k=4,h=5)
plot_heatmap(tpm_os,"OSG.heatmap.pdf",k=2,h=2.5)
plot_heatmap(tpm_osp,"OSG.partners.heatmap.pdf",h=6)
####
all_gg = rbind(data.frame(melt(tpm_sfp),type="TSG",index="AG"),data.frame(melt(tpm_partners),type="Partners",index="AG"),data.frame(melt(tpm_ts),type="TSG",index="TE"),data.frame(melt(tpm_tsp),type="Partners",index="TE"),data.frame(melt(tpm_os),type="TSG",index="OV"),data.frame(melt(tpm_osp),type="Partners",index="OV"))
all_gg[,6:7]=matrix(unlist(strsplit(as.character(all_gg[,2]),split="\\.",perl=T)),ncol=2,byrow = T)
all_gg$type=factor(all_gg$type,levels=c("TSG","Partners"))
all_gg$V6 = factor(all_gg$V6, levels=c("dme","dsi","dya","dan","dps"))
all_gg$V7 = factor(all_gg$V7, levels=c("AG","TE","OV"))
cols <- brewer.pal(12,"Paired")[c(2,4,6)]
p <- ggplot() + geom_boxplot(data=all_gg,aes(x=V6,y=log10(value),fill=V7),outlier.size = 0.5) + facet_grid(type~index,scales="free_y") + theme_bw(12)+ scale_fill_manual(values=rev(cols)) + xlab("") + ylab("log10(TPM)")+coord_cartesian(ylim=c(-2,6))
ggsave("boxplotTPMSfpsvsPartners.pdf",width = 10,height = 5)
by(all_gg[,c("value","type","V6","V7")],all_gg[,"index"],function(x){
  by(x[,c("value","V6","V7")],x[,"type"],function(y){
    by(y[,c("value","V7")],y[,"V6"],function(z){
      zag = subset(z,V7=="AG")[,1]
      zte = subset(z,V7=="TE")[,1]
      zov = subset(z,V7=="OV")[,1]
      p=c(wilcox.test(zag,zte)$p.value,wilcox.test(zag,zov)$p.value,wilcox.test(zte,zov)$p.value)
      return(ifelse(p<=0.001,"***",ifelse(p<=0.01,"**",ifelse(p<=0.05,"*","ns"))))
    })
})
})


########### Fold changes of tissue specific genes vs their partners

## testis specific genes
others = setdiff(setdiff(rownames(fc_all),tsg),ts_partners)
fc_ts = fc_all[intersect(rownames(fc_all),tsg),]
fc_others = fc_all[others,]
fc_tsp = fc_all[intersect(rownames(fc_all),ts_partners),]

fc_gg = rbind(data.frame(melt(fc_ts),type="TSGs"),data.frame(melt(fc_tsp),type="PPI partners"),data.frame(melt(fc_others),type="Other genes"))
fc_gg$variable = factor(fc_gg$variable,levels=c("AG","TE","OV"))
fc_gg$type = factor(fc_gg$type,levels=c("TSGs","PPI partners","Other genes"))
cols <- c("#808A87",brewer.pal(12,"Paired")[c(1,2)])
# p<-ggplot() + geom_boxplot(data=fc_gg,aes(x=variable,y=value,color=type),outlier.size=0.4) + theme_bw(10) + ylab("log2FC (dme/dps)")+ theme(axis.text.x=element_text(angle=45,hjust=1)) + scale_color_manual(values=rev(cols)) + coord_cartesian(ylim = c(-10,20))
# ggsave("TSG_ppi_evolution2.pdf",width = 3,height = 3)

fc_gg$type = factor(fc_gg$type, levels=c("Other genes","PPI partners","TSGs"))
fc_medians = melt(tapply(fc_gg[,2],fc_gg[,c(1,3)],function(x){median(x,na.rm=T)}))
p<-ggplot()+ geom_vline(data=fc_medians,aes(xintercept=value,color=type),linetype="dashed") +stat_ecdf(data=fc_gg,aes(x=value,color=type)) + facet_grid(variable~.) + theme_bw(10) + ylab("Density")+ theme(axis.text.x=element_text(angle=45,hjust=1), panel.grid=element_blank()) + scale_color_manual(values=cols) + xlim(c(-5,5))+ xlab("log2 FC (dme/dps)")
ggsave("TSG_ppi_evolution_ecdf.pdf",width = 3,height = 3)


by(fc_gg[,2:3],fc_gg[,1],function(x){
  x1 = subset(x,type=="TSGs")[,1]
  x2 = subset(x,type=="PPI partners")[,1]
  x3 = subset(x,type=="Other genes")[,1]
  p = c(wilcox.test(x1,x3)$p.value,wilcox.test(x2,x3)$p.value,wilcox.test(x1,x2)$p.value)
  #return(list(c(median(x1,na.rm = T),median(x2,na.rm = T),median(x3,na.rm=T)),ifelse(p>0.05,"ns",ifelse(p<=0.001,"***",ifelse(p<=0.01,"**","*")))))
  return(list(c(median(x1,na.rm = T),median(x2,na.rm = T),median(x3,na.rm=T)),p))
})
# fc_gg[, 1]: AG
# [1] -2.7640394 -0.2393230 -0.2693305
# fc_gg[, 1]: TE
# [1] -0.1332851  0.1442526  0.3244643
# fc_gg[, 1]: OV
# [1] -0.48405327  0.06234274  0.03391189
######ovary specific genes
others = setdiff(setdiff(rownames(fc_all),osg),os_partners)
fc_os = fc_all[intersect(rownames(fc_all),osg),]
fc_others = fc_all[others,]
fc_osp = fc_all[intersect(rownames(fc_all),os_partners),]

fc_gg = rbind(data.frame(melt(fc_os),type="OSGs"),data.frame(melt(fc_osp),type="PPI partners"),data.frame(melt(fc_others),type="Other genes"))
fc_gg$variable = factor(fc_gg$variable,levels=c("AG","TE","OV"))
fc_gg$type = factor(fc_gg$type,levels=c("OSGs","PPI partners","Other genes"))
cols <- c("#808A87",brewer.pal(12,"Paired")[c(3,4)])
# p<-ggplot() + geom_boxplot(data=fc_gg,aes(x=variable,y=value,color=type),outlier.size=0.4) + theme_bw(10) + ylab("log2FC (dme/dps)")+ theme(axis.text.x=element_text(angle=45,hjust=1)) + scale_color_manual(values=rev(cols))+ coord_cartesian(ylim = c(-10,20))
# ggsave("OSG_ppi_evolution2.pdf",width = 3.6,height = 3.2)


fc_gg$type = factor(fc_gg$type, levels=c("Other genes","PPI partners","OSGs"))
fc_medians = melt(tapply(fc_gg[,2],fc_gg[,c(1,3)],function(x){median(x,na.rm=T)}))
p<-ggplot()+ geom_vline(data=fc_medians,aes(xintercept=value,color=type),linetype="dashed") +stat_ecdf(data=fc_gg,aes(x=value,color=type)) + facet_grid(variable~.) + theme_bw(10) + ylab("Density")+ theme(axis.text.x=element_text(angle=45,hjust=1), panel.grid=element_blank()) + scale_color_manual(values=cols) + xlim(c(-5,5))+ xlab("log2 FC (dme/dps)")
ggsave("OSG_ppi_evolution_ecdf.pdf",width = 3,height = 3)

by(fc_gg[,2:3],fc_gg[,1],function(x){
  x1 = subset(x,type=="OSGs")[,1]
  x2 = subset(x,type=="PPI partners")[,1]
  x3 = subset(x,type=="Other genes")[,1]
  p = c(wilcox.test(x1,x3)$p.value,wilcox.test(x2,x3)$p.value,wilcox.test(x1,x2)$p.value)
  #return(list(c(median(x1,na.rm = T),median(x2,na.rm = T),median(x3,na.rm=T)),ifelse(p>0.05,"ns",ifelse(p<=0.001,"***",ifelse(p<=0.01,"**","*")))))
  return(list(c(median(x1,na.rm = T),median(x2,na.rm = T),median(x3,na.rm=T)),p))
})
# fc_gg[, 1]: AG
# [1] -1.9377961 -0.2879079 -0.3756119
# fc_gg[, 1]: TE
# [1] -0.8295149  0.3137246  0.2648479
# fc_gg[, 1]: OV
# [1] -0.117593051  0.164844839  0.008579594

# fc_gg = rbind(data.frame(logFC=fc_osp[,1],type="Partners in AG"),data.frame(logFC=fc_others[,1],type="Others in AG"),data.frame(logFC=fc_osp[,2],type="Partners in TE"),data.frame(logFC=fc_others[,2],type="Others in TE"),data.frame(logFC=fc_os[,3],type="OSG in OV"),data.frame(logFC=fc_others[,3],type="Others in OV"))
# fc_gg$type=factor(fc_gg$type,levels=c("Partners in AG","Others in AG","Partners in TE","Others in TE","OSG in OV","Others in OV"))
# cols <- brewer.pal(12,"Paired")[c(1,2,3,4,5,6)]
# p <- ggplot(fc_gg,aes(x=type,y=logFC,color=type)) + geom_boxplot(position=position_dodge(0.2),outlier.size = 0.5) + theme_bw(12) + ylab("log2FC (dme/dps)")+ theme(axis.text.x=element_text(angle=45,hjust=1)) + scale_color_manual(values=rev(cols))+ylim(c(-10,25))
# ggsave("OSG_ppi_evolution.pdf",width = 3.6,height = 3.5)
# 
# ixes = combn(c("Partners in AG","Others in AG","Partners in TE","Others in TE","OSG in OV","Others in OV"),2)
# for(i in 1:ncol(ixes)){
#   x1 = subset(fc_gg,type==ixes[1,i])[,1]
#   x2 = subset(fc_gg,type==ixes[2,i])[,1]
#   cat(ixes[,i],": ")
#   p = wilcox.test(x1,x2)$p.value
#   print(ifelse(p<=0.001,"***",ifelse(p<=0.01,"**",ifelse(p<=0.05,"*","n.s."))))
# }
# Partners in AG Others in AG : [1] "n.s."
# Partners in AG Partners in TE : [1] "n.s."
# Partners in AG Others in TE : [1] "n.s."
# Partners in AG OSG in OV : [1] "n.s."
# Partners in AG Others in OV : [1] "n.s."
# Others in AG Partners in TE : [1] "n.s."
# Others in AG Others in TE : [1] "n.s."
# Others in AG OSG in OV : [1] "n.s."
# Others in AG Others in OV : [1] "***"
# Partners in TE Others in TE : [1] "n.s."
# Partners in TE OSG in OV : [1] "n.s."
# Partners in TE Others in OV : [1] "*"
# Others in TE OSG in OV : [1] "n.s."
# Others in TE Others in OV : [1] "**"
# OSG in OV Others in OV : [1] "n.s."
by(fc_gg[,1],fc_gg[,2],function(x){median(x,na.rm=T)})
# fc_gg[, 2]: Partners in AG
# [1] 0.09492316
# ---------------------------------------------------------------------- 
#   fc_gg[, 2]: non-partners in AG
# [1] 0.007602144
# ---------------------------------------------------------------------- 
#   fc_gg[, 2]: Partners in TE
# [1] 0.07464907
# ---------------------------------------------------------------------- 
#   fc_gg[, 2]: non-partners in TE
# [1] -0.008618342
# ---------------------------------------------------------------------- 
#   fc_gg[, 2]: OSG in OV
# [1] -0.1387228
# ---------------------------------------------------------------------- 
#   fc_gg[, 2]: non-OSG in OV
# [1] -0.02304853
##################
##################
# correlation analysis
sfp = read.table("wigby2020.txt")
tsg = read.table("testis_specific_genes.txt",stringsAsFactors = F)
osg = read.table("ovary_specific_genes.txt",stringsAsFactors = F)
others = setdiff(setdiff(setdiff(rownames(tpm_m),sfp),fbgn2sym[tsg[,1],1]),fbgn2sym[osg[,1],1])

ppi = read.table("physical_interactions.tab")
ppi2 = cbind(fbgn2sym[ppi[,1],1],fbgn2sym[ppi[,2],1])
ppi2 = ppi2[-unique(which(is.na(ppi2),arr.ind=T)[,1]),]

cal_corr <- function(genelist,ppi2,tpm_m,ix1=c("dme.AG","dsi.AG","dya.AG","dan.AG","dps.AG"),ix2=c("dme.OV","dsi.OV","dya.OV","dan.OV","dps.OV")){
  r=vector()
  genelist = intersect(genelist,rownames(tpm_m))
  for (i in genelist){
    #browser()
    partners = intersect(setdiff(unlist(ppi2[ppi2[,1]==i|ppi2[,2]==i,]),i),rownames(tpm_m))
    if(length(partners)!=0){
      tpm_sfp = tpm_m[i,ix1]
      tpm_partners = tpm_m[partners,ix2]
      if(!is.null(nrow(tpm_partners))){
        r=c(r,apply(tpm_partners,1,function(a){cor(tpm_sfp,a,  method = "spearman")}))}
      else{r=c(r,cor(tpm_sfp,tpm_partners,method = "spearman"))}
    }
  }
  return(r)
}

x1 = rbind(data.frame(r=cal_corr(sfp[,1],ppi2,tpm_m,ix1=c("dme.AG","dsi.AG","dya.AG","dan.AG","dps.AG"),ix2=c("dme.AG","dsi.AG","dya.AG","dan.AG","dps.AG")),type="Sfps",ix="AG"),data.frame(r=cal_corr(others,ppi2,tpm_m,ix1=c("dme.AG","dsi.AG","dya.AG","dan.AG","dps.AG"),ix2=c("dme.AG","dsi.AG","dya.AG","dan.AG","dps.AG")),type="Others",ix="AG"))

x2 = rbind(data.frame(r=cal_corr(sfp[,1],ppi2,tpm_m,ix1=c("dme.AG","dsi.AG","dya.AG","dan.AG","dps.AG"),ix2=c("dme.TE","dsi.TE","dya.TE","dan.TE","dps.TE")),type="Sfps",ix="TE"),data.frame(r=cal_corr(others,ppi2,tpm_m,ix1=c("dme.AG","dsi.AG","dya.AG","dan.AG","dps.AG"),ix2=c("dme.TE","dsi.TE","dya.TE","dan.TE","dps.TE")),type="Others",ix="TE"))

x3 = rbind(data.frame(r=cal_corr(sfp[,1],ppi2,tpm_m,ix1=c("dme.AG","dsi.AG","dya.AG","dan.AG","dps.AG"),ix2=c("dme.OV","dsi.OV","dya.OV","dan.OV","dps.OV")),type="Sfps",ix="OV"),data.frame(r=cal_corr(others,ppi2,tpm_m,ix1=c("dme.AG","dsi.AG","dya.AG","dan.AG","dps.AG"),ix2=c("dme.OV","dsi.OV","dya.OV","dan.OV","dps.OV")),type="Others",ix="OV"))

x = rbind(x1,x2,x3)
x$type = factor(x$type, levels=c("Sfps","Others"))
x$ix = factor(x$ix,levels=c("AG","TE","OV"))
cols <- c(brewer.pal(9,"Set1")[9],brewer.pal(12,"Paired")[6])
p1 <- ggplot() + geom_boxplot(data=x,aes(x=ix,y=r,color=type)) + xlab("Organ that PPI partners expressed") + ylab("Spearman's rho\nSfps ~ PPI partners") + theme_bw(10) + ylim(-1.0,1.25) + scale_color_manual(values=rev(cols))+theme(panel.grid=element_blank())
ggsave("Sfps_PPI_partners_correlation.pdf",width = 2.8, height = 2.4)
by(x[,1:2],x[,3],function(x){
  wilcox.test(x[,1]~x[,2])
})
##################
y1 = rbind(data.frame(r=cal_corr(fbgn2sym[tsg[,1],1],ppi2,tpm_m,ix1=c("dme.TE","dsi.TE","dya.TE","dan.TE","dps.TE"),ix2=c("dme.AG","dsi.AG","dya.AG","dan.AG","dps.AG")),type="TSGs",ix="AG"),data.frame(r=cal_corr(others,ppi2,tpm_m,ix1=c("dme.TE","dsi.TE","dya.TE","dan.TE","dps.TE"),ix2=c("dme.AG","dsi.AG","dya.AG","dan.AG","dps.AG")),type="Others",ix="AG"))

y2 = rbind(data.frame(r=cal_corr(fbgn2sym[tsg[,1],1],ppi2,tpm_m,ix1=c("dme.TE","dsi.TE","dya.TE","dan.TE","dps.TE"),ix2=c("dme.TE","dsi.TE","dya.TE","dan.TE","dps.TE")),type="TSGs",ix="TE"),data.frame(r=cal_corr(others,ppi2,tpm_m,ix1=c("dme.TE","dsi.TE","dya.TE","dan.TE","dps.TE"),ix2=c("dme.TE","dsi.TE","dya.TE","dan.TE","dps.TE")),type="Others",ix="TE"))

y3 = rbind(data.frame(r=cal_corr(fbgn2sym[tsg[,1],1],ppi2,tpm_m,ix1=c("dme.TE","dsi.TE","dya.TE","dan.TE","dps.TE"),ix2=c("dme.OV","dsi.OV","dya.OV","dan.OV","dps.OV")),type="TSGs",ix="OV"),data.frame(r=cal_corr(others,ppi2,tpm_m,ix1=c("dme.TE","dsi.TE","dya.TE","dan.TE","dps.TE"),ix2=c("dme.OV","dsi.OV","dya.OV","dan.OV","dps.OV")),type="Others",ix="OV"))

y = rbind(y1,y2,y3)
y$type = factor(y$type,levels=c("TSGs","Others"))
y$ix = factor(y$ix,levels=c("AG","TE","OV"))
cols <- c(brewer.pal(9,"Set1")[9],brewer.pal(12,"Paired")[2])
p2 <- ggplot() + geom_boxplot(data=y,aes(x=ix,y=r,color=type)) + xlab("Organ that PPI partners expressed") + ylab("Spearman's rho\nTSGs ~ PPI partners") + theme_bw(10) + ylim(-1.0,1.25) + scale_color_manual(values=rev(cols))+theme(panel.grid=element_blank())
ggsave("TSG_PPI_partners_correlation.pdf",width = 2.8, height = 2.4)
by(y[,1:2],y[,3],function(x){
  wilcox.test(x[,1]~x[,2])
})
##############
z1 = rbind(data.frame(r=cal_corr(fbgn2sym[osg[,1],1],ppi2,tpm_m,ix1=c("dme.OV","dsi.OV","dya.OV","dan.OV","dps.OV"),ix2=c("dme.AG","dsi.AG","dya.AG","dan.AG","dps.AG")),type="OSGs",ix="AG"),data.frame(r=cal_corr(others,ppi2,tpm_m,ix1=c("dme.OV","dsi.OV","dya.OV","dan.OV","dps.OV"),ix2=c("dme.AG","dsi.AG","dya.AG","dan.AG","dps.AG")),type="Others",ix="AG"))

z2 = rbind(data.frame(r=cal_corr(fbgn2sym[osg[,1],1],ppi2,tpm_m,ix1=c("dme.OV","dsi.OV","dya.OV","dan.OV","dps.OV"),ix2=c("dme.TE","dsi.TE","dya.TE","dan.TE","dps.TE")),type="OSGs",ix="TE"),data.frame(r=cal_corr(others,ppi2,tpm_m,ix1=c("dme.OV","dsi.OV","dya.OV","dan.OV","dps.OV"),ix2=c("dme.TE","dsi.TE","dya.TE","dan.TE","dps.TE")),type="Others",ix="TE"))

z3 = rbind(data.frame(r=cal_corr(fbgn2sym[osg[,1],1],ppi2,tpm_m,ix1=c("dme.OV","dsi.OV","dya.OV","dan.OV","dps.OV"),ix2=c("dme.OV","dsi.OV","dya.OV","dan.OV","dps.OV")),type="OSGs",ix="OV"),data.frame(r=cal_corr(others,ppi2,tpm_m,ix1=c("dme.OV","dsi.OV","dya.OV","dan.OV","dps.OV"),ix2=c("dme.OV","dsi.OV","dya.OV","dan.OV","dps.OV")),type="Others",ix="OV"))

z = rbind(z1,z2,z3)
z$type = factor(z$type,levels=c("OSGs","Others"))
z$ix = factor(z$ix,levels = c("AG","TE","OV"))
cols <- c(brewer.pal(9,"Set1")[9],brewer.pal(12,"Paired")[4])
p3 <- ggplot() + geom_boxplot(data=z,aes(x=ix,y=r,color=type)) + xlab("Organ that PPI partners expressed") + ylab("Spearman's rho\nOSGs ~ PPI partners") + theme_bw(10) + ylim(-1.0,1.25) + scale_color_manual(values=rev(cols))+theme(panel.grid=element_blank())
ggsave("OSG_PPI_partners_correlation.pdf",width = 2.8, height = 2.4)

by(z[,1:2],z[,3],function(x){
  wilcox.test(x[,1]~x[,2])
})

#######
library(viridis)
library(ggridges)
ggplot(data=x, aes(x = r, y = ix, fill = type,colour = type)) +geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01, gradient_lwd = 1.) + scale_color_manual(values = cols)+
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.05, 0)) + 
  scale_fill_viridis(name = "type", option = "C",discrete = T,alpha = 0.06) +
  labs(title = 'Correlation between Sfps and PPI partners',x="Spearman's correlation") +
  theme_ridges(font_size = 10, grid = TRUE) + theme(axis.title.y = element_blank())
ggsave("Sfps_rho_3d.pdf",width = 5,height = 3)

