library(ggplot2)
#library(eulerr)
library(reshape)
#library(UpSetR)
library(ComplexHeatmap)
#library(ClusterGVis)
options(stringsAsFactors = F)
#################
#setwd("C:/Users/Administrator/OneDrive/10.shangrui/ConsDivEvolution")
setwd("/Users/aimeidai/OneDrive/10.shangrui/ConsDivEvolution")
dp_ta =  read.table("dps_AGvsTE.ConsDiv.txt")
dp_to =  read.table("dps_TEvsOV.ConsDiv.txt")
dp_ao =  read.table("dps_AGvsOV.ConsDiv.txt")

da_ta =  read.table("dan_AGvsTE.ConsDiv.txt")
da_to =  read.table("dan_TEvsOV.ConsDiv.txt")
da_ao =  read.table("dan_AGvsOV.ConsDiv.txt")

dy_ta =  read.table("dya_AGvsTE.ConsDiv.txt")
dy_to =  read.table("dya_TEvsOV.ConsDiv.txt")
dy_ao =  read.table("dya_AGvsOV.ConsDiv.txt")

ds_ta =  read.table("dsi_AGvsTE.ConsDiv.txt")
ds_to =  read.table("dsi_TEvsOV.ConsDiv.txt")
ds_ao =  read.table("dsi_AGvsOV.ConsDiv.txt")

genes = unique(c(rownames(dp_ta),rownames(da_ta),rownames(dy_ta),rownames(ds_ta)))
ta = cbind(dp=dp_ta[genes,"type"],da=da_ta[genes,"type"],dy=dy_ta[genes,"type"],ds=ds_ta[genes,"type"])
rownames(ta)=genes
ta[which(is.na(ta),arr.ind = T)]=0

ta_ix=apply(ta,1,function(x){
  if(sum(x==0)==4){return("ns")}
  else if(abs(sum(x==1)-sum(x==-1))<=1&(sum(x==1)!=0&sum(x==-1)!=0)){return("ambiguous")}
  else if((sum(x==0)==3|sum(x==0)==2)&sum(x)!=0){return("weak")}
  else if ((abs(sum(x))==4|abs(sum(x))==3|abs(sum(x))==2&sum(x!=0)==4)){return("strong")}
  return("opps")
})

genes = unique(c(rownames(dp_to),rownames(da_to),rownames(dy_to),rownames(ds_to)))
to = cbind(dp=dp_to[genes,"type"],da=da_to[genes,"type"],dy=dy_to[genes,"type"],ds=ds_to[genes,"type"])
rownames(to)=genes
to[which(is.na(to),arr.ind = T)]=0
to_ix=apply(to,1,function(x){
  if(sum(x==0)==4){return("ns")}
  else if(abs(sum(x==1)-sum(x==-1))<=1&(sum(x==1)!=0&sum(x==-1)!=0)){return("ambiguous")}
  else if((sum(x==0)==3|sum(x==0)==2)&sum(x)!=0){return("weak")}
  else if ((abs(sum(x))==4|abs(sum(x))==3|abs(sum(x))==2&sum(x!=0)==4)){return("strong")}
  return("opps")
})

genes = unique(c(rownames(dp_ao),rownames(da_ao),rownames(dy_ao),rownames(ds_ao)))
ao = cbind(dp=dp_ao[genes,"type"],da=da_ao[genes,"type"],dy=dy_ao[genes,"type"],ds=ds_ao[genes,"type"])
rownames(ao)=genes
ao[which(is.na(ao),arr.ind = T)]=0

ao_ix=apply(ao,1,function(x){
  if(sum(x==0)==4){return("ns")}
  else if(abs(sum(x==1)-sum(x==-1))<=1&(sum(x==1)!=0&sum(x==-1)!=0)){return("ambiguous")}
  else if((sum(x==0)==3|sum(x==0)==2)&sum(x)!=0){return("weak")}
  else if ((abs(sum(x))==4|abs(sum(x))==3|abs(sum(x))==2&sum(x!=0)==4)){return("strong")}
  return("opps")
})


require(RColorBrewer)
library(circlize)
#col_fun = colorRamp2(c(-3,0,3), c("blue", "white", "red"))
genes=unique(c(rownames(ta),rownames(ao),rownames(to)))
all = data.frame(cbind(AGvsTE=rowSums(ta)[genes],AGvsOV=rowSums(ao)[genes],TEvsOV=rowSums(to)[genes]))
rownames(all)=genes
all[which(is.na(all),arr.ind = T)]=0
# all$typeAGvsTE=ta_ix[rownames(all)]
# all$typeAGvsOV=ao_ix[rownames(all)]
# all$typeTEvsOV=to_ix[rownames(all)]
# all[which(is.na(all),arr.ind = T)]="ns"
# write.table(all,"all_coevolve_scores.txt",sep="\t",row.names = T,col.names = T)

pv = apply(ifelse(all[,1:3]<2&all[,1:3]>-2,0,ifelse(all[,1:3]>=2,1,-1)),1,function(x){paste(x,collapse="")})
write.table(pv,"coevolve_state_by_species.txt",sep="\t",row.names = T,col.names = F,quote=F)

########
# colors = structure(c(brewer.pal(12, "Paired")[c(5,6)],brewer.pal(9,"YlOrRd")[6],brewer.pal(11, "BrBG")[5],brewer.pal(11,"RdBu")[6:10]), names = c(2,4,3,1,0,-1,-2,-3,-4))
# colors=colors[c(2,3,1,4:9)]
# s=Heatmap(all[,1:3],name = "State", row_km =6, show_row_names = FALSE,show_row_dend =F, show_column_names = TRUE,cluster_columns=F,show_column_dend =F,col = colors)
# pdf("transition_heatmap.pdf",height = 5,width=1.8)
# set.seed(123)
# s=draw(s)
# dev.off()
# 
# rcl.s <- melt(row_order(s))
# rcl.s$geneName = rownames(all)[rcl.s[,1]]
# write.table(rcl.s[,2:3],"coevolve_state_by_species.txt",sep="\t",row.names = F,col.names = F,quote=F)



##### strongly evovled
all_s = list()
states=names(sort(table(pv),decreasing = T))[1:11]
for(i in states){
  all_s[[i]] = all[names(pv[pv==i]),1:3]
}
all_s[["Others"]] = all[names(pv[!pv%in%states]),1:3]

colors = structure(c(brewer.pal(12, "Paired")[c(5,6)],brewer.pal(9,"YlOrRd")[6],brewer.pal(11, "BrBG")[5],brewer.pal(11,"RdBu")[6:10]), names = c(2,4,3,1,0,-1,-2,-3,-4))
colors=colors[c(2,3,1,4:9)]

all_heatmap= list()
for(i in setdiff(names(all_s),"n.s.")){
  all_heatmap[[i]]=Heatmap(all_s[[i]],name = "State", row_km =1, show_row_names = FALSE,show_row_dend =F, show_column_names = TRUE,cluster_columns=F,show_column_dend =F,col = colors)# + rowAnnotation(link = anno_mark(at = which(rownames(all_s[[i]])%in%acps[,3]), labels = acps[acps[,3]%in%rownames(all_s[[i]]),1], labels_gp = gpar(fontsize = 4), padding = unit(1, "mm")))
  # pdf(paste("transition_heatmap",i,".pdf",sep=""),height=2,width=3)
  # draw(all_heatmap[[i]])
  # dev.off()
}
s=all_heatmap[[1]]%v%all_heatmap[[2]]

for(i in 3:length(all_heatmap)){
  s = s%v%all_heatmap[[i]]
}

pdf("transition_heatmap.pdf",height = 6.5,width=1.8)
set.seed(123)
s=draw(s)
dev.off()

rcl.s = melt(pv)
sum(rcl.s$value %in% c("010","100","001","111","110","101","011"))/nrow(rcl.s)
sum(rcl.s$value %in% c("010","100","001"))/nrow(rcl.s)
sum(rcl.s$value %in% c("111","110","101","011"))/nrow(rcl.s)
sum(rcl.s$value %in% c("-1-1-1","-1-10","-10-1","0-1-1"))/nrow(rcl.s)
sum(rcl.s$value %in% c("0-10","-100","00-1"))/nrow(rcl.s)
sum(rcl.s$value %in% c("000"))/nrow(rcl.s)
##########################
require(RColorBrewer)
cols <- c(rev(brewer.pal(9, "Set1"))[c(1,8,9)],brewer.pal(8, "Dark2")[1:8])
rcl.s = read.table("coevolve_state_by_species.txt",colClasses = "character")
states = names(sort(table(rcl.s$V2),decreasing=T))[1:11]
rcl.s$type = ifelse(rcl.s$V2%in%states,rcl.s$V2,"Others")
acps = read.table("wigby2020.txt")
rownames(acps)=acps[,3]
rcl.s$ifACP=ifelse(rcl.s$V1 %in% acps[,3],"Yes","No")
x=table(rcl.s[,c(3,4)])
r= subset(data.frame(x,ratio=x[,2]/(x[,1]+x[,2])),ifACP=="Yes")
r$total = x[r[,1],2]+x[r[,1],1]
r$type = factor(r$type,levels=c(states,"Others"))
p <- ggplot(subset(r,ifACP=="Yes")) + geom_bar(aes(x=type,y=ratio),stat="identity",fill="white",color="black")+geom_text(aes(x=type,y=ratio+0.002,label=paste(Freq,"/",total,sep="")),size=1.6)+ theme_bw(10) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="right", legend.title=element_blank()) + labs(title="",y="Proportion of Acps", x="Coevolution state") 
ggsave("ACPs_distribution.pdf",width=3,height=2.5)
#######################
x0 = x["000",]
for(i in c(states,"Others")){
  if (i == "000"){next}
  x1 = x[i,]
  m = cbind(x0,x1)
  cat(i,":",fisher.test(m)$p.value,"\t")
}
#010 : 0.2469496 	100 : 0.4774118 	001 : 0.6218072 	111 : 0.1514787 	110 : 0.1691973 	101 : 0.1573232 	-100 : 0.8071061 	00-1 : 0.7977928 	011 : 0.04825052 	0-10 : 0.0740356 	Others : 0.4364274 
###########################
### testis biased genes
tsg = read.table("../Sfp_SBG/testis_specific_genes.txt")
osg = read.table("../Sfp_SBG/ovary_specific_genes.txt")

rcl.s$ifTSG = ifelse(rcl.s$V1 %in% tsg[,1],"Yes","No")
rcl.s$ifOSG = ifelse(rcl.s$V1 %in% osg[,1],"Yes","No")

##### TSG
x=table(rcl.s[,c(3,5)])
r= subset(data.frame(x,ratio=x[,2]/(x[,1]+x[,2])),ifTSG=="Yes")
r$total = x[r[,1],2]+x[r[,1],1]
r$type = factor(r$type,levels=c(states,"Others"))
p <- ggplot(subset(r,ifTSG=="Yes")) + geom_bar(aes(x=type,y=ratio),stat="identity",fill="white",color="black")+geom_text(aes(x=type,y=ratio+0.008,label=paste(Freq,"/",total,sep="")),size=1.6)+ theme_bw(10) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="right", legend.title=element_blank()) + labs(title="",y="Proportion of TSGs", x="Coevolution state") 
ggsave("TSGs_distribution.pdf",width=3,height=2.5)

x0 = x["000",]
for(i in c(states,"Others")){
  if (i == "000"){next}
  x1 = x[i,]
  m = cbind(x0,x1)
  cat(i,": ",fisher.test(m)$p.value,"\n")
}
# 010 :  3.292323e-09 
# 100 :  0.09613695 
# 001 :  0.3631113 
# 111 :  0.002266209 
# 110 :  1.904852e-08 
# 101 :  3.053812e-07 
# -100 :  2.365187e-08 
# 00-1 :  0.002495242 
# 011 :  0.259453 
# 0-10 :  0.0002350902 
# Others :  3.779766e-07 
###### OSG
x=table(rcl.s[,c(3,6)])
r= subset(data.frame(x,ratio=x[,2]/(x[,1]+x[,2])),ifOSG=="Yes")
r$total = x[r[,1],2]+x[r[,1],1]
r$type = factor(r$type,levels=c(states,"Others"))
p <- ggplot(subset(r,ifOSG=="Yes")) + geom_bar(aes(x=type,y=ratio),stat="identity",fill="white",color="black")+geom_text(aes(x=type,y=ratio+0.0005,label=paste(Freq,"/",total,sep="")),size=1.6)+ theme_bw(10) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5), plot.title = element_text(hjust = 0.5), panel.grid=element_blank(), legend.position="right", legend.title=element_blank()) + labs(title="",y="Proportion of OSGs", x="Coevolution state") 
ggsave("OSGs_distribution.pdf",width=3,height=2.5)

x0 = x["000",]
for(i in c(states,"Others")){
  if (i == "000"){next}
  x1 = x[i,]
  m = cbind(x0,x1)
  cat(i,": ",fisher.test(m)$p.value,"\n")
}
# 010 :  0.001282268 
# 100 :  0.6880281 
# 001 :  0.2655702 
# 111 :  0.2415449 
# 110 :  0.5196444 
# 101 :  1 
# -100 :  0.08885544 
# 00-1 :  0.6256741 
# 011 :  0.4353238 
# 0-10 :  0.2610144 
# Others :  0.134637 

############################
#setwd("C:/Users/Administrator/OneDrive/10.shangrui/coevolution")
#setwd("/Users/aimeidai/OneDrive/10.shangrui/coevolution")
ag = read.table("../coevolution/AG_cluster.txt",row.names = 1)
te = read.table("../coevolution/TE_cluster.txt",row.names = 1)
ov = read.table("../coevolution/OV_cluster.txt",row.names = 1)
rcl.s$TE = te[rcl.s[,1],1]
rcl.s$AG = ag[rcl.s[,1],1]
rcl.s$OV = ov[rcl.s[,1],1]

x_te = table(rcl.s[,c(3,4)])
r_te = melt(x_te/rowSums(x_te))

x_ag = table(rcl.s[,c(3,5)])
r_ag = melt(x_ag/rowSums(x_ag))

x_ov = table(rcl.s[,c(3,6)])
r_ov = melt(x_ov/rowSums(x_ov))

all = merge(r_te,r_ag,by.x=c(1,2),by.y=c(1,2))
all = merge(all,r_ov,by.x=c(1,2),by.y=c(1,2))
colnames(all)=c("cluster","type","TE","AG","OV")
all_gg = melt(all,id=c("cluster","type"))
all_gg$cluster=factor(all_gg$cluster,levels=c(states,"Others"))
all_gg$type = factor(all_gg$type,levels=1:6)

require(RColorBrewer)
colors = c(brewer.pal(12, "Paired")[c(5,6)],brewer.pal(11,"RdBu")[c(2,6,8,10)])
names(colors) = c(4,5,6,3,2,1)
#cols <- c(rev(brewer.pal(9, "Set1"))[c(1,8,9)],brewer.pal(8, "Dark2")[1:8])
p <- ggplot(all_gg) + geom_bar(aes(x=cluster,y=value,fill=type),stat="identity") + facet_wrap(~variable) + scale_fill_manual(values=colors)+ theme(axis.text.x=element_text(angle=90,hjust=1,vjust = 0.5))+ylab("Fraction of genes")
ggsave("correlationBetweenTwoCoevolution.pdf",width=6,height = 4)
################################################








ta_cons = data.frame(row.names=unique(c(rownames(subset(dp_ta,type==1)),rownames(subset(da_ta,type==1)),rownames(subset(dy_ta,type==1)),rownames(subset(ds_ta,type==1)))))
ta_cons$gene_name=rownames(ta_cons)
ta_cons$dps = ifelse(ta_cons$gene_name%in%rownames(subset(dp_ta,type==1)),1,0)
ta_cons$dan = ifelse(ta_cons$gene_name%in%rownames(subset(da_ta,type==1)),1,0)
ta_cons$dya = ifelse(ta_cons$gene_name%in%rownames(subset(dy_ta,type==1)),1,0)
ta_cons$dsi = ifelse(ta_cons$gene_name%in%rownames(subset(ds_ta,type==1)),1,0)

ta_div = data.frame(row.names=unique(c(rownames(subset(dp_ta,type==2)),rownames(subset(da_ta,type==2)),rownames(subset(dy_ta,type==2)),rownames(subset(ds_ta,type==2)))))
ta_div$gene_name=rownames(ta_div)
ta_div$dps = ifelse(ta_div$gene_name%in%rownames(subset(dp_ta,type==2)),1,0)
ta_div$dan = ifelse(ta_div$gene_name%in%rownames(subset(da_ta,type==2)),1,0)
ta_div$dya = ifelse(ta_div$gene_name%in%rownames(subset(dy_ta,type==2)),1,0)
ta_div$dsi = ifelse(ta_div$gene_name%in%rownames(subset(ds_ta,type==2)),1,0)

pdf("TEvsAG_upsetR.pdf",width = 3.5, height = 2)
upset(ta_cons, nsets = 4, sets=c("dps","dan","dya","dsi"), keep.order = T, mb.ratio = c(0.5, 0.5), order.by = "degree", decreasing = FALSE)
upset(ta_div, nsets = 4, sets=c("dps","dan","dya","dsi"), keep.order = T, mb.ratio = c(0.5, 0.5), order.by = "degree", decreasing = FALSE)
dev.off()
##############################
to_cons = data.frame(row.names=unique(c(rownames(subset(dp_to,type==1)),rownames(subset(da_to,type==1)),rownames(subset(dy_to,type==1)),rownames(subset(ds_to,type==1)))))
to_cons$gene_name=rownames(to_cons)
to_cons$dps = ifelse(to_cons$gene_name%in%rownames(subset(dp_to,type==1)),1,0)
to_cons$dan = ifelse(to_cons$gene_name%in%rownames(subset(da_to,type==1)),1,0)
to_cons$dya = ifelse(to_cons$gene_name%in%rownames(subset(dy_to,type==1)),1,0)
to_cons$dsi = ifelse(to_cons$gene_name%in%rownames(subset(ds_to,type==1)),1,0)

to_div = data.frame(row.names=unique(c(rownames(subset(dp_to,type==2)),rownames(subset(da_to,type==2)),rownames(subset(dy_to,type==2)),rownames(subset(ds_to,type==2)))))
to_div$gene_name=rownames(to_div)
to_div$dps = ifelse(to_div$gene_name%in%rownames(subset(dp_to,type==2)),1,0)
to_div$dan = ifelse(to_div$gene_name%in%rownames(subset(da_to,type==2)),1,0)
to_div$dya = ifelse(to_div$gene_name%in%rownames(subset(dy_to,type==2)),1,0)
to_div$dsi = ifelse(to_div$gene_name%in%rownames(subset(ds_to,type==2)),1,0)


pdf("TEvsOV_upsetR.pdf",width = 3.5, height = 2)
upset(to_cons, nsets = 4, sets=c("dps","dan","dya","dsi"), keep.order = T, mb.ratio = c(0.5, 0.5), order.by = "degree", decreasing = FALSE)
upset(to_div, nsets = 4, sets=c("dps","dan","dya","dsi"), keep.order = T, mb.ratio = c(0.5, 0.5), order.by = "degree", decreasing = FALSE)
dev.off()
#############################
ao_cons = data.frame(row.names=unique(c(rownames(subset(dp_ao,type==1)),rownames(subset(da_ao,type==1)),rownames(subset(dy_ao,type==1)),rownames(subset(ds_ao,type==1)))))
ao_cons$gene_name=rownames(ao_cons)
ao_cons$dps = ifelse(ao_cons$gene_name%in%rownames(subset(dp_ao,type==1)),1,0)
ao_cons$dan = ifelse(ao_cons$gene_name%in%rownames(subset(da_ao,type==1)),1,0)
ao_cons$dya = ifelse(ao_cons$gene_name%in%rownames(subset(dy_ao,type==1)),1,0)
ao_cons$dsi = ifelse(ao_cons$gene_name%in%rownames(subset(ds_ao,type==1)),1,0)

ao_div = data.frame(row.names=unique(c(rownames(subset(dp_ao,type==2)),rownames(subset(da_ao,type==2)),rownames(subset(dy_ao,type==2)),rownames(subset(ds_ao,type==2)))))
ao_div$gene_name=rownames(ao_div)
ao_div$dps = ifelse(ao_div$gene_name%in%rownames(subset(dp_ao,type==2)),1,0)
ao_div$dan = ifelse(ao_div$gene_name%in%rownames(subset(da_ao,type==2)),1,0)
ao_div$dya = ifelse(ao_div$gene_name%in%rownames(subset(dy_ao,type==2)),1,0)
ao_div$dsi = ifelse(ao_div$gene_name%in%rownames(subset(ds_ao,type==2)),1,0)


pdf("AGvsOV_upsetR.pdf",width = 3.5, height = 2)
upset(ao_cons, nsets = 4, sets=c("dps","dan","dya","dsi"), keep.order = T, mb.ratio = c(0.5, 0.5), order.by = "degree", decreasing = FALSE)
upset(ao_div, nsets = 4, sets=c("dps","dan","dya","dsi"), keep.order = T, mb.ratio = c(0.5, 0.5), order.by = "degree", decreasing = FALSE)
dev.off()

##################
#load("../ATAC-seq/cons.div.RData")
library(ComplexHeatmap)
require(RColorBrewer)
set.seed(123)
colors = structure(c('#56B4E9', "lightgrey", '#E69F00'), names = c("-1", "0", "1"))
cols2 <- brewer.pal(9, "Set1")[c(1,2,4)]
# dps
ix=2
genes = unique(c(ta_cons[rowSums(ta_cons[,2:5])==1&ta_cons[,ix]==1,1],ta_div[rowSums(ta_div[,2:5])==1&ta_div[,ix]==1,1],to_cons[rowSums(to_cons[,2:5])==1&to_cons[,ix]==1,1],to_div[rowSums(to_div[,2:5])==1&to_div[,ix]==1,1],ao_cons[rowSums(ao_cons[,2:5])==1&ao_cons[,ix]==1,1],ao_div[rowSums(ao_div[,2:5])==1&ao_div[,ix]==1,1]))
type_matrix=matrix(nrow=length(genes),ncol=3,dimnames = list(genes,c("TEvsAG","TEvsOV","AGvsOV")))
type_matrix[,"TEvsAG"]=dp_ta[genes,3]
type_matrix[,"TEvsOV"]=dp_to[genes,3]
type_matrix[,"AGvsOV"]=dp_ao[genes,3]
type_matrix[which(is.na(type_matrix),arr.ind=T)]=0
type_matrix[which(type_matrix==2,arr.ind=T)]=-1
# p=pheatmap(type_matrix,show_rownames = F,color = colorRampPalette(c('#56B4E9', "white", '#E69F00'))(10),cutree_rows = 5)
ix = apply(type_matrix,1,function(x){
  ix = ifelse(!any(x==-1),1,ifelse(!any(x==1),-1,0))
  if(ix==1){
    if(sum(x)==1){return(10)}else{return(11)}
  }
  if(ix==-1){
    if(sum(x)==-1){return(-30)}else{return(-31)}
  }
  if(ix==0){
    #browser()
    if(x[1]==1){return(21)}
    if(x[2]==1){return(22)}
    if(x[3]==1){return(23)}
  }
  return(ix)
})
cNo=data.frame(dps=melt(table(ix))[,2])
rownames(cNo)=melt(table(ix))[,1]
m1 = Heatmap(type_matrix[ix==11,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m2 = Heatmap(type_matrix[ix==10,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=3,show_row_dend = F,gap=unit(0,"npc"))
m3 = Heatmap(type_matrix[ix==21,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m4 = Heatmap(type_matrix[ix==22,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m5 = Heatmap(type_matrix[ix==23,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m6 = Heatmap(type_matrix[ix==-30,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m7 = Heatmap(type_matrix[ix==-31,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m = m1%v%m2%v%m3%v%m4%v%m5%v%m6%v%m7

pdf("dp.pheatmap.type.transit.pdf",width = 2.5,height = 6)
draw(m)
dev.off()

save(m1,m2,m3,m4,m5,m6,m7,ix,type_matrix,file="dps_transit.RData")

## dan
ix=3
genes = unique(c(ta_cons[rowSums(ta_cons[,2:5])==1&ta_cons[,ix]==1,1],ta_div[rowSums(ta_div[,2:5])==1&ta_div[,ix]==1,1],to_cons[rowSums(to_cons[,2:5])==1&to_cons[,ix]==1,1],to_div[rowSums(to_div[,2:5])==1&to_div[,ix]==1,1],ao_cons[rowSums(ao_cons[,2:5])==1&ao_cons[,ix]==1,1],ao_div[rowSums(ao_div[,2:5])==1&ao_div[,ix]==1,1]))
type_matrix=matrix(nrow=length(genes),ncol=3,dimnames = list(genes,c("TEvsAG","TEvsOV","AGvsOV")))
type_matrix[,"TEvsAG"]=da_ta[genes,3]
type_matrix[,"TEvsOV"]=da_to[genes,3]
type_matrix[,"AGvsOV"]=da_ao[genes,3]
type_matrix[which(is.na(type_matrix),arr.ind=T)]=0
type_matrix[which(type_matrix==2,arr.ind=T)]=-1
# p=pheatmap(type_matrix,show_rownames = F,color = colorRampPalette(c('#56B4E9', "white", '#E69F00'))(10),cutree_rows = 5)
ix = apply(type_matrix,1,function(x){
  ix = ifelse(!any(x==-1),1,ifelse(!any(x==1),-1,0))
  if(ix==1){
    if(sum(x)==1){return(10)}else{return(11)}
  }
  if(ix==-1){
    if(sum(x)==-1){return(-30)}else{return(-31)}
  }
  if(ix==0){
    #browser()
    if(x[1]==1){return(21)}
    if(x[2]==1){return(22)}
    if(x[3]==1){return(23)}
  }
  return(ix)
})
table(ix)
cNo$dan=melt(table(ix))[,2]
m1 = Heatmap(type_matrix[ix==11,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m2 = Heatmap(type_matrix[ix==10,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=3,show_row_dend = F,gap=unit(0,"npc"))
m3 = Heatmap(type_matrix[ix==21,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m4 = Heatmap(type_matrix[ix==22,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m5 = Heatmap(type_matrix[ix==23,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m6 = Heatmap(type_matrix[ix==-30,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m7 = Heatmap(type_matrix[ix==-31,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m = m1%v%m2%v%m3%v%m4%v%m5%v%m6%v%m7

pdf("da.pheatmap.type.transit.pdf",width = 2,height = 6)
draw(m)
dev.off()
save(m1,m2,m3,m4,m5,m6,m7,ix,type_matrix,file="dan_transit.RData")

## dya
ix=4
genes = unique(c(ta_cons[rowSums(ta_cons[,2:5])==1&ta_cons[,ix]==1,1],ta_div[rowSums(ta_div[,2:5])==1&ta_div[,ix]==1,1],to_cons[rowSums(to_cons[,2:5])==1&to_cons[,ix]==1,1],to_div[rowSums(to_div[,2:5])==1&to_div[,ix]==1,1],ao_cons[rowSums(ao_cons[,2:5])==1&ao_cons[,ix]==1,1],ao_div[rowSums(ao_div[,2:5])==1&ao_div[,ix]==1,1]))
type_matrix=matrix(nrow=length(genes),ncol=3,dimnames = list(genes,c("TEvsAG","TEvsOV","AGvsOV")))
type_matrix[,"TEvsAG"]=dy_ta[genes,3]
type_matrix[,"TEvsOV"]=dy_to[genes,3]
type_matrix[,"AGvsOV"]=dy_ao[genes,3]
type_matrix[which(is.na(type_matrix),arr.ind=T)]=0
type_matrix[which(type_matrix==2,arr.ind=T)]=-1
# p=pheatmap(type_matrix,show_rownames = F,color = colorRampPalette(c('#56B4E9', "white", '#E69F00'))(10),cutree_rows = 5)
ix = apply(type_matrix,1,function(x){
  ix = ifelse(!any(x==-1),1,ifelse(!any(x==1),-1,0))
  if(ix==1){
    if(sum(x)==1){return(10)}else{return(11)}
  }
  if(ix==-1){
    if(sum(x)==-1){return(-30)}else{return(-31)}
  }
  if(ix==0){
    #browser()
    if(x[1]==1){return(21)}
    if(x[2]==1){return(22)}
    if(x[3]==1){return(23)}
  }
  return(ix)
})
table(ix)
cNo$dya=melt(table(ix))[,2]
m1 = Heatmap(type_matrix[ix==11,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m2 = Heatmap(type_matrix[ix==10,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=3,show_row_dend = F,gap=unit(0,"npc"))
m3 = Heatmap(type_matrix[ix==21,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m4 = Heatmap(type_matrix[ix==22,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m5 = Heatmap(type_matrix[ix==23,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m6 = Heatmap(type_matrix[ix==-30,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m7 = Heatmap(type_matrix[ix==-31,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m = m1%v%m2%v%m3%v%m4%v%m5%v%m6%v%m7

pdf("dy.pheatmap.type.transit.pdf",width = 2,height = 6)
draw(m)
dev.off()
save(m1,m2,m3,m4,m5,m6,m7,ix,type_matrix,file="dya_transit.RData")

## dsi
ix=5
genes = unique(c(ta_cons[rowSums(ta_cons[,2:5])==1&ta_cons[,ix]==1,1],ta_div[rowSums(ta_div[,2:5])==1&ta_div[,ix]==1,1],to_cons[rowSums(to_cons[,2:5])==1&to_cons[,ix]==1,1],to_div[rowSums(to_div[,2:5])==1&to_div[,ix]==1,1],ao_cons[rowSums(ao_cons[,2:5])==1&ao_cons[,ix]==1,1],ao_div[rowSums(ao_div[,2:5])==1&ao_div[,ix]==1,1]))
type_matrix=matrix(nrow=length(genes),ncol=3,dimnames = list(genes,c("TEvsAG","TEvsOV","AGvsOV")))
type_matrix[,"TEvsAG"]=ds_ta[genes,3]
type_matrix[,"TEvsOV"]=ds_to[genes,3]
type_matrix[,"AGvsOV"]=ds_ao[genes,3]
type_matrix[which(is.na(type_matrix),arr.ind=T)]=0
type_matrix[which(type_matrix==2,arr.ind=T)]=-1
# p=pheatmap(type_matrix,show_rownames = F,color = colorRampPalette(c('#56B4E9', "white", '#E69F00'))(10),cutree_rows = 5)
ix = apply(type_matrix,1,function(x){
  ix = ifelse(!any(x==-1),1,ifelse(!any(x==1),-1,0))
  if(ix==1){
    if(sum(x)==1){return(10)}else{return(11)}
  }
  if(ix==-1){
    if(sum(x)==-1){return(-30)}else{return(-31)}
  }
  if(ix==0){
    #browser()
    if(x[1]==1){return(21)}
    if(x[2]==1){return(22)}
    if(x[3]==1){return(23)}
  }
  return(ix)
})
table(ix)
cNo$dsi=melt(table(ix))[,2]
m1 = Heatmap(type_matrix[ix==11,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m2 = Heatmap(type_matrix[ix==10,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=3,show_row_dend = F,gap=unit(0,"npc"))
m3 = Heatmap(type_matrix[ix==21,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m4 = Heatmap(type_matrix[ix==22,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m5 = Heatmap(type_matrix[ix==23,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m6 = Heatmap(type_matrix[ix==-30,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m7 = Heatmap(type_matrix[ix==-31,],name = "Transition", show_row_names = FALSE, show_column_names = TRUE,col = colors, cluster_columns = F,row_km=1,show_row_dend = F)
m = m1%v%m2%v%m3%v%m4%v%m5%v%m6%v%m7

pdf("ds.pheatmap.type.transit.pdf",width = 2,height = 6)
draw(m)
dev.off()
save(m1,m2,m3,m4,m5,m6,m7,ix,type_matrix,file="dsi_transit.RData")

#GO


#cols2 <- brewer.pal(9, "Set1")[c(1,2,4)]
cols1 = c('#E69F00', '#56B4E9',brewer.pal(9, "Set3")[3])
names(cols1)=c("1","3","2")
cNo1 = rbind(colSums(cNo[1:2,]),colSums(cNo[3:4,]),colSums(cNo[5:7,]))
rownames(cNo1)=c("3","1","2")
gg=melt(sweep(cNo1,2,colSums(cNo1),"/"))
p <- ggplot(gg) + geom_bar(aes(x=X2,y=value,fill=as.factor(X1)),stat="identity",color="white") + scale_fill_manual(values=cols1)+ theme_classic(12) + xlab("") + ylab("Proportion of genes")
ggsave("proportion_three_types_by_species.pdf",width = 3.1,height = 3)

cNo2 = as.matrix(sweep(cNo[5:7,],2,colSums(cNo[5:7,]),"/"))
rownames(cNo2)=c("TEvsAG","TEvsOV","AGvsOV")
gg = melt(cNo2)
cols2 <- brewer.pal(9, "Set1")[c(1,2,4)]
p <- ggplot(gg) + geom_bar(aes(x=X2,y=value,fill=X1),stat="identity",color="white") + scale_fill_manual(values=cols2)+ theme_classic(12) + xlab("") + ylab("Proportion of genes")
ggsave("proportion_three_comparisons_by_species.pdf",width = 3,height = 3)
############################ Functional annotation
library("clusterProfiler")
library("org.Dm.eg.db")
library("ggplot2")
xx1 <- as.list(org.Dm.egFLYBASE2EG)
xx2 <- as.list(org.Dm.egFLYBASECG)
#### GO enrichment
go_enrich <- function(dme,outfile){
  # browser()
  genelist = as.character(xx1[dme[,1]])
  ego_all=data.frame()
  for (i in unique(dme[,2])){
    genes = as.character(xx1[dme[dme[,2]==i,1]])
    if(length(genes)<=3){next}
    cat(i,length(genes),"\n")
    ego <- enrichGO(gene          = genes,
                    universe      = genelist,
                    OrgDb         = org.Dm.eg.db,
                    ont           = "BP",
                    pAdjustMethod = "fdr",
                    pvalueCutoff  = 1,
                    qvalueCutoff  = 1,
                    readable      = TRUE)
    if(nrow(ego)!=0&&!is.null(ego)){
      ego_all = rbind(ego_all,data.frame(ego,moduleColor=i))
    }
  }
  if(nrow(ego_all)==0){return("no enrichment!")}
  ego_all$FoldEnrichment = apply(ego_all[,c("GeneRatio","BgRatio")],1,function(x){
    #browser()
    x1 = as.numeric(unlist(strsplit(x[1],"/")))
    x2 = as.numeric(unlist(strsplit(x[2],"/")))
    r = (x1[1]/x1[2])/(x2[1]/x2[2])
  })
  write.table(ego_all,outfile,row.names = F, col.names = T, sep="\t",quote=F)
}

dme=subset(rcl.s,TE%in%c(6))[,c(1,3)];outfile="TE.State.GO.enrich.txt"
go_enrich(dme,outfile)
dme=subset(rcl.s,OV%in%c(4,5,6))[,c(1,3)];outfile="OV.State.GO.enrich.txt"
go_enrich(dme,outfile)
dme=subset(rcl.s,AG%in%c(4,5,6))[,c(1,3)];outfile="AG.State.GO.enrich.txt"
go_enrich(dme,outfile)

go = subset(rbind(data.frame(read.table("TE.State.GO.enrich.txt",sep="\t",quote="",header = T),organ="TE"),
            data.frame(read.table("OV.State.GO.enrich.txt",sep="\t",quote="",header = T),organ="OV"),
data.frame(read.table("AG.State.GO.enrich.txt",sep="\t",quote="",header = T),organ="AG")),qvalue<=0.05)

p<- ggplot(subset(go,moduleColor=="-100")) + geom_point(aes(x=Description,y=-log10(qvalue),color=as.character(moduleColor),size=FoldEnrichment),stat="identity")  +  theme(axis.text.x=element_text(angle=0,hjust=1,size=1.5), plot.title = element_text(hjust = 0.5))+ coord_flip() + scale_color_manual(values = cols)+ theme_minimal(10)+facet_wrap(~interaction(organ,moduleColor),scales = "free_y",ncol=3)
ggsave("State-100_GOenrich.pdf",width = 12,height = 2)
p<- ggplot(subset(go,moduleColor=="101")) + geom_point(aes(x=Description,y=-log10(qvalue),color=as.character(moduleColor),size=FoldEnrichment),stat="identity")  +  theme(axis.text.x=element_text(angle=0,hjust=1,size=1.5), plot.title = element_text(hjust = 0.5))+ coord_flip() + scale_color_manual(values = cols)+ theme_minimal(10)+facet_wrap(~interaction(organ,moduleColor),scales = "free_y",ncol=2)
ggsave("State101_GOenrich.pdf",width = 12,height = 3)


##################################
total_genes= 9339
tevsag = read.table("dps_TEvsAG.ConsDiv.txt")
tevsag$te = te[rownames(tevsag),1]
tevsag$ag = ag[rownames(tevsag),1]

x1 = data.frame(table(tevsag[,c(3,5)]),tissue="TE")
colnames(x1)=c("type","cluster","Freq","tissue")
#x1$pro = x1$Freq/sum(x1$Freq)
x2 = data.frame(table(tevsag[,c(3,6)]),tissue="AG")
colnames(x2)=c("type","cluster","Freq","tissue")
#x2$pro = x2$Freq/sum(x2$Freq)
x=rbind(x1,x2)
ctot = tapply(x[,3],x[,c(1,4)],sum)
ctot = melt(ctot);colnames(ctot)=c("type","tissue","Freq");ctot$cluster="0"
x=rbind(ctot,x)
ntot = tapply(x[,3],x[,c(2,4)],sum)

x$Prop = apply(x,1,function(a){
 # browser()
  a1=as.numeric(a[3])/ntot[a[2],a[4]]
  #a2=ntot[a[2],a[4]]/total_genes
  return(a1)
})
#x$lOR = log2(x$Prop)

p <- ggplot(data=x,aes(x=cluster,y=Prop,fill=as.character(type))) + geom_bar(stat="identity") + facet_grid(~tissue)+ theme_classic(10) + scale_fill_manual(values=c("grey",'#E69F00', '#56B4E9'))+ xlab("Clusters") + ylab("Proportion of genes")+ theme(legend.position="none")+ labs(linetype=NULL)
ggsave("RatioOfConsDivbyCluster_TEvsAG.pdf",width=2.5,height=2.5)

tevsov = read.table("dps_TEvsOV.ConsDiv.txt")
tevsov$te = te[rownames(tevsov),1]
tevsov$ov = ov[rownames(tevsov),1]

x1 = data.frame(table(tevsov[,c(3,5)]),tissue="TE")
colnames(x1)=c("type","cluster","Freq","tissue")
#x1$pro = x1$Freq/sum(x1$Freq)
x2 = data.frame(table(tevsov[,c(3,6)]),tissue="OV")
colnames(x2)=c("type","cluster","Freq","tissue")
#x2$pro = x2$Freq/sum(x2$Freq)
x=rbind(x1,x2)
ctot = tapply(x[,3],x[,c(1,4)],sum)
ctot = melt(ctot);colnames(ctot)=c("type","tissue","Freq");ctot$cluster="0"
x=rbind(ctot,x)
ntot = tapply(x[,3],x[,c(2,4)],sum)
x$Prop = apply(x,1,function(a){
  # browser()
  a1=as.numeric(a[3])/ntot[a[2],a[4]]
  #a2=ntot[a[2],a[4]]/total_genes
  return(a1)
})
#x$lOR = log2(x$Prop)

p <- ggplot(data=x,aes(x=cluster,y=Prop,fill=as.character(type))) + geom_bar(stat="identity") + facet_grid(~tissue)+ theme_classic(10) + scale_fill_manual(values=c("grey",'#E69F00', '#56B4E9'))+ xlab("Clusters") + ylab("Proportion of genes")+ theme(legend.position="none")+ labs(linetype=NULL) 
ggsave("RatioOfConsDivbyCluster_TEvsOV.pdf",width=2.5,height=2.5)


agvsov = read.table("dps_AGvsOV.ConsDiv.txt")
agvsov$ag = ag[rownames(agvsov),1]
agvsov$ov = ov[rownames(agvsov),1]

x1 = data.frame(table(agvsov[,c(3,5)]),tissue="AG")
colnames(x1)=c("type","cluster","Freq","tissue")
#x1$pro = x1$Freq/sum(x1$Freq)
x2 = data.frame(table(agvsov[,c(3,6)]),tissue="OV")
colnames(x2)=c("type","cluster","Freq","tissue")
#x2$pro = x2$Freq/sum(x2$Freq)
x=rbind(x1,x2)
ctot = tapply(x[,3],x[,c(1,4)],sum)
ctot = melt(ctot);colnames(ctot)=c("type","tissue","Freq");ctot$cluster="0"
x=rbind(ctot,x)
ntot = tapply(x[,3],x[,c(2,4)],sum)
x$Prop = apply(x,1,function(a){
  # browser()
  a1=as.numeric(a[3])/ntot[a[2],a[4]]
  #a2=ntot[a[2],a[4]]/total_genes
  return(a1)
})
#x$lOR = log2(x$Prop)

p <- ggplot(data=x,aes(x=cluster,y=Prop,fill=as.character(type))) + geom_bar(stat="identity") + facet_grid(~tissue)+ theme_classic(10) + scale_fill_manual(values=c("grey",'#E69F00', '#56B4E9'))+ xlab("Clusters") + ylab("Proportion of genes")+ theme(legend.position="none")+ labs(linetype=NULL) 
ggsave("RatioOfConsDivbyCluster_AGvsOV.pdf",width=2.5,height=2.5)







#######################
agvsov = read.table("dps_AGvsOV.ConsDiv.txt")
agvsov$ag = ag[rownames(agvsov),1]
agvsov$ov = ov[rownames(agvsov),1]

x1 = data.frame(table(agvsov[,c(3,5)]),tissue="AG")
colnames(x1)=c("type","cluster","Freq","tissue")
#x1$pro = x1$Freq/sum(x1$Freq)
x2 = data.frame(table(agvsov[,c(3,6)]),tissue="OV")
colnames(x2)=c("type","cluster","Freq","tissue")
#x2$pro = x2$Freq/sum(x2$Freq)
x=rbind(x1,x2)
ctot = tapply(x[,3],x[,c(2,4)],sum)
x$Prop = apply(x,1,function(a){
  #browser()
  as.numeric(a[3])/ctot[a[2],a[4]]
})
p <- ggplot(data=x,aes(x=cluster,y=Prop,fill=as.character(type))) + geom_bar(stat="identity") + facet_grid(~tissue)+ theme_classic(10) + scale_fill_manual(values=c("grey",'#E69F00', '#56B4E9'))+ xlab("Clusters") + ylab("Ratio")+ theme(legend.position="none")+ labs(linetype=NULL) + scale_y_continuous(limits=c(0,0.35))
ggsave("ProportionOfConsDivbyCluster_AGvsOV.pdf",width=2,height=2.5)

tevsov = read.table("dps_TEvsOV.ConsDiv.txt")
tevsov$te = te[rownames(tevsov),1]
tevsov$ov = ov[rownames(tevsov),1]

x1 = data.frame(table(tevsov[,c(3,5)]),tissue="TE")
colnames(x1)=c("type","cluster","Freq","tissue")
#x1$pro = x1$Freq/sum(x1$Freq)
x2 = data.frame(table(tevsov[,c(3,6)]),tissue="OV")
colnames(x2)=c("type","cluster","Freq","tissue")
#x2$pro = x2$Freq/sum(x2$Freq)
x=rbind(x1,x2)
ctot = tapply(x[,3],x[,c(1,4)],sum)

x$Prop = apply(x,1,function(a){
  #browser()
  as.numeric(a[3])/ctot[a[1],a[4]]
})
p <- ggplot(data=x,aes(x=cluster,y=Prop,fill=as.character(type))) + geom_bar(stat="identity", position=position_dodge()) + facet_grid(type~tissue)+ theme_classic(10) + scale_fill_manual(values=c("grey",'#E69F00', '#56B4E9'))+ xlab("Clusters") + ylab("Ratio")+ theme(legend.position="none")+ labs(linetype=NULL) + scale_y_continuous(limits=c(0,0.35))
ggsave("ProportionOfConsDivbyCluster_TEvsOV.pdf",width=2,height=2.5)
