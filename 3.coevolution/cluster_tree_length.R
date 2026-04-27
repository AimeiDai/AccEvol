library(ggplot2)
library(ape)
options(stringsAsFactors = F)
setwd("/Users/aimeidai/OneDrive/10.shangrui/remove_batch_effect/")
st = read.table("sample_table_rc.txt",header = T)
expr = read.table("CPM_adjusted_5species.txt")
colnames(expr)=st[,3]


setwd("/Users/aimeidai/OneDrive/10.shangrui/coevolution")
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

tree_len = data.frame()

#### testis
st_te = subset(st,merged_tissue2=="TE")
expr_te = do.call("cbind",tapply(st_te[,3],st_te[,4],function(x){
  #browser()
  if(length(x)==1){
    return(expr[,x])
  }
  y = apply(expr[,x],1,mean)
  return(y)
}))

cte = read.table("TE_cluster.txt")
for (i in unique(cte[,2])){
  ci = subset(cte,V2==i)[,1]
  ei = expr_te[ci,]
  tree_len = rbind(tree_len,data.frame(totTreeLen=cal_tl_bt(ei,times = 1000),cluster=i,tissue="TE"))
}
#tree_len = rbind(tree_len,data.frame(totTreeLen=cal_tl_bt(expr_te,times = 1000),cluster=0,tissue="TE"))


#### ovary
st_te = subset(st,merged_tissue2=="OV")
expr_te = do.call("cbind",tapply(st_te[,3],st_te[,4],function(x){
  #browser()
  if(length(x)==1){
    return(expr[,x])
  }
  y = apply(expr[,x],1,mean)
  return(y)
}))

cte = read.table("OV_cluster.txt")
for (i in unique(cte[,2])){
  ci = subset(cte,V2==i)[,1]
  ei = expr_te[ci,]
  tree_len = rbind(tree_len,data.frame(totTreeLen=cal_tl_bt(ei,times = 1000),cluster=i,tissue="OV"))
}
#tree_len = rbind(tree_len,data.frame(totTreeLen=cal_tl_bt(expr_te,times = 1000),cluster=0,tissue="OV"))


#### accessory gland
st_te = subset(st,merged_tissue2=="AG")
expr_te = do.call("cbind",tapply(st_te[,3],st_te[,4],function(x){
  #browser()
  if(length(x)==1){
    return(expr[,x])
  }
  y = apply(expr[,x],1,mean)
  return(y)
}))

cte = read.table("AG_cluster.txt")
for (i in unique(cte[,2])){
  ci = subset(cte,V2==i)[,1]
  ei = expr_te[ci,]
  tree_len = rbind(tree_len,data.frame(totTreeLen=cal_tl_bt(ei,times = 1000),cluster=i,tissue="AG"))
}
#tree_len = rbind(tree_len,data.frame(totTreeLen=cal_tl_bt(expr_te,times = 1000),cluster=0,tissue="AG"))

#####
p <- ggplot(data=tree_len,aes(x=as.factor(cluster), y= totTreeLen)) + geom_boxplot(outlier.shape = NA) + facet_grid(.~tissue) + theme_bw(12) + theme(axis.text.x=element_text(angle=0,hjust=0.5)) + xlab("Clusters") + ylab("Total tree length")
ggsave("clustersVSTotTreeLen.pdf",width = 4.5, height = 2.5)
by(tree_len[,1],tree_len[,2:3],median)
