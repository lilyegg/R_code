#' 利用指定基因进行一致性聚类，返回分类结果
#'
#' @param exp :表达谱数据
#' @param genelist ：基因列表，如m6A
#' @param pItem :0.8
#' @param pFeature :0.8
#' @param clusterAlg ：pam,hc
#' @param distance :euclidean,spearman,pearson
#'
#' @return：cluster1<-conClust[[k]]$consensusClass
#' @export
#'
#' @examples
#' cluster<-consensus(exp,genelist,pItem = 0.8,pFeature = 0.8,clusterAlg = "pam",distance = 'spearman')
consensus<-function(exp,genelist,pItem = 0.8,pFeature = 0.8,clusterAlg = "pam",distance = 'spearman'){
    library(ConsensusClusterPlus)
    getOptK <- function(conClust, minCls = 2, maxCls = 10) {
    #   最佳分类数
    Kvec = minCls: maxCls
    x1 = 0.1
    x2 = 0.9 # threshold defining the intermediate sub-interval

    PAC = rep(NA, length(Kvec))
    names(PAC) = paste("K=", Kvec, sep = "") # from 2 to maxK
    for(i in Kvec){
       M = conClust[[i]][["consensusMatrix"]]
       Fn = ecdf(M[lower.tri(M)])
       PAC[i-1] = Fn(x2) - Fn(x1)
    }
    optK = Kvec[which.min(PAC)]
    return(optK)
    }

   emt<-exp
   conClust <- ConsensusClusterPlus(
     as.matrix(emt[which(rownames(emt) %in% genelist),]),
     maxK = 5,
     reps = 1000,
     pItem = pItem,
     pFeature = pFeature,
     clusterAlg = clusterAlg,
     distance = distance,
     corUse = "complete.obs",
     seed = 123456,
     plot = "png",
     title = "./final",
     finalLinkage = 'ward.D',
     innerLinkage = 'ward.D',
     writeTable = FALSE
    )
    k<-getOptK(conClust, minCls = 2, maxCls = 5)
    message(sprintf("最佳分类数：%d", k))
    if(k>2){k<-2}
    cluster1<-conClust[[k]]$consensusClass
    return(cluster1)
}



#' 绘制分组样本的KM曲线
#'
#' @param cluster ：分组信息，第一列样本，第二列分组
#' @param clinFile ：临床信息，第一列样本，后面列包含time_os（months），status_os
#' @param fileName :KM曲线图名称
#' @param colors :线条配色
#'
#' @return：3张KM曲线图
#' @export
#'
#' @examples：getKM(cluster,clinFile,fileName="Figure1",colors=colorRampPalette(brewer.pal(5,"Set3")))
getKM<-function(cluster,clinFile,fileName,colors){
  library(survival, quietly = TRUE)
  library(survminer, quietly = TRUE)
  library(RColorBrewer)
  library(colorspace)

  clin1<-read.table(file=clinFile, header=TRUE,sep="\t",as.is=TRUE)
  samples<-intersect(cluster[,1],clin1[,1])
  group<-cluster[match(samples,cluster[,1]),]
  clin1<-clin1[match(samples,clin1[,1]),]
  survm6AclusterData<-data.frame()
  survm6AclusterData <- data.frame(
    time = as.numeric(clin1$time_os),
    status = as.numeric(clin1$status_os),
    group = group[,2]
  )
  #rownames(survm6AclusterData)<-clin1[,1]
  kmplot<-survfit(Surv(time,status)~group,data=survm6AclusterData)
  surv_diff <- survdiff(Surv(time,status)~group,data=survm6AclusterData)
  p.val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  message(sprintf("pvalue:%f",p.val))

  survm6Acluster <- ggsurvplot(
    kmplot,
    font.legend=14,xlab="OS Time (months)",legend.title="cluster",pval=signif(p.val,2),palette = colors,
    font.y=c(20),font.x=c(20),font.tickslab =c(18),pval.size=7,risk.table = TRUE,
    ggtheme = theme_survminer()
  )


  pdf(paste0(fileName,".pdf"), width = 10, height = 8)
  print(survm6Acluster,newpage=FALSE)
  dev.off()
  tiff(paste0(fileName,"-72ppi.tif"), width = 10, height = 8, units = "in", res = 72)
  print(survm6Acluster)
  dev.off()
  tiff(paste0(fileName,"-300ppi.tif"), width = 10, height = 8, units = "in", res = 300)
  print(survm6Acluster)
  dev.off()
}



