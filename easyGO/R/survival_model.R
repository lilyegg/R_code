#' 对候选基因进行单因素cox分析
#'
#' @param exp 表达谱
#' @param candidate 候选基因
#' @param clinFile 临床信息文件，包含sample time_os status_os
#' @param pvalue 0.05
#' @param fileName 输出结果文件
#'
#' @return 筛选的与预后显著相关的基因
#' @export
#'
#' @examples
#' singleCox(exp,candidate,clinFile,pvalue,fileName)
singleCox<-function(exp,candidate,clinFile,pvalue,fileName){
  library(survival, quietly = TRUE)
  library(survminer, quietly = TRUE)
  candidate<-candidate
  mrna_datexpr_cancer<-exp
  clin<-read.table(file=clinFile, header=TRUE,sep="\t",as.is=TRUE)

  gg<-t(mrna_datexpr_cancer[which(rownames(mrna_datexpr_cancer) %in% candidate),])
  exp<-data.frame(sample=rownames(gg),gg)
  merged<-merge(clin,exp,by="sample")
  merged_data<-merged[,-1]
  rownames(merged_data)<-merged[,1]
  colnames(merged_data)<-gsub("[.]","-",colnames(merged_data))

  coxResult<-vector()
  tt<-c()
  for(i in 1:length(candidate)){
    value<-merged_data[,which(colnames(merged_data) %in% candidate[i])]
    if(sum(value)>0){
      tt<-c(tt,candidate[i])
      group=ifelse(value>median(value,na.rm=TRUE),1,0)
      dat <- list(time=merged_data$time_os,
                  status=merged_data$status_os,
                  Risk=group)
      fmla <- as.formula("Surv(time, status) ~Risk")
      cox <- coxph(fmla, data = dat)
      re=c(summary(cox)[[7]][5],summary(cox)[[7]][2],summary(cox)[[8]][3],summary(cox)[[8]][4])
      coxResult<-rbind(coxResult,re)
    }
  }
  colnames(coxResult)<-c('p_value','HR','Low 95%CI','High 95%CI')
  rownames(coxResult)<-tt
  write.table(coxResult,file=paste0(fileName,"_coxResult.txt"),quote=FALSE,sep="\t",row.names=TRUE,col.names=TRUE)
  cox_gene<-rownames(coxResult)[which(coxResult[,1]<pvalue)]
  message(sprintf("cox_gene:%d",length(cox_gene)))
  result<-list(coxResult=coxResult,cox_gene=cox_gene)
  return(result)
}


#' 根据单个基因表达谱绘制KM曲线
#'
#' @param exp 表达谱
#' @param plotgene 要画图的基因
#' @param clinFile 临床信息文件，包含sample time_os status_os
#' @param fileName 图片文件名
#' @param color 图片中曲线颜色，有默认值
#'
#' @return p.val
#' @export
#'
#' @examples
#' getKMgene(exp,plotgene,clinFile,fileName,color=c("#E7B800","#2E9FDF"))
getKMgene<-function(exp,plotgene,clinFile,fileName,color=c("#E7B800","#2E9FDF")){
  library(survival, quietly = TRUE)
  library(survminer, quietly = TRUE)
  mrna_datexpr_cancer<-exp
  clin<-read.table(file=clinFile, header=TRUE,sep="\t",as.is=TRUE)

  gg<-t(mrna_datexpr_cancer)
  exp<-data.frame(sample=rownames(gg),gg)
  merged<-merge(clin,exp,by="sample")
  merged_data<-merged[,-1]
  rownames(merged_data)<-merged[,1]
  colnames(merged_data)<-gsub("[.]","-",colnames(merged_data))


  index<-which(colnames(merged_data) %in% plotgene)
  value<-as.numeric(merged_data[,index])
  Risk<-ifelse(value>median(value,na.rm=TRUE),"High_exp","Low_exp")
  test1 <- list(time=merged_data$time_os,
                status=merged_data$status_os,
                Exp=Risk)
  kmplot<-survfit(Surv(time,status)~Exp,data=test1)
  surv_diff <- survdiff(Surv(time,status)~Exp,data=test1)
  p.val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)

  tiff(paste0(fileName,"-72ppi.tif"),width=6,height=6,units="in",res=72)
  j<-ggsurvplot(kmplot, data = test1,conf.int=TRUE,legend.title = paste0(plotgene," (median value ",round(median(value,na.rm=TRUE),2),")"),legend.labs = c("High_exp", "Low_exp"),palette = color,font.legend=14,xlab="OS Time (months)",pval=TRUE,ggtheme = theme_bw())
  print (j)
  dev.off()

  tiff(paste0(fileName,"-300ppi.tif"),width=6,height=6,units="in",res=300)
  j<-ggsurvplot(kmplot, data = test1,conf.int=TRUE,legend.title = paste0(plotgene," (median value ",round(median(value,na.rm=TRUE),2),")"),legend.labs = c("High_exp", "Low_exp"),palette = color,font.legend=14,xlab="OS Time (months)",pval=TRUE,ggtheme = theme_bw())
  print (j)
  dev.off()

  pdf(paste0(fileName,".pdf"),width=6,height=6)
  j<-ggsurvplot(kmplot, data = test1,conf.int=TRUE,legend.title = paste0(plotgene," (median value ",round(median(value,na.rm=TRUE),2),")"),legend.labs = c("High_exp", "Low_exp"),palette = color,font.legend=14,xlab="OS Time (months)",pval=TRUE,ggtheme = theme_bw())
  print (j,newpage=FALSE)
  dev.off()
  return(p.val)

}

#' 对筛选的基因绘制森林图
#'
#' @param coxResult singleCox分析结果
#' @param fileName 图片名称
#' @param xticks x轴区间
#'
#' @return 森林图
#' @export
#'
#' @examples plotForest(coxResult,xticks,fileName)
plotForest<-function(coxResult,xticks,fileName){
  test<-coxResult
  sample<-data.frame(Gene=rownames(test),test,stringsAsFactors=FALSE)
  tabletext1<-as.character(sample[,1]) ##gene name
  tabletext2<-round(as.numeric(sample[,2]),3) ##pvalue
  tabletext3<-paste(round(as.numeric(sample[,3]),2),round(as.numeric(sample[,4]),2),sep="(")
  tabletext4<-paste(tabletext3,round(as.numeric(sample[,5]),2),sep="-")
  tabletext5<-paste0(tabletext4,sep=")") ###HR
  tabletext<-cbind(tabletext1,tabletext2,tabletext5)
  tabletext<-rbind(c("Gene","Pvalue","Hazard Ratio (95% CI)"),c("","",""),tabletext)
  #tabletext<-rbind(c("Gene","Pvalue","Hazard Ratio (95% CI)"),tabletext)
  library("forestplot")
  pdf(paste0(fileName,".pdf"),width=10,height=8)
  forestplot(labeltext=tabletext,
             mean = c(NA,NA,round(sample[,3],2)),##
             lower = c(NA,NA,round(sample[,4],2)),##
             upper = c(NA,NA,round(sample[,5],2)),#
             boxsize = 0.3,##
             graph.pos=2,#
             graphwidth = unit(0.4,"npc"),#
             fn.ci_norm="fpDrawCircleCI",#
             col=fpColors(box="#00BFC4", lines="black", zero = "black"),#
             lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#
             zero=1,#
             lwd.zero=2,#
             grid=T,
             lwd.xaxis=2,#
             title="Hazard Ratio",
             xlab="HR value",#
             clip=c(-Inf,6),
             xticks = xticks,#
             colgap = unit(0.5,"cm"),#
             txt_gp = fpTxtGp(ticks = gpar(cex=1.3),xlab  = gpar(cex = 1.5)),
             new_page = FALSE
  )
  dev.off()
  tiff(paste0(fileName,"-72ppi.tif"),width=10,height=8,units="in",res=72)
  forestplot(labeltext=tabletext,
             mean = c(NA,NA,round(sample[,3],2)),##
             lower = c(NA,NA,round(sample[,4],2)),##
             upper = c(NA,NA,round(sample[,5],2)),#
             boxsize = 0.3,##
             graph.pos=2,#
             graphwidth = unit(0.4,"npc"),#
             fn.ci_norm="fpDrawCircleCI",#
             col=fpColors(box="#00BFC4", lines="black", zero = "black"),#
             lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#
             zero=1,#
             lwd.zero=2,#
             grid=T,
             lwd.xaxis=2,#
             title="Hazard Ratio",
             xlab="HR value",#
             clip=c(-Inf,6),
             xticks = xticks,#
             colgap = unit(0.5,"cm"),#
             txt_gp = fpTxtGp(ticks = gpar(cex=1.3),xlab  = gpar(cex = 1.5)),
             new_page = FALSE
  )
  dev.off()
  tiff(paste0(fileName,"-300ppi.tif"),width=10,height=8,units="in",res=300)
  forestplot(labeltext=tabletext,
             mean = c(NA,NA,round(sample[,3],2)),##
             lower = c(NA,NA,round(sample[,4],2)),##
             upper = c(NA,NA,round(sample[,5],2)),#
             boxsize = 0.3,##
             graph.pos=2,#
             graphwidth = unit(0.4,"npc"),#
             fn.ci_norm="fpDrawCircleCI",#
             col=fpColors(box="#00BFC4", lines="black", zero = "black"),#
             lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#
             zero=1,#
             lwd.zero=2,#
             grid=T,
             lwd.xaxis=2,#
             title="Hazard Ratio",
             xlab="HR value",#
             clip=c(-Inf,6),
             xticks = xticks,#
             colgap = unit(0.5,"cm"),#
             txt_gp = fpTxtGp(ticks = gpar(cex=1.3),xlab  = gpar(cex = 1.5)),
             new_page = FALSE
  )
  dev.off()

}




#' lasso分析
#'
#' @param merged_data
#' @param cox_gene 基因列表
#' @param fileName 图片名称
#' @param fileName1 文件名称
#'
#' @return coef_gene
#' @export
#'
#' @examples lasso(merged_data,cox_gene,fileName,fileName1)
lasso<-function(merged_data,cox_gene,fileName,fileName1){
  library(glmnet)
  set.seed(10101)
  index<-which(merged_data$time_os>0)
  exp_losso <- as.matrix(merged_data[index,cox_gene])
  cli_losso <- as.matrix(data.frame(time=merged_data$time_os[index],status=merged_data$status_os[index]))
  fit <- glmnet(exp_losso,cli_losso, family = 'cox')
  #pdf('./result/Figure/Figure7_A.pdf',width = 8,height = 6)
  #plot(fit, xvar = 'lambda')
  #dev.off()
  cvfit = cv.glmnet(exp_losso,cli_losso, family = 'cox',nfolds=20)
  pdf(paste0(fileName,'_cvfit.pdf'),width = 8,height = 8)
  plot(cvfit)
  dev.off()
  tiff(paste0(fileName,'_cvfit-72ppi.tif'),width = 8,height = 8,units="in",res=72)
  plot(cvfit)
  dev.off()
  tiff(paste0(fileName,'_cvfit-300ppi.tif'),width = 8,height = 8,units="in",res=300)
  plot(cvfit)
  dev.off()
  #基因筛选
  coef.min = coef(cvfit,s='lambda.min')
  coef.min
  coef.min.out <- coef.min[which(coef.min != 0),]
  coef_gene <- data.frame(gene=names(coef.min.out),coef=coef.min.out) # 3 genes
  write.table(coef_gene,paste0(fileName1,"_lasso_coef_gene.txt"),col.names = T,row.names = F,quote = F,sep = '\t')

  pdf(paste0(fileName,'_coef.pdf'),width = 8,height = 8)
  bar<-barplot(coef_gene[,2],names.arg=coef_gene[,1],xlab="",ylab="",col="red",main="coefficients",border="red",horiz=TRUE,width=0.7)
  abline(v=0)
  text(0.1,bar,round(coef_gene[,2],4))
  dev.off()
  tiff(paste0(fileName,'_coef-72ppi.tiff'),width = 8,height = 8,units="in",res=72)
  bar<-barplot(coef_gene[,2],names.arg=coef_gene[,1],xlab="",ylab="",col="red",main="coefficients",border="red",horiz=TRUE,width=0.7)
  abline(v=0)
  text(0.1,bar,round(coef_gene[,2],4))
  dev.off()
  tiff(paste0(fileName,'_coef-300ppi.tiff'),width = 8,height = 8,units="in",res=300)
  bar<-barplot(coef_gene[,2],names.arg=coef_gene[,1],xlab="",ylab="",col="red",main="coefficients",border="red",horiz=TRUE,width=0.7)
  abline(v=0)
  text(0.1,bar,round(coef_gene[,2],4))
  dev.off()
  return(coef_gene)

}


#' 根据风险得分进行分类，然后绘制KM曲线
#'
#' @param merged_data 包含time_os，status_os
#' @param riskScore 第一列样本，第二列风险得分
#' @param fileName1 文件名
#' @param fileName2 图片名
#' @param color KM曲线颜色
#'
#' @return 图和文件,以及dat
#' @export
#'
#' @examples
#' getKMriskscore(merged_data,riskScore,fileName1,fileName2,color=c("firebrick3","skyblue"))
#' dat<-getKMriskscore(merged_data=merged_data,riskScore=riskScore,fileName1="./result/Table/SupplementaryTable14",fileName2="./result/Figure/Figure11D",color=c("firebrick3","skyblue"))
#' plotThree(dat=dat,sigData=sigData,fileName="./result/Figure/Figure11")
#' roc1(dat,marker,fileName)
#' roc2(dat,marker,fileName)

getKMriskscore<-function(merged_data,riskScore,fileName1,fileName2,color=c("firebrick3","skyblue")){
  library("maxstat")
  library(survival, quietly = TRUE)
  library(survminer, quietly = TRUE)
  library(RColorBrewer)
  library(colorspace)
  test1 <- data.frame(time=merged_data$time_os,
                      status=merged_data$status_os,
                      Risk=riskScore[,2])
  res_cut <- surv_cutpoint(test1, time = "time", event = "status",variables = "Risk")
  cutoff<-summary(res_cut)$cutpoint  ## -0.2745616
  message(sprintf("cutoff:%f",cutoff))
  Risk=ifelse(riskScore[,2]>cutoff,'High','Low')   ####
  test1$group<-Risk
  dat<-data.frame(sample=riskScore[,1],test1)
  write.table(cbind(riskScore,Risk),file=paste0(fileName1,"_riskscore_group.txt"),quote=FALSE,sep="\t",row.names=TRUE,col.names=FALSE)

  test1 <- list(time=merged_data$time_os,
                status=merged_data$status_os,
                Risk=Risk)

  kmplot<-survfit(Surv(time,status)~Risk,data=test1)
  surv_diff <- survdiff(Surv(time,status)~Risk,data=test1)
  p.val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  #p<-ggsurvplot(kmplot, data = test1,font.legend=14,xlab="Time (days)",legend.title="Risk",pval=TRUE,font.y=c(20),font.x=c(20),font.tickslab =c(18),pval.size=7,palette = c("#B2DF8A","#A6CEE3"),risk.table = TRUE)
  p<-ggsurvplot(kmplot, data = test1,font.legend=14,xlab="OS (months)",legend.title="Risk",pval=TRUE,font.y=c(20),font.x=c(20),font.tickslab =c(18),
                pval.size=7,palette = color,risk.table = TRUE,ggtheme = theme_bw())

  tiff(paste0(fileName2,"-72ppi.tif"),width=8,height=6,units="in",res=72)
  print(p)
  dev.off()
  tiff(paste0(fileName2,"-300ppi.tif"),width=8,height=6,units="in",res=300)
  print(p)
  dev.off()
  pdf(paste0(fileName2,".pdf"),width=8,height=6)
  print(p,newpage=FALSE)
  dev.off()
  return(dat)
}




#' 绘制三联图
#'
#' @param dat 包括sample time status Risk group
#' @param sigData 表达谱
#' @param fileName 图名
#'
#' @return 三联图ABC
#' @export
#'
#' @examples plotTree(dat,sigData,fileName)
plotThree<-function(dat,sigData,fileName){
  count1<-length(which(dat$group=='Low'))
  count2<-length(which(dat$group=='High'))
  riskScore1<-dat$Risk
  tiff(paste0(fileName,"A-72ppi.tif"),width=8,height=4,units="in",res=72)
  par(mar=c(5,5,4,2))
  plot(1:length(riskScore1),riskScore1[order(riskScore1)],col=c(rep("skyblue",count1),rep("firebrick3",count2)),xaxt="n",xlab = NA,ylab="Risk Score",cex.lab=1.5)
  legend("bottomright",legend=c('Low','High'),col =c("skyblue","firebrick3"),pch=19)
  dev.off()
  tiff(paste0(fileName,"A-300ppi.tif"),width=8,height=4,units="in",res=300)
  par(mar=c(5,5,4,2))
  plot(1:length(riskScore1),riskScore1[order(riskScore1)],col=c(rep("skyblue",count1),rep("firebrick3",count2)),xaxt="n",xlab = NA,ylab="Risk Score",cex.lab=1.5)
  legend("bottomright",legend=c('Low','High'),col =c("skyblue","firebrick3"),pch=19)
  dev.off()
  pdf(paste0(fileName,"A.pdf"),width=8,height=4)
  par(mar=c(5,5,4,2))
  plot(1:length(riskScore1),riskScore1[order(riskScore1)],col=c(rep("skyblue",count1),rep("firebrick3",count2)),xaxt="n",xlab = NA,ylab="Risk Score",cex.lab=1.5)
  legend("bottomright",legend=c('Low','High'),col =c("skyblue","firebrick3"),pch=19)
  dev.off()

  tiff(paste0(fileName,"B-72ppi.tif"),width=8,height=4,units="in",res=72)
  par(mar=c(5,5,4,2))
  plot(1:length(riskScore1),dat[order(riskScore1),"time"],col=ifelse(dat[order(riskScore1),"status"]==1,'firebrick3','skyblue'),pch=19,xaxt="n",xlab = NA,ylab="survival time",cex.lab=1.5)
  legend("topright",legend = c("LIVING","DECEASED"),col = c("skyblue","firebrick3"),pch=19)
  dev.off()
  tiff(paste0(fileName,"B-300ppi.tif"),width=8,height=4,units="in",res=300)
  par(mar=c(5,5,4,2))
  plot(1:length(riskScore1),dat[order(riskScore1),"time"],col=ifelse(dat[order(riskScore1),"status"]==1,'firebrick3','skyblue'),pch=19,xaxt="n",xlab = NA,ylab="survival time",cex.lab=1.5)
  legend("topright",legend = c("LIVING","DECEASED"),col = c("skyblue","firebrick3"),pch=19)
  dev.off()
  pdf(paste0(fileName,"B.pdf"),width=8,height=4)
  par(mar=c(5,5,4,2))
  plot(1:length(riskScore1),dat[order(riskScore1),"time"],col=ifelse(dat[order(riskScore1),"status"]==1,'firebrick3','skyblue'),pch=19,xaxt="n",xlab = NA,ylab="survival time",cex.lab=1.5)
  legend("topright",legend = c("LIVING","DECEASED"),col = c("skyblue","firebrick3"),pch=19)
  dev.off()

  gg<-t(sigData)[,order(riskScore1)]
  min<-round(min(gg))
  max<-round(max(gg))
  bk1=unique(seq(min,0,length=50))
  bk2=unique(seq(0.001,max,length=50))
  bk<-c(bk1,bk2)
  bk<-seq(min,quantile(gg)[4],length=100)
  color = c(colorRampPalette(c("navy","white"))(length(bk1)),colorRampPalette(c("white","firebrick3"))(length(bk2)))
  library(pheatmap)
  par(mar=c(5,5,4,2))
  p<-pheatmap(t(sigData)[,order(riskScore1)],cluster_cols = F,cluster_rows = T,show_colnames = F,breaks=bk,color= colorRampPalette(c("navy", "white", "firebrick3"))(100))
  tiff(paste0(fileName,"C-72ppi.tif"),width=8,height=4,units="in",res=72)
  print(p)
  dev.off()
  tiff(paste0(fileName,"C-300ppi.tif"),width=8,height=4,units="in",res=300)
  print(p)
  dev.off()
  pdf(paste0(fileName,"C.pdf"),width=8,height=4)
  print(p)
  dev.off()
}



#' 绘制ROC曲线
#'
#' @param dat 包括sample time status Risk group
#' @param marker 1或者0，1表示利用风险得分，0表示利用分组
#' @param fileName 图名
#'
#' @return 1，3，5年的roc曲线
#' @export
#'
#' @examples
#' roc1(dat,marker,fileName)
roc1<-function(dat,marker,fileName){
  library("survivalROC")
  if(marker==1){
    marker<-dat$Risk
  }else{marker<-as.numeric(dat$group=="High")}
  dat<-as.data.frame(dat)
  fit <- survivalROC(Stime =dat$time, status = dat$status, marker = marker, predict.time = 1, method = "KM")
  message(sprintf("1 month AUC: %f",fit$AUC) ) ##0.6740844
  optimalCutoff <- fit$cut.values[which.max(fit$TP - fit$FP)]
  ROC.1<-fit

  fit <- survivalROC(Stime =dat$time, status = dat$status, marker = marker, predict.time = 6, method = "KM")
  message(sprintf("6 months AUC: %f",fit$AUC) )
  optimalCutoff <- fit$cut.values[which.max(fit$TP - fit$FP)]
  ROC.2<-fit

  fit <- survivalROC(Stime =dat$time, status = dat$status, marker = marker, predict.time = 12, method = "KM")
  message(sprintf("1 year AUC: %f",fit$AUC) )
  optimalCutoff <- fit$cut.values[which.max(fit$TP - fit$FP)]
  ROC.3<-fit

  fit <- survivalROC(Stime =dat$time, status = dat$status, marker = marker, predict.time = 36, method = "KM")
  message(sprintf("3 years AUC: %f",fit$AUC) )
  optimalCutoff <- fit$cut.values[which.max(fit$TP - fit$FP)]
  ROC.4<-fit

  fit <- survivalROC(Stime =dat$time, status = dat$status, marker = marker, predict.time = 60, method = "KM")
  message(sprintf("5 years AUC: %f",fit$AUC) )
  optimalCutoff <- fit$cut.values[which.max(fit$TP - fit$FP)]
  ROC.5<-fit


  #"#FC8D62","#E78AC3","#8DA0CB","#66C2A5"
  tiff(paste0(fileName,"-72ppi.tif"),width=8,height=6,units="in",res=72)
  par(mar=c(5,5,4,2))
  plot(ROC.3$FP,ROC.3$TP,type='l',xlab='False positive rate',ylab='True postive rate',main='Time-dependent ROC curve',col="#FC8D62",lwd = 2,cex.lab=1.5)
  lines(ROC.4$FP,ROC.4$TP,type='l',col="#E78AC3",lwd=2)
  lines(ROC.5$FP,ROC.5$TP,type='l',col="#8DA0CB",lwd=2)
  abline(0,1)
  legend(0.6,0.25,c(paste0("AUC at 1 year: ",signif(ROC.3$AUC,2)),paste0("AUC at 3 years: ",signif(ROC.4$AUC,2)),paste0("AUC at 5 years: ",signif(ROC.5$AUC,2))),bty="n",lty=c(1,1,1),lwd=c(2,2,2),col=c("#FC8D62","#E78AC3","#8DA0CB","#66C2A5"),text.col =c("#FC8D62","#E78AC3","#8DA0CB","#66C2A5","#A6D854"))
  dev.off()
  tiff(paste0(fileName,"-300ppi.tif"),width=8,height=6,units="in",res=300)
  par(mar=c(5,5,4,2))
  plot(ROC.3$FP,ROC.3$TP,type='l',xlab='False positive rate',ylab='True postive rate',main='Time-dependent ROC curve',col="#FC8D62",lwd = 2,cex.lab=1.5)
  lines(ROC.4$FP,ROC.4$TP,type='l',col="#E78AC3",lwd=2)
  lines(ROC.5$FP,ROC.5$TP,type='l',col="#8DA0CB",lwd=2)
  abline(0,1)
  legend(0.6,0.25,c(paste0("AUC at 1 year: ",signif(ROC.3$AUC,2)),paste0("AUC at 3 years: ",signif(ROC.4$AUC,2)),paste0("AUC at 5 years: ",signif(ROC.5$AUC,2))),bty="n",lty=c(1,1,1),lwd=c(2,2,2),col=c("#FC8D62","#E78AC3","#8DA0CB","#66C2A5"),text.col =c("#FC8D62","#E78AC3","#8DA0CB","#66C2A5","#A6D854"))
  dev.off()
  pdf(paste0(fileName,".pdf"),width=8,height=6)
  par(mar=c(5,5,4,2))
  plot(ROC.3$FP,ROC.3$TP,type='l',xlab='False positive rate',ylab='True postive rate',main='Time-dependent ROC curve',col="#FC8D62",lwd = 2,cex.lab=1.5)
  lines(ROC.4$FP,ROC.4$TP,type='l',col="#E78AC3",lwd=2)
  lines(ROC.5$FP,ROC.5$TP,type='l',col="#8DA0CB",lwd=2)
  abline(0,1)
  legend(0.6,0.25,c(paste0("AUC at 1 year: ",signif(ROC.3$AUC,2)),paste0("AUC at 3 years: ",signif(ROC.4$AUC,2)),paste0("AUC at 5 years: ",signif(ROC.5$AUC,2))),bty="n",lty=c(1,1,1),lwd=c(2,2,2),col=c("#FC8D62","#E78AC3","#8DA0CB","#66C2A5"),text.col =c("#FC8D62","#E78AC3","#8DA0CB","#66C2A5","#A6D854"))
  dev.off()
}


#' 绘制单线条的ROC曲线
#'
#' @param dat 包括sample time status Risk group
#' @param marker "high"或者"low"
#' @param fileName 图名
#'
#' @return ROC曲线
#' @export
#'
#' @examples
#' roc2(dat,marker,fileName)
roc2<-function(dat,marker,fileName){
  library(pROC)
  data<-data.frame(status=as.numeric(dat$status),score=dat$Risk,group=as.numeric(dat$group==marker))
  a<-plot.roc(status~group,data)
  b<-plot.roc(status~score,data)
  maxvalue<-max(a$auc,b$auc)
  message(sprintf("max auc:%f",maxvalue))
  if(maxvalue==a$auc){
    data1<-data.frame(status=as.numeric(dat$status),type=as.numeric(dat$group==marker) )
  }else{
    data1<-data.frame(status=as.numeric(dat$status),type=dat$Risk)
  }
  tiff(paste0(fileName,"-72ppi.tif"),width=8,height=6,units="in",res=72)
  plot.roc(status~type,data1,col="red",main="ROC curve")
  legend(0.5,0.15,paste0("AUC : ",signif(maxvalue,2)),bty="n",lty=c(1,1,1),lwd=c(2,2,2),col=c("#FC8D62","#66C2A5","#A6D854"),text.col =c("#FC8D62","#E78AC3","#8DA0CB","#66C2A5","#A6D854"))
  dev.off()
  tiff(paste0(fileName,"-300ppi.tif"),width=8,height=6,units="in",res=300)
  plot.roc(status~type,data1,col="red",main="ROC curve")
  legend(0.5,0.15,paste0("AUC : ",signif(maxvalue,2)),bty="n",lty=c(1,1,1),lwd=c(2,2,2),col=c("#FC8D62","#66C2A5","#A6D854"),text.col =c("#FC8D62","#E78AC3","#8DA0CB","#66C2A5","#A6D854"))
  dev.off()
  pdf(paste0(fileName,".pdf"),width=8,height=6)
  plot.roc(status~type,data1,col="red",main="ROC curve")
  legend(0.5,0.15,paste0("AUC : ",signif(maxvalue,2)),bty="n",lty=c(1,1,1),lwd=c(2,2,2),col=c("#FC8D62","#66C2A5","#A6D854"),text.col =c("#FC8D62","#E78AC3","#8DA0CB","#66C2A5","#A6D854"))
  dev.off()
}





#' KM曲线临床稳定性分析
#'
#' @param dat 包括sample time status Risk group
#' @param subsamples 临床亚群样本名称列表
#' @param label 临床特征名称"age: >60" "Stage: I-II"
#' @param fileName 图名
#'
#' @return KM曲线图
#' @export
#'
#' @examples getKMsubclin(dat,subsamples,label,fileName)
getKMsubclin<-function(dat,subsamples,label,fileName){
  index<-which(dat[,1] %in% subsamples)
  tmp<-dat[index,]
  test<-list(time=tmp$time,status=tmp$status,score=tmp$Risk,group=tmp$group)
  kmplot<-survfit(Surv(time,status)~group,data=test)
  surv_diff <- survdiff(Surv(time,status)~group,data=test)
  p.val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  p<-ggsurvplot(kmplot, data = test,font.legend=14,xlab="OS (months)",legend.title="group",pval=TRUE,font.y=c(20),font.x=c(20),font.tickslab =c(18),
                pval.size=7,palette = c("#B2DF8A","#A6CEE3"),risk.table = FALSE,ggtheme = theme_bw(),title=label,font.title=c(18))
  tiff(paste0(fileName,"-72ppi.tif"),width=8,height=6,units="in",res=72)
  print(p)
  dev.off()
  tiff(paste0(fileName,"-300ppi.tif"),width=8,height=6,units="in",res=300)
  print(p)
  dev.off()
  pdf(paste0(fileName,".pdf"),width=8,height=6)
  print(p,newpage=FALSE)
  dev.off()
}

#' 风险得分在临床特征间的分布箱线图
#'
#' @param feature 临床特征，Risk group
#' @param clintype 临床特征名字
#' @param fileName 图名
#'
#' @return 图
#' @export
#'
#' @examples clinBox(feature,clintype,fileName)
clinBox<-function(feature,clintype,fileName){
  gg<-feature[,c(clintype,"Risk")]
  index<-which(gg[,1] %in% c(NA,""))
  if(length(index)>0){gg<-gg[-index,]}
  levels<-names(table(gg[,1]))
  my_comparisons<-list()
  k<-1
  for(i in 1:(length(levels)-1)){
    for(j in (i+1):length(levels)){
      if(i!=j){
        my_comparisons[[k]]<-c(levels[i],levels[j])
        k<-k+1
      }
    }
  }
  p<-ggboxplot(gg,x=clintype,y='Risk',fill=clintype)+stat_compare_means(comparisons=my_comparisons,label = "p.signif",method = "wilcox.test")+
    xlab('')+ylab('risk score')+theme(text=element_text(size=15))

  tiff(paste0(fileName,"-72ppi.tif"),width=8,height=6,units="in",res=72)
  print(p)
  dev.off()
  tiff(paste0(fileName,"-300ppi.tif"),width=8,height=6,units="in",res=300)
  print(p)
  dev.off()
  pdf(paste0(fileName,".pdf"),width=8,height=6)
  print(p)
  dev.off()

}






#' 针对单个基因绘制一年的ROC曲线
#'
#' @param dat sample status os Risk group
#' @param marker 1:风险得分 0:Low
#' @param fileName 图片名字
#'
#' @return 图片
#' @export
#'
#' @examples roc3(dat,marker,fileName)
roc3<-function(dat,marker,fileName){
  library("survivalROC")
  if(marker==1){
    marker<-dat$Risk
  }else{marker<-as.numeric(dat$group=="High")}
  dat<-as.data.frame(dat)
  fit <- survivalROC(Stime =dat$time, status = dat$status, marker = marker, predict.time = 1, method = "KM")
  message(sprintf("1 month AUC: %f",fit$AUC) ) ##0.6740844
  optimalCutoff <- fit$cut.values[which.max(fit$TP - fit$FP)]
  ROC.1<-fit

  fit <- survivalROC(Stime =dat$time, status = dat$status, marker = marker, predict.time = 6, method = "KM")
  message(sprintf("6 months AUC: %f",fit$AUC) )
  optimalCutoff <- fit$cut.values[which.max(fit$TP - fit$FP)]
  ROC.2<-fit

  fit <- survivalROC(Stime =dat$time, status = dat$status, marker = marker, predict.time = 12, method = "KM")
  message(sprintf("1 year AUC: %f",fit$AUC) )
  optimalCutoff <- fit$cut.values[which.max(fit$TP - fit$FP)]
  ROC.3<-fit

  fit <- survivalROC(Stime =dat$time, status = dat$status, marker = marker, predict.time = 36, method = "KM")
  message(sprintf("3 years AUC: %f",fit$AUC) )
  optimalCutoff <- fit$cut.values[which.max(fit$TP - fit$FP)]
  ROC.4<-fit

  fit <- survivalROC(Stime =dat$time, status = dat$status, marker = marker, predict.time = 60, method = "KM")
  message(sprintf("5 years AUC: %f",fit$AUC) )
  optimalCutoff <- fit$cut.values[which.max(fit$TP - fit$FP)]
  ROC.5<-fit

  if(ROC.1$AUC>0.6){
    tiff(paste0(fileName,"-72ppi.tif"),width=8,height=6,units="in",res=72)
    par(mar=c(5,5,4,2))
    plot(ROC.1$FP,ROC.1$TP,type='l',xlab='False positive rate',ylab='True postive rate',main='Time-dependent ROC curve',col="#FC8D62",lwd = 2,cex.lab=1.5)
    abline(0,1)
    legend(0.6,0.25,c(paste0("AUC ",signif(ROC.1$AUC,2))),bty="n",lty=c(1,1,1),lwd=c(2,2,2),col=c("#FC8D62","#E78AC3","#8DA0CB","#66C2A5"),text.col =c("#FC8D62","#E78AC3","#8DA0CB","#66C2A5","#A6D854"))
    dev.off()
    tiff(paste0(fileName,"-300ppi.tif"),width=8,height=6,units="in",res=300)
    par(mar=c(5,5,4,2))
    plot(ROC.1$FP,ROC.1$TP,type='l',xlab='False positive rate',ylab='True postive rate',main='Time-dependent ROC curve',col="#FC8D62",lwd = 2,cex.lab=1.5)
    abline(0,1)
    legend(0.6,0.25,c(paste0("AUC ",signif(ROC.1$AUC,2))),bty="n",lty=c(1,1,1),lwd=c(2,2,2),col=c("#FC8D62","#E78AC3","#8DA0CB","#66C2A5"),text.col =c("#FC8D62","#E78AC3","#8DA0CB","#66C2A5","#A6D854"))
    dev.off()
    pdf(paste0(fileName,".pdf"),width=8,height=6)
    par(mar=c(5,5,4,2))
    plot(ROC.1$FP,ROC.1$TP,type='l',xlab='False positive rate',ylab='True postive rate',main='Time-dependent ROC curve',col="#FC8D62",lwd = 2,cex.lab=1.5)
    abline(0,1)
    legend(0.6,0.25,c(paste0("AUC ",signif(ROC.1$AUC,2))),bty="n",lty=c(1,1,1),lwd=c(2,2,2),col=c("#FC8D62","#E78AC3","#8DA0CB","#66C2A5"),text.col =c("#FC8D62","#E78AC3","#8DA0CB","#66C2A5","#A6D854"))
    dev.off()
  }
  return(signif(ROC.1$AUC,2))
}



