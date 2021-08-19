#' WGCNA模块挖掘
#'
#' @param exp 基因表达谱
#' @param genelist 基因列表
#' @param fileName 图片文件名称前缀
#' @param fileName1 文本文件名称前缀
#'
#' @return result<-list(datexpr=datexpr,moduleColors=moduleColors,pow=pow)
#' @export
#'
#' @examples
#' myWGCNA(exp,genelist,fileName="./result/Figure/Figure5",fileName1="./result/Table/SupplementaryTable")
#' mrna_datexpr_cancer<-Data1[[1]][,Data1[[2]][which(Data1[[2]][,2]=="Tumor"),1]]
#' mylist<-unique(c(m6A,cd_m6A_geneset1,cd_m6A_geneset2))
#' aa<-myWGCNA(exp=mrna_datexpr_cancer,genelist=mylist,fileName="./result/Figure/Figure5",fileName1="./result/Table/SupplementaryTable5")


myWGCNA<-function(exp,genelist,fileName,fileName1){
  library(WGCNA)
  exp1<-exp[which(rownames(exp) %in% mylist),]
  datexpr<-as.data.frame(t(exp1))
  gsg = goodSamplesGenes(datexpr, verbose = 3);
  gsg$allOK
  sampleTree = hclust(dist(datexpr), method = "average")
  clust = cutreeStatic(sampleTree, cutHeight = 130, minSize = 10)
  keepSamples = (clust==1 | clust==2)
  #datexpr = datexpr[keepSamples, ]
  nGenes = ncol(datexpr)
  nSamples = nrow(datexpr)
  powers = c(c(1:10), seq(from = 12, to=20, by=1))
  sft = pickSoftThreshold(datexpr, powerVector = powers, verbose = 5)

  pdf(paste0(fileName,"A.pdf"),width=12,height=10)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  abline(h=0.90,col="red")
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()

  tiff(paste0(fileName,"A-72ppi.tif"),width=12,height=10, units = "in", res = 72)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  abline(h=0.90,col="red")
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()

  tiff(paste0(fileName,"A-300ppi.tif"),width=12,height=10, units = "in", res = 300)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  abline(h=0.90,col="red")
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()

  pow<-sft$powerEstimate
  message(sprintf("pow:%d",pow))
  net = blockwiseModules(datexpr, power = pow, maxBlockSize = 7000,
                         TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "FPKM-TOM",
                         verbose = 3)
  mergedColors = labels2colors(net$colors)

  pdf(paste0(fileName,"B.pdf"),width=12,height=10)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      groupLabels = c("Module colors",
                                      "GS.weight"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)

  dev.off()
  tiff(paste0(fileName,"B-72ppi.tif"),width=12,height=10, units = "in", res = 72)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      groupLabels = c("Module colors",
                                      "GS.weight"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  tiff(paste0(fileName,"B-300ppi.tif"),width=12,height=10, units = "in", res = 300)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      groupLabels = c("Module colors",
                                      "GS.weight"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()

  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs
  geneTree = net$dendrograms[[1]]
  genes<-cbind(colnames(datexpr), moduleColors)
  colnames(genes)<-c("gene","module")
  genes<-genes[order(genes[,2]),]
  write.table(genes,file=paste0(fileName1,"_module_gene.txt"),quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
  result<-list(datexpr=datexpr,moduleColors=moduleColors,pow=pow)
  return(result)

}


#' 筛选模块特征
#'
#' @param lastinput ：myWGCNA返回的结果
#' @param traitDt ：临床特征表格
#' @param fileName1 ：文件名称
#' @param fileName2 ：图片名称
#'
#' @return 图片和文本
#' @export
#'
#' @examples
#' selectModule(lastinput,traitDt,fileName1="./result/Table/SupplementaryTable",fileName2="./result/Figure/Figure5")
#' traitDt<-as.matrix(clin_infor[match(rownames(aa$datexpr),rownames(clin_infor)),c(21,22,52,53,82)])
#' traitDt[,ncol(traitDt)]<-plyr::revalue(traitDt[,ncol(traitDt)],c("m6A.clusterA"=1,"m6A.clusterB"=2))
#' selectModule(lastinput=aa,traitDt=traitDt,fileName1="./result/Table/SupplementaryTable5",fileName2="./result/Figure/Figure5")

selectModule<-function(lastinput,traitDt,fileName1,fileName2){
  library(WGCNA)
  MEs0 = moduleEigengenes(lastinput$datexpr, lastinput$moduleColors)$eigengenes
  MEsFemale = orderMEs(MEs0)
  modTraitCor = cor(MEsFemale, traitDt, use = "p")
  nSamples<-nrow(datexpr)
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
  mylist<-list()
  for(i in 1:ncol(traitDt)){
    mylist[[i]] = paste(signif(modTraitCor[,i], 3), "\n(", signif(modTraitP[,i], 1), ")", sep = "")
  }
  textMatrix1<-do.call(cbind,mylist)
  write.table(modTraitCor,file=paste(fileName1,"_WGCNA_module_trait_cor.txt"),quote=FALSE,sep="\t",col.names=NA,row.names=TRUE)
  write.table(modTraitP,file=paste(fileName1,"_WGCNA_module_trait_pvalue.txt"),quote=FALSE,sep="\t",col.names=NA,row.names=TRUE)

  pdf(paste0(fileName2,"C.pdf"),width=12,height=10)
  par(mar=c(12,10,5,3))
  label1 = rownames(modTraitP)
  labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitDt),
                 yLabels = label1,
                 cex.lab = 1,
                 ySymbols = label1, colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix1, setStdMargins = FALSE,
                 cex.text = 1, zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  tiff(paste0(fileName2,"C-72ppi.tif"),width=12,height=10, units = "in", res = 72)
  par(mar=c(12,10,5,3))
  label1 = rownames(modTraitP)
  labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitDt),
                 yLabels = label1,
                 cex.lab = 1,
                 ySymbols = label1, colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix1, setStdMargins = FALSE,
                 cex.text = 1, zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  tiff(paste0(fileName2,"C-300ppi.tif"),width=12,height=10, units = "in", res = 300)
  par(mar=c(12,10,5,3))
  label1 = rownames(modTraitP)
  labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitDt),
                 yLabels = label1,
                 cex.lab = 1,
                 ySymbols = label1, colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix1, setStdMargins = FALSE,
                 cex.text = 1, zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
}





#' 计算基因网显著性
#'
#' @param trait 临床特征，选择其中一列
#' @param lastinput myWGCNA输出结果
#' @param h 虚线高度
#' @param fileName1 文本名字
#' @param fileName2 图片名字
#'
#' @return 图和文本
#' @export
#'
#' @examples
#' GS(trait=traitDt[,3],lastinput,h=0.8,fileName1,fileName1)
GS<-function(trait=traitDt[,3],lastinput,h,fileName1,fileName2){
  moduleColors<-lastinput$moduleColors
  GS1=as.numeric (cor (as.numeric(trait), lastinput$datexpr, use="p"))
  GeneSignificance=abs (GS1)
  ModuleSignificance=tapply (GeneSignificance, moduleColors, mean, na.rm=T)
  geneSignificance<-GeneSignificance
  colors<-moduleColors
  means1 = as.vector(tapply(geneSignificance, colors, mean, na.rm = TRUE))
  se1 = as.vector(tapply(geneSignificance, colors, stdErr))
  pdf(paste0(fileName2,"D.pdf"),width=12,height=10)
  par(las="2",mar=c(6,4,2,2))
  barplot(means1, names.arg = names(table(colors)), ylim=c (0,0.15),col = names(table(colors)), ylab ="Gene Significance", main = "Gene significance across modules")
  addErrorBars(as.vector(means1), as.vector(1.96 * se1), two.side = TRUE)
  abline(h=h,lty = 3) ##figure14
  dev.off()
  tiff(paste0(fileName2,"D-72ppi.tif"),width=12,height=10, units = "in", res = 72)
  par(las="2",mar=c(6,4,2,2))
  barplot(means1, names.arg = names(table(colors)), ylim=c (0,0.15),col = names(table(colors)), ylab ="Gene Significance", main = "Gene significance across modules")
  addErrorBars(as.vector(means1), as.vector(1.96 * se1), two.side = TRUE)
  abline(h=h,lty = 3) ##figure14
  dev.off()
  tiff(paste0(fileName2,"D-300ppi.tif"),width=12,height=10, units = "in", res = 300)
  par(las="2",mar=c(6,4,2,2))
  barplot(means1, names.arg = names(table(colors)), ylim=c (0,0.15),col = names(table(colors)), ylab ="Gene Significance", main = "Gene significance across modules")
  addErrorBars(as.vector(means1), as.vector(1.96 * se1), two.side = TRUE)
  abline(h=h,lty = 3) ##figure14
  dev.off()
  write.table(ModuleSignificance,file=paste(fileName1,"_WGCNA_ModuleSignificance.txt"),quote=FALSE,sep="\t",col.names=NA,row.names=TRUE)
}




#' 筛选模块hub基因
#'
#' @param cut1 横坐标阈值
#' @param cut2 纵坐标阈值
#' @param module 模型名称
#' @param lastinput myWGCNA输出结果
#' @param fileName1 输出图片名称
#'
#' @return 模块hub基因
#' @export
#'
#' @examples
#' getModuleHub(cut1,cut2,module = "turquoise",lastinput,fileName1)
getModuleHub<-function(cut1,cut2,module = "turquoise",lastinput,fileName1){
  cut1<-cut1
  cut2<-cut2
  module = module
  moduleColors<-lastinput$moduleColors
  datexpr<-lastinput$datexpr
  nSamples<-nrow(datexpr)
  MEs0 = moduleEigengenes(datexpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datexpr, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))



  column = match(module, modNames);
  moduleGenes1 = moduleColors==module;
  verboseScatterplot(abs(geneModuleMembership[moduleGenes1, column]),
                     abs(geneTraitSignificance[moduleGenes1, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module,pch=19)
  x = abs(geneModuleMembership[moduleGenes1, column])
  y = abs(geneTraitSignificance[moduleGenes1, 1])
  cor = signif(cor(x,y), 2)
  p<-signif(corPvalueStudent(cor(x,y),nSamples),2)
  main<-paste0("Module membership vs. gene significance\n","cor=",cor,", p=",p)
  pdf(paste0(fileName1,".pdf"),width=8,height=10)
  plot(abs(geneModuleMembership[moduleGenes1, column]),
       abs(geneTraitSignificance[moduleGenes1, 1]), main = main, xlab = paste("Module Membership in", module, "module"),
       ylab = "Gene significance", cex.axis = 1.2, cex.lab = 1.2,
       cex.main = 1.2, col = module,
       pch = 19)
  abline(v=cut1, col = "gray60")
  abline(h=cut2, col = "gray60")
  dev.off()
  tiff(paste0(fileName1,"-72ppi.tif"),width=8,height=10, units = "in", res = 72)
  plot(abs(geneModuleMembership[moduleGenes1, column]),
       abs(geneTraitSignificance[moduleGenes1, 1]), main = main, xlab = paste("Module Membership in", module, "module"),
       ylab = "Gene significance", cex.axis = 1.2, cex.lab = 1.2,
       cex.main = 1.2, col = module,
       pch = 19)
  abline(v=cut1, col = "gray60")
  abline(h=cut2, col = "gray60")
  dev.off()
  tiff(paste0(fileName1,"-300ppi.tif"),width=8,height=10, units = "in", res = 72)
  plot(abs(geneModuleMembership[moduleGenes1, column]),
       abs(geneTraitSignificance[moduleGenes1, 1]), main = main, xlab = paste("Module Membership in", module, "module"),
       ylab = "Gene significance", cex.axis = 1.2, cex.lab = 1.2,
       cex.main = 1.2, col = module,
       pch = 19)
  abline(v=cut1, col = "gray60")
  abline(h=cut2, col = "gray60")
  dev.off()
  x<-abs(geneModuleMembership[moduleGenes1, column])
  y<-abs(geneTraitSignificance[moduleGenes1, 1])
  genes1<-rownames(geneModuleMembership)[moduleGenes1]
  index<-which(x>cut1 & y>cut2)
  hub1<-genes1[index]
  return(result=list(genes=genes1,hub=hub1))
}














