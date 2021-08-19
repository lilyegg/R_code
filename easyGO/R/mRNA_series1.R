#' 处理下载的表达谱数据
#'
#' @param inputfile xena下载的表达谱数据
#' @param gene_map_info_v22 /Users/denglili/Desktop/database/gencode.v22.annotation.gene.probeMap
#'
#' @return list result<-list(exp=mrna_datexpr,pdata=pdata)
#' @export
#'
#' @examples getexp(inputfile="TCGA-LIHC.htseq_fpkm.tsv",gene_map_info_v22=
#' "/Users/denglili/Desktop/database/gencode.v22.annotation.gene.probeMap")
getExp<-function(inputfile,gene_map_info_v22){
	  data<-read.delim(file=inputfile,header=TRUE,sep="\t",as.is=TRUE)
    colnames(data)<-gsub("[.]","-",colnames(data))
    gene_map_info_v22 <- read.delim(gene_map_info_v22,header = T,stringsAsFactors = F)
    TCGA_primary <- merge(data,gene_map_info_v22[,c('id','gene')],by.x = 'Ensembl_ID',by.y = 'id')
    TCGA_primary_gene <- TCGA_primary[,c('gene',names(TCGA_primary)[grepl('TCGA',names(TCGA_primary))])]
    library(limma)
    rt2 = TCGA_primary_gene
    rt2 = as.matrix(rt2)
    rownames(rt2) = rt2[,1]
    exp2 = rt2[,2:ncol(rt2)]
    dimnames2=list(rownames(exp2),colnames(exp2))
    data2 = matrix(as.numeric(as.matrix(exp2)),nrow = nrow(exp2),dimnames = dimnames2)
    mrna_datexpr = avereps(data2)
    data<-read.table("/Users/denglili/Desktop/database/allgene.txt",header=FALSE,sep="\t",as.is=TRUE)
    coding<-data[which(data[,2]=="protein_coding"),7]
    mrna_datexpr<-mrna_datexpr[which(rownames(mrna_datexpr) %in% coding),]
    ########筛选样本
    type<-unlist(lapply(colnames(mrna_datexpr),function(b){strsplit(b,split="-")[[1]][4]}))
    index<-which(type %in% c("01A","11A"))
    mrna_datexpr<-mrna_datexpr[,index]
    type<-unlist(lapply(colnames(mrna_datexpr),function(b){strsplit(b,split="-")[[1]][4]}))
    type<-ifelse(type=="01A","Tumor","Normal")
    pdata<-data.frame(sample=colnames(mrna_datexpr),type=type)
    result<-list(exp=mrna_datexpr,pdata=pdata)
    message(sprintf("Tumor: %d,Normal: %d", length(type[which(type=="Tumor")]),length(type[which(type=="Normal")])))
    return(result)
}

#' 根据表达谱和样本分类信息进行差异基因分析
#'
#' @param exp ：表达谱文件
#' @param pdata ：样本分类信息
#' @param caseLabel ：分组编号如Tumor
#' @param controlLabel ：分组编号如Normal
#' @param fc :差异倍数，如2或者1.5
#' @param pvalue :是否使用pvalue,如果不使用，则pvalue=NULL
#' @param adjP :阈值，一般为0.05
#' @param fileName :nrDEG输出文件名
#'
#' @return :差异基因列表
#' @export
#'
#' @examples getDEGs(exp=result[[1]],pdata=result[[2]],caseLabel="Tumor",controlLabel="Normal"
#' ,fc=1.5,pvalue=NULL,adjP=0.05,fileName="SupplementaryTable_diffgene.txt")
#'
getDEGs<-function(exp,pdata,caseLabel,controlLabel,fc,pvalue=NULL,adjP,fileName){
	  library(limma)
	  samples<-intersect(colnames(exp),pdata[,1])
	  datexpr<-exp[,match(samples,colnames(exp))]
	  pdata<-pdata[match(samples,pdata[,1]),]
	  group_list<-factor(pdata[,2])
    design <- model.matrix(~0+group_list)
    colnames(design) <- levels(group_list)
    rownames(design) <- colnames(datexpr)

    # 差异表达分析 - 创建差异比较矩阵
    contrasts <- paste0(caseLabel, "-", controlLabel)
    contrastMatrix <- makeContrasts(
      contrasts = contrasts,
      levels = levels(group_list)
    )
    fit <- lmFit(datexpr, design)
    fit2 <- contrasts.fit(fit, contrastMatrix)
    fit2 <- eBayes(fit2,0.05)
    tempOutput = topTable(fit2, adjust="fdr", sort.by="B", number=nrow(datexpr))
    nrDEG = na.omit(tempOutput)

    tTag1<-as.data.frame(nrDEG)
    cutoff<-fc
    if(is.null(pvalue)){
    	tTag1$change<-as.factor(ifelse(tTag1$adj.P.Val < adjP & abs(tTag1$logFC) > log(cutoff,base=2),ifelse(tTag1$logFC > log(cutoff,base=2),'up','down'),'not-signficant'))
    }else{
    	tTag1$change<-as.factor(ifelse(tTag1$P.Value < pvalue & abs(tTag1$logFC) > log(cutoff,base=2),ifelse(tTag1$logFC > log(cutoff,base=2),'up','down'),'not-signficant'))
    }
    diffgene<-rownames(tTag1)[which(tTag1$change %in% c("up","down"))]
    message(sprintf("diffgene:%d",length(diffgene)))
    write.table(tTag1,file=fileName,quote=FALSE,sep="\t",col.names=NA,row.names=TRUE)
    return(diffgene)
}


#' 绘制火山图
#'
#' @param nrDEG ：getDEGs生成的差异基因分析文件
#' @param fileName ：火山图文件名
#' @param width ：宽度
#' @param height ：长度
#' @param label ：是否增加top10基因的标签
#'
#' @return
#' @export
#'
#' @examples huoshantu(nrDEG="aa.txt",fileName="Figure12A",width=8,height=8,label=TRUE)
huoshantu<-function(nrDEG,fileName,width,height,label=TRUE){
  library(ggplot2)
  library(ggrepel)
  nrDEG<-read.table(file=nrDEG,header=TRUE,sep="\t",as.is=TRUE)
  p <- ggplot(nrDEG, aes(x = logFC, y = -log10(adj.P.Val), color = change)) +
    geom_point(alpha = 0.6, size = 1.2) + theme_bw(base_size = 15) +
    scale_color_manual(name = "", values = c("#DFC03F","#4376BE", "black"), limits = c("up", "down", "not-signficant"))
  #####加标签
  up <- subset(nrDEG, change == 'up')
  up <- up[order(up$adj.P.Val), ][1:10, ]
  down <- subset(nrDEG, change == 'down')
  down <- down[order(down$adj.P.Val), ][1:10, ]
  data = rbind(up, down)
  data$gene=rownames(data)
  p1 <- p + theme(legend.position = 'right') +
    geom_text_repel(data = rbind(up, down), aes(x = logFC, y = -log10(adj.P.Val), label = gene),
                    size = 3,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE)
  #####保存图片
  if(label==TRUE){
    print(p1)
  }else{
    print(p)
  }
  ggsave(paste0(fileName,'.pdf'),device = 'pdf',width=width,height=height)
  ggsave(paste0(fileName,'-72ppi.tif'),device = 'tiff',width=width,height=height,dpi = 72)
  ggsave(paste0(fileName,'-300ppi.tif'),device = 'tiff',width=width,height=height,dpi = 300)
}


#' 绘制热图
#'
#' @param gg 绘图数据
#' @param anno_col colbar列注释
#' @param ann_colors 列注释的颜色
#' @param bk 分段
#' @param color 热图颜色
#' @param fileName 热图保存文件名
#' @param showrownames 是否显示行名T
#' @param showcolnames  是否显示列名T
#' @param width
#' @param height
#'
#' @return 热图
#' @export
#'
#' @examples phmap(gg,anno_col=NULL,ann_colors=NULL,bk=seq(-2.5,2.5,length=100),color=colorRampPalette(c("navy", "white", "firebrick3"))(100),fileName="./resut/Figure/Figure")
phmap<-function(gg,anno_col=NULL,ann_colors=NULL,bk=seq(-2.5,2.5,length=100),color=colorRampPalette(c("navy", "white", "firebrick3"))(100),fileName,showrownames=T,showcolnames=T,width,height,cellwidth = NA, cellheight=NA){
  library("pheatmap")
  p<-pheatmap(gg,clustering_method = "ward.D",scale="row",breaks=bk,color=color,cluster_cols = F,
              show_rownames =showrownames,show_colnames =showcolnames,fontsize_row=10,annotation_col = anno_col,annotation_colors = ann_colors,cellwidth = cellwidth, cellheight=cellheight)
  pdf(paste0(fileName,".pdf"), width = width, height = height)
  print(p)
  dev.off()
  tiff(paste0(fileName,"-72ppi.tif"), width = width, height = height, units = "in", res = 72)
  print(p)
  dev.off()
  tiff(paste0(fileName,"-300ppi.tif"), width = width, height = height, units = "in", res = 300)
  print(p)
  dev.off()
}

#' PCA分析，不带圈圈
#'
#' @param exp 表达谱数据
#' @param pdata 分组信息
#' @param fileName 图片文件名称
#' @param color 分组颜色，NULL
#'
#' @return PCA图
#' @export
#'
#' @examples
#' pca1(exp,pdata,fileName,color=c("#DFC03F","#4376BE"))
pca1<-function(exp,pdata,fileName,color=c("#DFC03F","#4376BE")){
  library(ggplot2)
  samples<-intersect(colnames(exp),pdata[,1])
  pdata<-pdata[match(samples,pdata[,1]),]
  exp<-exp[,match(samples,colnames(exp))]
  data.pca <- prcomp(t(exp),retx=T,scale=T,center=T)
  a <- summary(data.pca)
  tmp <- a$importance
  pro1 <- as.numeric(sprintf("%.3f",tmp[2,1]))*100
  pro2 <- as.numeric(sprintf("%.3f",tmp[2,2]))*100
  pc = as.data.frame(a$x)
  pc1<-pc$PC1
  pc2<-pc$PC2
  pc3<-pc$PC3
  type<-pdata[,2]
  gg1<-data.frame(x=pc1,y=pc2,z=pc3,group=as.factor(type))
  rownames(gg1)<-rownames(pc)
  if(is.null(color)){color<-c("#DFC03F","#4376BE")}
  a<-ggplot(gg1,aes(x=x,y=y,colour=group))+geom_point()+scale_color_manual(values = c("#DFC03F","#4376BE"))
  b<-a+theme (panel.background=element_blank(), panel.border=element_rect(colour="black", fill=NA))
  c<-b+labs(x = "PCA1", y = "PCA2")
  d<-c+theme(strip.background=element_rect(fill=NA),strip.text=element_text(size=rel(1.2)))
  e<-d+labs(colour="group")
  f<-e+theme_bw()+theme(text=element_text(size=15))

  pdf(paste0(fileName,".pdf"), width = 8, height = 8)
  print(f)
  dev.off()
  tiff(paste0(fileName,"-72ppi.tif"), width = 8, height = 8, units = "in", res = 72)
  print(f)
  dev.off()
  tiff(paste0(fileName,"-300ppi.tif"), width = 8, height = 8, units = "in", res = 300)
  print(f)
  dev.off()
}


#' PCA分析，带圈圈
#'
#' @param exp 表达谱数据
#' @param pdata 分组信息
#' @param fileName 图片文件名称
#'
#' @return PAC图
#' @export
#'
#' @examples
#' pca2(exp,pdata,fileName)
pca2<-function(exp,pdata,fileName){
  library(ggplot2)
  library(ggbiplot)
  samples<-intersect(colnames(exp),pdata[,1])
  pdata<-pdata[match(samples,pdata[,1]),]
  exp<-exp[,match(samples,colnames(exp))]

  data<-t(exp) ###data行是样本列是基因
  data.pca <- prcomp(data,retx=T,scale=T,center=T)
  ff<-ggbiplot(data.pca, obs.scale = 1, var.scale = 1,groups = pdata[,2], ellipse = TRUE,var.axes = F)

  pdf(paste0(fileName,".pdf"), width = 8, height = 8)
  print(ff)
  dev.off()
  tiff(paste0(fileName,"-72ppi.tif"), width = 8, height = 8, units = "in", res = 72)
  print(ff)
  dev.off()
  tiff(paste0(fileName,"-300ppi.tif"), width = 8, height = 8, units = "in", res = 300)
  print(ff)
  dev.off()
}

#' 绘制韦恩图1
#'
#' @param x list列表
#' @param fileName 图片名字
#'
#' @return 韦恩图
#' @export
#'
#' @examples
#' myvenn(x=list1,fileName="./result/Figure/Figure-300ppi.tif")
myvenn<-function(x,fileName){
  library("VennDiagram")
  color<-c("red", "blue", "yellow")[1:length(x)]
  venn.plot <- venn.diagram(x,filename=fileName,
                            resolution = 300,units ="in",
                            height = 10, width = 10,
                            col = "transparent",
                            fill = color,
                            alpha = 0.5,
                            cex = 2.5,#内标签的字体大小
                            fontfamily = "serif",
                            fontface = "bold",
                            cat.default.pos = "outer",#设置标签在圆外面
                            #cat.col = c("darkred", "darkblue", "darkgreen"),
                            cat.cex = 2,#外标签的字体大小
                            cat.fontfamily = "serif",
                            cat.dist = c(0.05, 0.05, 0.05)[1:length(x)],#相对圆圈的位置
                            cat.pos = c(-20,20,180)[1:length(x)]  #相对12点方向旋转的角度)
                            )

}
