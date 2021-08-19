

#' ggplot2绘制小提琴图
#'
#' @param gg 两列 type vaule(any names)
#' @param ylab y轴标签
#' @param fileName 图的名字
#'
#' @return 图
#' @export
#'
#' @examples violin_plot(gg,ylab,fileName)
violin_plot<-function(gg,ylab,fileName){
  library("ggplot2")
  library("ggpubr")
  levels<-names(table(gg$type))
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
  ###########箱式图
  p<-ggplot(gg, aes_string(x = "type", y = colnames(gg)[2],fill="type"))+ geom_violin(trim = FALSE)+ geom_boxplot(width = 0.2)+
    scale_fill_manual(values=c( "#E69F00", "#56B4E9"))+theme_classic()+theme(text=element_text(size=15))+
    stat_compare_means(comparisons=my_comparisons,label = "p.signif",method = "t.test")+
    ylab(ylab)+xlab("")+labs(title = colnames(gg)[2])+theme(legend.position = 'none')

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

#' 相似性散点图
#'
#' @param x x轴数据
#' @param y y轴数据
#' @param fileName 文件名
#' @param title 图片上展示的main title
#'
#' @return 图
#' @export
#'
#' @examples ggcors(x,y,fileName,title)
ggcors<-function(x,y,fileName,title){
  library(ggplot2)
  gg<-data.frame(x=x,y=y)
  p2<-ggplot(data =gg, mapping = aes(x =x, y = y)) + geom_point(alpha=0.8,colour="#FC8D62")+
    geom_smooth(method="lm",se=FALSE)+theme_bw()+
    stat_cor(data=gg,method="pearson")+xlab("Methylation value")+ylab("Expression value")+
    theme(text=element_text(size=15))+labs(title = title)
  pdf(paste0(fileName,".pdf"), width = 8, height = 8)
  print(p2)
  dev.off()
  tiff(paste0(fileName,"-72ppi.tif"), width = 8, height = 8, units = "in", res = 72)
  print(p2)
  dev.off()
  tiff(paste0(fileName,"-300ppi.tif"), width = 8, height = 8, units = "in", res = 300)
  print(p2)
  dev.off()
}

