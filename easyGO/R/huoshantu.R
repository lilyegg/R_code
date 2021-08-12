#' Title
#'
#' @param nrDEG : limma差异基因分析结果，增加change一列
#' @param fileName ：输出文件名字
#' @param width ：宽度
#' @param height ：长度
#' @param label :火山图上是否标记top10的基因
#'
#' @return 返回一张带标签的火山图
#' @export
#'
#' @examples
huoshantu<-function(nrDEG,fileName,width,height,label=TRUE){
  library(ggplot2)
  library(ggrepel)
  p <- ggplot(nrDEG, aes(x = logFC, y = -log10(adj.P.Val), color = change)) +
    geom_point(alpha = 0.6, size = 1.2) + theme_bw(base_size = 15) +
    scale_color_manual(name = "", values = c("red", "green", "black"), limits = c("up", "down", "not-signficant"))
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



