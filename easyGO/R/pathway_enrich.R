#' Title
#' 对基因进行功能富集分析
#' @param genelist 输入基因列表
#' @param fileName 输出结果文件
#' @param width 保存图片的宽度
#' @param height 保存图片的高度
#'
#' @return 气泡图和富集结果csv文件
#' @export
#'
#' @examples pathway_enrich(genelit=m6A,fileName="aa",width=10,height=10)
pathway_enrich<-function(genelist,fileName,width,height){
	library(clusterProfiler)
    test1 = bitr(genelist, fromType="SYMBOL",toType="ENTREZID", OrgDb="org.Hs.eg.db")
    ego_MF <- enrichGO(gene = test1$ENTREZID,OrgDb = org.Hs.eg.db,ont = "MF", pAdjustMethod = "BH",readable = TRUE)
    ego_BP <- enrichGO(gene = test1$ENTREZID,OrgDb = org.Hs.eg.db,ont = "BP", pAdjustMethod = "BH",readable = TRUE)
    ego_CC <- enrichGO(gene = test1$ENTREZID,OrgDb = org.Hs.eg.db,ont = "CC", pAdjustMethod = "BH",readable = TRUE)
    kk <- enrichKEGG(gene = test1$ENTREZID,organism = 'hsa')
    write.csv(as.data.frame(ego_MF),paste0(fileName,"_GOMF-enrich.csv"),row.names =FALSE)
    write.csv(as.data.frame(ego_BP),paste0(fileName,"_GOBP-enrich.csv"),row.names =FALSE)
    write.csv(as.data.frame(ego_CC),paste0(fileName,"_GOCC-enrich.csv"),row.names =FALSE)
    write.csv(as.data.frame(kk),paste0(fileName,"_KEGG-enrich.csv"),row.names =FALSE)

    pdf(paste0(fileName,"_GOBP.pdf"),width=width,height=height)
    dotplot(ego_BP,title="EnrichmentGO_BP_dot")#点图，按富集的数从大到小的
    dev.off()
    tiff(paste0(fileName,"_GOBP-72ppi.tiff"),width=width,height=height,units="in",res=72)
    dotplot(ego_BP,title="EnrichmentGO_BP_dot")#点图，按富集的数从大到小的
    dev.off()
    tiff(paste0(fileName,"_GOBP-300ppi.tiff"),width=width,height=height,units="in",res=300)
    dotplot(ego_BP,title="EnrichmentGO_BP_dot")#点图，按富集的数从大到小的
    dev.off()
    pdf(paste0(fileName,"_GOMF.pdf"),width=width,height=height)
    dotplot(ego_MF,title="EnrichmentGO_MF_dot")#点图，按富集的数从大到小的
    dev.off()
    tiff(paste0(fileName,"_GOMF-72ppi.tiff"),width=width,height=height,units="in",res=72)
    dotplot(ego_MF,title="EnrichmentGO_MF_dot")#点图，按富集的数从大到小的
    dev.off()
    tiff(paste0(fileName,"_GOMF-300ppi.tiff"),width=width,height=height,units="in",res=300)
    dotplot(ego_MF,title="EnrichmentGO_MF_dot")#点图，按富集的数从大到小的
    dev.off()
    pdf(paste0(fileName,"_GOCC.pdf"),width=width,height=height)
    dotplot(ego_CC,title="EnrichmentGO_CC_dot")#点图，按富集的数从大到小的
    dev.off()
    tiff(paste0(fileName,"_GOCC-72ppi.tif"),width=width,height=height,units="in",res=72)
    dotplot(ego_CC,title="EnrichmentGO_CC_dot")#点图，按富集的数从大到小的
    dev.off()
    tiff(paste0(fileName,"_GOCC-300ppi.tif"),width=width,height=height,units="in",res=300)
    dotplot(ego_CC,title="EnrichmentGO_CC_dot")#点图，按富集的数从大到小的
    dev.off()
    pdf(paste0(fileName,"_KEGG.pdf"),width=width,height=height)
    dotplot(kk,title="Enrichment_KEGG_dot")#点图，按富集的数从大到小的
    dev.off()
    tiff(paste0(fileName,"_KEGG-72ppi.tiff"),width=width,height=height,units="in",res=72)
    dotplot(kk,title="Enrichment_KEGG_dot")#点图，按富集的数从大到小的
    dev.off()
    tiff(paste0(fileName,"_KEGG-300ppi.tiff"),width=width,height=height,units="in",res=300)
    dotplot(kk,title="Enrichment_KEGG_dot")#点图，按富集的数从大到小的
    dev.off()
}

