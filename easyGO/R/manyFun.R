#' 计算所有基因跟m6A基因的相关性
#'
#' @param cancerExp :表达谱数据
#' @param genelist ：基因列表
#'
#' @return :list(cors=cors,pvalues=pvalues)
#' @export
#'
#' @examples getCors(cancerExp=exp,genelist=m6A)
getCors<-function(cancerExp,genelist){
	pvalues<-vector()
    cors<-vector()
    for(i in 1:length(genelist)){
    tmp<-apply(cancerExp[setdiff(rownames(cancerExp),genelist),],1,function(x){ct<-cor.test(x,cancerExp[genelist[i],],method = "pearson")
        return(c(ct$estimate, ct$p.value))})
    tmp<-t(tmp)
    pvalues<-cbind(pvalues,tmp[,2])
    cors<-cbind(cors,tmp[,1])
    }
    colnames(pvalues)<-genelist
    colnames(cors)<-genelist
    pvalues[is.na(pvalues)]<-1
    cors[is.na(cors)]<-0
    result<-list(cors=cors,pvalues=pvalues)
    return(result)

}


#' 筛选与m6A显著相关的基因
#'
#' @param result getCors返回的结果
#' @param cor 相关性阈值
#' @param cor_percent 比例
#' @param pvalue 显著性阈值
#' @param pvalue_percent 比例
#'
#' @return result<-list(cors=cors,pvalues=pvalues,genelist=cd_m6A_geneset1)
#' @export
#'
#' @examples
#' getCors_sig(result,cor=0.3,cor_percent=0.5,pvalue=0.05,pvalue_percent=1)
getCors_sig<-function(result,cor,cor_percent,pvalue,pvalue_percent){
    cors<-result$cors
    pvalues<-result$pvalues
    #########筛选显著相关基因
    cc<-apply(cors,1,function(b){length(which(abs(b)>cor))})
    candidate1<-names(cc)[which(cc>=ncol(cors)*cor_percent)]
    length(candidate1)
    pvalus_fileter<-pvalues[candidate1,]
    count<-apply(pvalus_fileter,1,function(x){length(which(x<pvalue))})
    cd_m6A_geneset1<-names(count)[which(count>=ncol(pvalues)*pvalue_percent)]
    message(sprintf("get genes : %d",length(cd_m6A_geneset1)))
    result<-list(cors=cors,pvalues=pvalues,genelist=cd_m6A_geneset1)
    return(result)
}


#' risk score与突变负荷、HRD、新抗原负荷和染色体不稳定、干性指数的相关性分析
#'
#' @param cancertype 癌症类型STAD
#' @param fileName1 图片名称
#' @param fileName2 文本名称
#'
#' @return 图片和文本
#' @export
#'
#' @examples hrd(cancertype,fileName1,fileName2)
hrd<-function(cancertype,fileName1,fileName2){
    hrd = system.file('extdata', 'HRD.txt', package = 'easyGO')
    public <- read.table(hrd,header=TRUE,sep="\t",as.is=TRUE)
    public_gbm<-public[which(public[,3]==cancertype),]
    public_gbm[,2]<-paste0(public_gbm[,2],"A")
    ann<-as.data.frame(dat)
    gg1<-merge(ann,public_gbm,by.x="sample",by.y="TCGA.sample.barcode")
    gg1[,1]<-unlist(lapply(gg1[,1],function(b){substr(b,1,15)}))

    snv<-system.file('extdata', 'TCGA_pMHC_SNV_sampleSummary_MC3_v0.2.8.CONTROLLED_170404.tsv', package = 'easyGO')
    SNV_Neoantigen<-read.table(snv,header=TRUE,sep="\t",as.is=TRUE)
    SNV_Neoantigen[,1]<-unlist(lapply(SNV_Neoantigen[,1],function(b){substr(b,1,15)}))

    pubinfor<-merge(gg1,SNV_Neoantigen,by.x="sample",by.y="barcode")

    library("openxlsx")
    si<-system.file('extdata', 'stemindex_mRNAsi_mDNAsi.xlsx', package = 'easyGO')
    mRNAsi<-read.xlsx(si,sheet=1)
    mRNAsi[,1]<-unlist(lapply(mRNAsi[,1],function(b){substr(b,1,15)}))
    mRNAsi<-mRNAsi[which(mRNAsi[,2]==cancertype),]

    pubinfor<-merge(pubinfor,mRNAsi,by.x="sample",by.y="TCGAlong.id")


    gg<-pubinfor[,c(4,5,9,10,46,33,43,44,45,57,60)]
    rownames(gg)<-pubinfor[,1]
    colnames(gg)[10]<-c("SNV_Neoantigen")


    tt<-data.frame(sample=rownames(gg),risk_score=gg[,1],risk_group=gg[,2])
    gg1<-as.matrix(gg)
    w_df1 = reshape2::melt(gg1[,-c(1,2)], varnames = c("sample", "type"))
    final<-merge(w_df1,tt,by="sample")
    write.table(final,file="final.txt",quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)

    final<-read.table(file="final.txt",header=TRUE,sep="\t",as.is=TRUE)
    my_comparisons <- list( c("High", "Low"))
    final$risk_group<-as.factor(as.character(final$risk_group))
    final$type<-factor(final$type,levels=c("mutLoad_nonsilent","mutLoad_silent","CNA_frac_altered","HRD_TAI","HRD_LST","HRD_LOH","HRD_Score","SNV_Neoantigen","mRNAsi"))
    final$value<-as.numeric(final$value)
    p1<-ggboxplot(final,x='risk_group',y='value',fill='risk_group')+xlab("")+ylab("")+
        stat_compare_means(comparisons = my_comparisons,label = "p.format",method = "wilcox.test")+
        scale_fill_manual(values=c( "#377EB8","#FF7F00"))+facet_wrap(~type,nrow=2,scales="free")


    final$risk_score<-as.numeric(final$risk_score)
    p2<-ggplot(data = final, mapping = aes(x = value, y = risk_score)) + geom_point(alpha=0.8,colour="#FC8D62")+
        geom_smooth(method="lm",se=FALSE)+theme_bw()+
        stat_cor(data=final,method="pearson")+
        theme(text=element_text(size=15))+facet_wrap(~type,nrow=2,scales="free")
    pdf(paste0(fileName1,".pdf"), width = 16, height = 10)
    print(p2)
    dev.off()
    tiff(paste0(fileName1,"-72ppi.tif"), width = 16, height = 10, units = "in", res = 72)
    print(p2)
    dev.off()
    tiff(paste0(fileName1,"-300ppi.tif"), width = 16, height = 10, units = "in", res = 300)
    print(p2)
    dev.off()

    write.table(gg,file=paste0(fileName2,"_publicinfor.txt"),quote=FALSE,sep="\t",col.names=NA,row.names=TRUE)

}






#' 高低风险分组免疫评分、基质评分、肿瘤纯度差异
#'
#' @param exp 表达谱
#' @param fileName1 文本文件
#' @param fileName2 图片名称
#' @param dat sample statu os Risk group
#'
#' @return 文本和图片
#' @export
#'
#' @examples myestimate(exp,dat,fileName1,fileName2)
myestimate<-function(exp,dat,fileName1,fileName2){
    write.table(exp,file="purityInput.txt",quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE)
    library(estimate)
    OvarianCancerExpr <- "purityInput.txt"
    filterCommonGenes(input.f=OvarianCancerExpr,output.f="genes.gct",id="GeneSymbol")
    estimateScore(input.ds = "genes.gct",output.ds=paste0(fileName1,"_estimate.txt"),platform="affymetrix")
    result <- read.table(paste0(fileName1,"_estimate.txt"),skip = 2,header = T,sep = "\t")
    colnames(result)<-gsub("[.]","-",colnames(result))
    rownames(result)=result[,1]
    result[1:3,1:3]
    scores=t(result[,3:ncol(result)])
    write.table(scores,paste0(fileName1,"_estimate.txt"),col.names = NA,row.names = T,sep = "\t",quote = F)

    scores<-read.table(paste0(fileName1,"_estimate.txt"),header=TRUE,sep="\t",as.is=TRUE)
    title<-scores[,1]
    scores<-scores[,-1]
    rownames(scores)<-title
    scores<-data.frame(sample=rownames(scores),scores)
    gg<-merge(scores,dat)
    my_comparisons<-list(c("High","Low"))
    p1<-ggviolin(gg,x='Risk',y='StromalScore',fill='Risk',palette=c("#FF7F00","#377EB8"),add=c("boxplot","jitter"),
                 add.params=list(fill="white"),legend="",ylab="Stromal Score",xlab="")+theme(text=element_text(size=18))+
        stat_compare_means(comparisons = my_comparisons, label = "p.format",bracket.size=1,method="wilcox.test")
    p2<-ggviolin(gg,x='Risk',y='ImmuneScore',fill='Risk',palette=c("#FF7F00","#377EB8"),add=c("boxplot","jitter"),
                 add.params=list(fill="white"),legend="",ylab="Immune Score",xlab="")+theme(text=element_text(size=18))+
        stat_compare_means(comparisons = my_comparisons, label = "p.format",bracket.size=1,method="wilcox.test")
    p3<-ggviolin(gg,x='Risk',y='TumorPurity',fill='Risk',palette=c("#FF7F00","#377EB8"),add=c("boxplot","jitter"),
                 add.params=list(fill="white"),legend="",ylab="Tumor Purity",xlab="")+theme(text=element_text(size=18))+
        stat_compare_means(comparisons = my_comparisons, label = "p.format",bracket.size=1,method="wilcox.test")

    pdf(paste0(fileName2,"A.pdf"),width=6,height=8)
    print(p1)
    dev.off()
    tiff(paste0(fileName2,"A-72ppi.tif"),width=6,height=8,units="in",res=72)
    print(p1)
    dev.off()
    tiff(paste0(fileName2,"A-300ppi.tif"),width=6,height=8,units="in",res=300)
    print(p1)
    dev.off()

    pdf(paste0(fileName2,"B.pdf"),width=6,height=8)
    print(p2)
    dev.off()
    tiff(paste0(fileName2,"B-72ppi.tif"),width=6,height=8,units="in",res=72)
    print(p2)
    dev.off()
    tiff(paste0(fileName2,"B-300ppi.tif"),width=6,height=8,units="in",res=300)
    print(p2)
    dev.off()

    pdf(paste0(fileName2,"C.pdf"),width=6,height=8)
    print(p3)
    dev.off()
    tiff(paste0(fileName2,"C-72ppi.tif"),width=6,height=8,units="in",res=72)
    print(p3)
    dev.off()
    tiff(paste0(fileName2,"C-300ppi.tif"),width=6,height=8,units="in",res=300)
    print(p3)
    dev.off()
}

#' GSVA分析及组间比较和绘制热图
#'
#' @param gmtFile gmt文件
#' @param exp 表达谱
#' @param gg 分组信息，第一列是样本名，第二列是cluster1，cluster2
#' @param fileName1 文本名称
#' @param fileName2 图片名称
#'
#' @return 表格和图片
#' @export
#'
#' @examples mygsva(gmtFile,exp,gg,fileName1,fileName2)
mygsva<-function(gmtFile,exp,gg,fileName1,fileName2){
    library("GSEABase")
    library(GSVA)
    library(pheatmap)
    geneSets <- getGmt(gmtFile)
    mydata = as.matrix(exp)
    res_es_TCGA_cancer <- gsva(mydata, geneSets, mx.diff=FALSE, verbose=FALSE,method='gsva')
    write.table(res_es_TCGA_cancer,file=past0(fileName1,"_GSVA_KEGG.txt"),quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE)


    res_es_TCGA_cancer<-res_es_TCGA_cancer[,match(gg[,1],colnames(res_es_TCGA_cancer))]
    wilcox_result <- apply(res_es_TCGA_cancer, 1, function(geneVec){
        result12 <- wilcox.test(geneVec[which(gg[,2] == "cluster1")], geneVec[which(gg[,2] == "cluster2")], paired = FALSE)
        #result13 <- wilcox.test(geneVec[which(sample_class == 1)], geneVec[which(sample_class == 3)], paired = FALSE)
        #result23 <- wilcox.test(geneVec[which(sample_class == 2)], geneVec[which(sample_class == 3)], paired = FALSE)
        return(c(result12$p.value))
    })

    #p值矫正
    wilcox_result_adjust <- p.adjust(wilcox_result, method = "bonferroni")
    pathway<-names(wilcox_result_adjust)[which(wilcox_result_adjust<0.05)]
    length(pathway)
    pathway_top20<-names(sort(wilcox_result_adjust))[1:20]
    write.table(cbind(wilcox_result=wilcox_result,wilcox_result_adjust=wilcox_result_adjust),
                file=paste0(fileName1,"_gsva_kegg.txt"),quote=FALSE,sep="\t",col.names=TRUE,row.names=TRUE)

    ###########画图
    #ann_colors = list(
    #    cluster =c("m6AclusterA"="#FF7F00","m6AclusterB"="#377EB8")
    #    #project=c("TCGA_LGG"="#3B7FB6","CGGA.mRNAseq_325"="#EC1122","CGGA.mRNAseq_693"="#855124")
    #)
    pdf(paste0(fileName2,".pdf"),height=8,width=12)
    pheatmap(res_es_TCGA_cancer[pathway_top20,sample_order],clustering_method = "ward.D",show_rownames =T,show_colnames =F,fontsize_row=10,
             annotation_col = annotation_col,cluster_cols = F,color = colors,scale="row")
    dev.off()
    tiff(paste0(fileName2,"-72ppi.tif"),height=6,width=10,units="in",res=72)
    pheatmap(res_es_TCGA_cancer[pathway_top20,sample_order],clustering_method = "ward.D",show_rownames =T,show_colnames =F,fontsize_row=10,
             annotation_col = annotation_col,cluster_cols = F,color = colors,scale="row")
    dev.off()
    tiff(paste0(fileName2,"-300ppi.tif"),height=6,width=10,units="in",res=300)
    pheatmap(res_es_TCGA_cancer[pathway_top20,sample_order],clustering_method = "ward.D",show_rownames =T,show_colnames =F,fontsize_row=10,
             annotation_col = annotation_col,cluster_cols = F,color = colors,scale="row")
    dev.off()

}



