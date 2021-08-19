#' 下载样本的临床信息
#'
#' @param cancer ："hnsc_tcga"
#' @param fileName ：输出文件名称
#'
#' @return：result<-list(clin_down=clin_down,survival=tcga_pdata_os)
#' @export
#'
#' @examples getClin(cancer="hnsc_tcga",fileName="./result/Table/SupplementaryTable1_clin.txt")
getClin<-function(cancer,fileName){
	library("cgdsr")
	mycgds = CGDS("http://www.cbioportal.org/")
    a<-getCancerStudies(mycgds)
    mycancerstudy = a[which(a[,1]==cancer),1]
    b<-getCaseLists(mycgds,mycancerstudy)
    mycaselist = b[1,1]
    clin_down<-getClinicalData(mycgds,mycaselist)
    rownames(clin_down)<-gsub("[.]","-",rownames(clin_down))
    rownames(clin_down)<-paste0(rownames(clin_down),"A")
    tcga_pdata_os<-data.frame(sample=rownames(clin_down),status_os=clin_down$OS_STATUS,time_os=clin_down$OS_MONTHS,status_dfs=clin_down$DFS_STATUS,time_dfs=clin_down$DFS_MONTHS)
    tcga_pdata_os[,2]<-gsub(":LIVING","",tcga_pdata_os[,2])
    tcga_pdata_os[,2]<-gsub(":DECEASED","",tcga_pdata_os[,2])
    tcga_pdata_os[,4]<-gsub(":DiseaseFree","",tcga_pdata_os[,4])
    tcga_pdata_os[,4]<-gsub(":Recurred/Progressed","",tcga_pdata_os[,4])
    result<-list(clin_down=clin_down,survival=tcga_pdata_os)
    write.table(tcga_pdata_os,fileName,quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
    return(result)
}

