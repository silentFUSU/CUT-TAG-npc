set.seed(1)  
.libPaths(c("/storage/zhangyanxiaoLab/suzhuojie/R/x86_64-redhat-linux-gnu-library/4.2",
            "/storage/zhangyanxiaoLab/suzhuojie/miniconda3/envs/R4.2/lib/R/library"))
setwd("/storage/zhangyanxiaoLab/suzhuojie/projects/npc_yulab/")
library(magrittr)
library(BSgenome)
library(GenomicFeatures)
library(GenomicRanges)
library(Rsubread)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
#library(DiffBind)
library(edgeR)
library(ChIPseeker)
library(clusterProfiler)
library(EnsDb.Hsapiens.v86)
library(tidyr)
edb <- EnsDb.Hsapiens.v86
antibodys <-c("1","2","4","5","6","8","13","14","15","16","18")
seqlevelsStyle(edb) <- "UCSC"
GO_database <- 'org.Hs.eg.db'
timepoint <- c("D0","D1","D2","D3","D6")
genes <- c("PAX6","POU5F1","NANOG","SOX1","SOX2","CBX2")
for(i in 1:11){
  antibody <- antibodys[i]
  tab = read.delim(paste0("data/Cut_tag/D0-D6/",antibody,"/",antibody,"_merged_peaks.counts"),skip=1)
  counts = tab[,c(7:16)]
  rownames(counts)= tab$Geneid
  colnames(counts) = paste(rep(c("D0","D1","D2","D3","D6"),each=2), rep(c(1,2),2),sep='_')
  group =c("D0","D0","D1","D1","D2","D2","D3","D3","D6","D6")
  design = model.matrix(~0+group)
  lvls = levels(factor(group))
  colnames(design) = lvls
  len = length(lvls)
  myContrasts = paste(lvls[2:len],lvls[1:(len-1)],sep='-')
  contrast.matrix = eval(as.call(c(as.symbol("makeContrasts"),as.list(myContrasts),levels=list(design))))
  y= DGEList(counts=counts,group=group)
  y =  calcNormFactors(y)
  y<-estimateCommonDisp(y)
  y<-estimateGLMTagwiseDisp(y,design)
  fit_tag = glmFit(y,design)
  lrt = glmLRT(fit_tag, contrast = contrast.matrix)
  qBH = p.adjust(lrt$table$PValue,method="BH")
  out = cbind(tab[,1:6],cpm(y),lrt$table,qBH)
  out[,c(7:24)] = sweep(out[,c(7:24)],1,out$Length,'/')*1000
  for (l in 1:length(myContrasts)){
    lrt = glmLRT(fit_tag, contrast = contrast.matrix[,l])
    out[,paste0("PValue.",myContrasts[l])] = lrt$table$PValue
  }
  mds <- as.data.frame(t(out[,c(7:24)]))
  mds <- mds[1:10,]
  mds_col <- mds %>%
    dist() %>%          
    cmdscale() %>%
    as_tibble()
  colnames(mds_col) <- c("Dim.1", "Dim.2")
  ggplot(data = mds_col, mapping = aes(x = Dim.1, y = Dim.2, colour = group)) + 
    geom_point(size = 3)+  
    geom_text_repel(
      aes(label = rownames(mds)),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"))+
    # 图例
    theme(legend.position = "bottom",panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme_bw()
  ggsave(paste0("result/cut_tag/D0-D6/to_SONG/",antibody,"/",antibody,"_mds.png"),width = 6,height = 5)
  peak <- readPeakFile(paste0("data/Cut_tag/D0-D6/",antibody,"/",antibody,"_merged_peaks.bed"))
  peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                           TxDb=edb, annoDb="org.Hs.eg.db")
  peakAnno<-as.data.frame(peakAnno)
  out$annotation <- peakAnno$annotation
  out$ens_id <- peakAnno$geneId
  out$symbol <- peakAnno$SYMBOL
  
  keep = which(rowSums(cpm(y)>1)>=2)
  out<-out[keep,]
  #write.table(out,paste0("result/cut_tag/D0-D6/",antibody,"/",antibody,".txt"),row.names=F,sep='\t',quote=F)
  out_gene<-out[-which(is.na(out$symbol)),]
  out_gene<-out_gene[which(str_detect(out_gene$annotation,"Promoter")),]
  for(j in 2:5){
    genelist_up <-bitr(out_gene$symbol[which(out_gene[,23+j] < 0.05 & out_gene[,15+j]>= 1)],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database) 
    genelist_up_GO <-enrichGO( genelist_up$ENTREZID,#GO富集分析
                               OrgDb = GO_database,
                               keyType = "ENTREZID",#设定读取的gene ID类型
                               ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                               pvalueCutoff = 0.05,#设定p值阈值
                               qvalueCutoff = 0.05,#设定q值阈值
                               readable = T)
    if(nrow(genelist_up_GO) > 0){ 
      dotplot(genelist_up_GO,font.size=15)
      ggsave(paste0("result/cut_tag/D0-D6/to_SONG/",antibody,"/",antibody,"_",timepoint[j],"_",timepoint[j-1],"_up.png"),width=9,height=7)
      gene_list <- genelist_up_GO@result
      gene_list <- genelist_up_GO[1:10,]
      gene_list<- separate(gene_list,geneID,paste0("gene",1:max(genelist_up_GO$Count)),"/")
      write.csv(gene_list,paste0("result/cut_tag/D0-D6/to_SONG/",antibody,"/",antibody,"_",timepoint[j],"_",timepoint[j-1],"_up_genelist.csv"))
      }
    genelist_down <-bitr(out_gene$symbol[which(out_gene[,23+j] < 0.05 & out_gene[,15+j] <= -1)],fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database) 
    genelist_down_GO <-enrichGO( genelist_down$ENTREZID,#GO富集分析
                               OrgDb = GO_database,
                               keyType = "ENTREZID",#设定读取的gene ID类型
                               ont = "BP",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
                               pvalueCutoff = 0.05,#设定p值阈值
                               qvalueCutoff = 0.05,#设定q值阈值
                               readable = T)
    if(nrow(genelist_down_GO) > 0){ 
      dotplot(genelist_down_GO,font.size=15)
      ggsave(paste0("result/cut_tag/D0-D6/to_SONG/",antibody,"/",antibody,"_",timepoint[j],"_",timepoint[j-1],"_down.png"),width=9,height=7)
      gene_list <- genelist_down_GO@result
      gene_list <- genelist_down_GO[1:10,]
      gene_list<- separate(gene_list,geneID,paste0("gene",1:max(genelist_down_GO$Count)),"/")
      write.csv(gene_list,paste0("result/cut_tag/D0-D6/to_SONG/",antibody,"/",antibody,"_",timepoint[j],"_",timepoint[j-1],"_down_genelist.csv"))
    }
  }
  for(k in 1:6){
    gene<-genes[k]
    out_gene_interest<-out_gene[which(str_detect(out_gene$symbol,gene)),]
    out_gene_interest<-out_gene_interest[,c(1:6,17:20,25:31)]
    write.csv(out_gene_interest,paste0("result/cut_tag/D0-D6/to_SONG/",antibody,"/",gene,".csv"))
  }
}

# out_gene_PAX6 <- out_gene[which(out_gene$`PValue.D3-D2`<0.01 & out_gene$logFC.D3.D2>1 & out_gene$Length >1000),]
# out_gene_PAX6 <- out_gene_PAX6[,c(1:6,19,27,31)]
