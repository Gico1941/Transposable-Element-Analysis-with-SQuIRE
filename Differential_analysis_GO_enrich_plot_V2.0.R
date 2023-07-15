#library(dplyr)
rm(list = ls())


require(DOSE)
library(readr)
library(fgsea)
library(DESeq2)
library(pathview)
library(gage)
library(gageData)
library("AnnotationDbi")
library("org.Mm.eg.db")
library('dplyr')
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(GSEABase)
library(enrichplot)
data(kegg.sets.mm)
data(sigmet.idx.mm)

#columns(org.Mm.eg.db)

kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]
#library(tidyverse)

#names <- rownames(gene)

#data[grepl('Art2b',rownames(data)),]
#seperate TE and gene original expression data from ccounttable.txt
#DATA <- read_tsv('Call/GSE211061_gene_TE_counttable.txt')
#DATA_gene <- DATA [-c(grep('\\|',DATA$gene_id)),]
#DATA_TE <- DATA [grep('\\|',DATA$gene_id),]
#DATA_gene$gene_id <- unlist(lapply(DATA_gene$gene_id,function(x) strsplit(x,',')[[1]][1]))
#DATA_DES <- DATA_gene[,1]
#DATA_DES$Description <- NA
#
#DATA_gene<- cbind(DATA_DES,DATA_gene[,-c(1)])
#colnames(DATA_gene)[1] <- 'Name'
#write_tsv(DATA_gene,'Call/GSE211061_gene_counttable.txt')
#write_tsv(DATA_TE,'Call/GSE211061_TE_counttable.txt')


c <-PCA_plot()

PCA_plot <- function(addr='Call',GSE='GSE211061'){
  samples <- read.csv(paste0(addr,'/',GSE,"_coldata.txt"), sep = "\t", row.names = 1)
  data <- read.csv(paste0(addr,'/',GSE,"_gene_TE_counttable.txt"), sep = "\t", row.names = 1)
  gene <- data[-grep('\\|',rownames(data)),]
  TE <-  data[grep('\\|',rownames(data)),]
  min_read <- 1
  PCA_sub_plot <- function(data){
    data <- data[apply(data,1,function(x){max(x)}) > min_read,]
    groups <- factor(c(rep("Female_LPS",3),rep("Male_LPS",3),rep("Female_CTRL",3),rep("Male_CTRL",3)))
    sampleInfo <- data.frame(groups,row.names=colnames(data))
    dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)
    vsd <- vst(dds, blind = FALSE)
    plotPCA(vsd,intgroup=("groups"))
  }
  
  c<- PCA_sub_plot(TE) #gene,TE,all
 #c<- lapply(list(data,gene,TE), PCA_sub_plot )
 return (c)
}







pathway_enrich <- function(addr ='Call',data_addr = 'DESeq2_RefSeq_only.txt',p_hold = 10**-6, FC_hold = 2,enrichpaths = 15,org='mmu'){
  
  data <- read.table(paste0(addr,'/',data_addr), sep = "\t", row.names = 1)
  #categorize dataset
  data$significant <- ifelse(data$padj<p_hold&abs(data$log2FoldChange)>=FC_hold,
                             ifelse(data$log2FoldChange>FC_hold,"Up","Down"),"Stable")
  
  data $Symbol <- unlist(lapply(row.names(data), function(x) strsplit(x,",")[[1]][1]))
  

  #add data entrezid
  data <- merge(data,bitr(data$Symbol,
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Mm.eg.db)
                ,by.x='Symbol',by.y='SYMBOL') 
  
  data <- data[complete.cases(data),]
  
  DEgenes<- data[which(data$significant !='Stable'),]
  
  

  #plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu", new.signature=FALSE)
  # plot multiple pathways (plots saved to disk and returns a throwaway list object)
  #tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="mmu"))
  
  
  
  
  ###kegg /GO and bubble plot with cluster profiler
  genelist = DEgenes$log2FoldChange
  names(genelist) <- DEgenes$ENTREZID
  genelist <- sort(genelist,decreasing = T)
  
  
  CP_Kegg <- gseKEGG( geneList = genelist,
                      organism = org,
                      nPerm = 1000,
                      minGSSize = 10,
                      pvalueCutoff = 0.05,
                      verbose = F) #%>% filter(row_number()<=10)
  

  ###### GO
  #CP_GO <- enrichGO(gene         = names(genelist),
                     #OrgDb         = org.Mm.eg.db,
                     #ont           = "all",
                    # pvalueCutoff  = 0.05,
                     #qvalueCutoff  = 0.05,readable = T)
  
  #summary(CP_GO)
  
  KE<-DOSE::setReadable(CP_Kegg,OrgDb = 'org.Mm.eg.db',keyType = 'ENTREZID')
  greaters<- KE@result[which(KE@result$enrichmentScore>0),][order(KE@result[which(KE@result$enrichmentScore>0),]$p.adjust)[1:enrichpaths],]
  lowers <- KE@result[which(KE@result$enrichmentScore<0),][order(KE@result[which(KE@result$enrichmentScore<0),]$p.adjust)[1:enrichpaths],]
  
  

  
  KEGG_plot <- function(x,y){ plot <- ggplot(x,aes(x=enrichmentScore,y=Description))+
    theme_classic()+theme(plot.margin = unit(c(1,1,1,1),"cm"))+
      geom_point(mapping=aes(size=setSize,color=p.adjust),alpha = 0.4)+
    scale_color_gradient(low='red',high='blue')
    ggsave(paste0(y,'_KEGG_enrich.pdf'),plot,dpi=1000)
  } 
  KEGG_plot(greaters,paste0(strsplit(data_addr,'.txt')[[1]][1],'greaters'))
  KEGG_plot(lowers,paste0(strsplit(data_addr,'.txt')[[1]][1],'lowers'))
}




volcano_plot <- function(addr ='Call',data_addr = 'DESeq2_RefSeq_only.txt',p_hold = 10**-6, FC_hold = 2){
  
  ##vlocanno scatter plot and add significance
  data <- read.table(paste0(addr,'/',data_addr), sep = "\t", row.names = 1)
  #categorize dataset
  data$significant <- ifelse(data$padj<p_hold&abs(data$log2FoldChange)>=FC_hold,
                             ifelse(data$log2FoldChange>FC_hold,"Up","Down"),"Stable")
  
  if (data_addr == 'DESeq2_TE_only.txt'){
    data $Symbol <- unlist(lapply(row.names(data), function(x) strsplit(x,"|")[[1]][4]))
   }else{
     data $Symbol <- unlist(lapply(row.names(data), function(x) strsplit(x,",")[[1]][1]))
    }
  
  data <- data[complete.cases(data),]

  plot <- ggplot(data,aes(x=log2FoldChange,y=-log10(pvalue)))+
    theme(plot.margin = unit(c(1,1,1,1),"cm"))+
    geom_point(aes(color = significant),size=2,alpha=0.4)+
    scale_color_manual(values=c('blue','grey','red'))+
    geom_text_repel(
      data = subset(data,padj<p_hold&abs(data$log2FoldChange)>=FC_hold),
      aes(label = Symbol),size = 5,
      box.padding = unit(0.35,"lines"),
      point.padding = unit(0.3,"lines"))+
    geom_vline(xintercept = c(-FC_hold,FC_hold),lty=4,col='black',lwd=0.8)+
    geom_hline(yintercept = -log10(p_hold),lty=4,col='black',lwd=0.8)+
    theme(legend.position = 'bottom')
  ggsave(paste0(data_addr,'Volcano_plot.pdf'),plot,dpi=1000)
  return(plot)
}



GO_enrich_and_GSEA_plot <- function(addr ='Call',data_addr = 'DESeq2_RefSeq_only.txt',p_hold = 10**-6, FC_hold = 2,enrichpaths = 15,org='mmu'){
  data <- read.table(paste0(addr,'/',data_addr), sep = "\t", row.names = 1)
  #categorize dataset
  data$significant <- ifelse(data$padj<p_hold&abs(data$log2FoldChange)>=FC_hold,
                             ifelse(data$log2FoldChange>FC_hold,"Up","Down"),"Stable")
  
  data $Symbol <- unlist(lapply(row.names(data), function(x) strsplit(x,",")[[1]][1]))
  
  
  #add data entrezid
  data <- merge(data,bitr(data$Symbol,
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Mm.eg.db)
                ,by.x='Symbol',by.y='SYMBOL') 
  
  data <- data[complete.cases(data),]
  
  DEgenes<- data[which(data$significant !='Stable'),]

  
  ###GO and GSEA plot
  genelist = DEgenes$log2FoldChange
  names(genelist) <- DEgenes$ENTREZID
  genelist <- sort(genelist,decreasing = T)
  
  gse <- gseGO(geneList=genelist, 
                 ont ="ALL", 
                 keyType = "ENTREZID",
                 nPerm = 10000, 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = org.Mm.eg.db)
  
  #write_csv(CP_GO,paste0(strsplit(data_addr,'.txt')[[1]][1],'_GO_enrich.csv'))
  return(gse) 
  
  #pdf("GO_enrich.pdf")
  #dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
  
  #dev.off()
  #pdf("GO_enrich__.pdf")
  #return(gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1))
  #dev.off()
 
}


c <-PCA_plot()

gene_DE <- pathway_enrich(addr ='Call',data_addr = 'DESeq2_RefSeq_only.txt',p_hold = 10**-6, FC_hold = 2,enrichpaths = 15)



gene_DE <- volcano_plot(addr ='Call',data_addr = 'DESeq2_RefSeq_only.txt',p_hold = 10**-6, FC_hold = 2)

TE_DE <- volcano_plot(data_addr = 'DESeq2_TE_only.txt',p_hold = 10**-6, FC_hold = 2)

############GO and plot
c <- GO_enrich_and_GSEA_plot(addr ='Call',data_addr = 'DESeq2_RefSeq_only.txt',p_hold = 10**-6, FC_hold = 2,enrichpaths = 15)
dotplot(c, showCategory=10, split=".sign") + facet_grid(.~.sign)

for (i in c(1:10)){
  pdf(file = paste0('result/GSEA_plot/',c$Description[i],"_GO_enrich_GSEA.pdf"),width = 15,height = 8)
  print(gseaplot2(c, by = "all", title = c$Description[i], geneSetID = i))
  dev.off()
}



print(gseaplot2(c, title = c$Description[i], geneSetID = i,pvalue_table = T))

#####kegg with gage

#data$entrez = mapIds(org.Mm.eg.db,
#keys=data$Symbol , 
#column="ENTREZID",
#keytype="SYMBOL",
#multiVals="first")

#foldchanges = data$log2FoldChange
#names(foldchanges) = data$entrez
#keggres = gage(foldchanges, gsets=kegg.sets.mm, same.dir=FALSE)
#lapply(keggres, head)
#keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater)%>% 
#tibble::as_tibble()%>% 
#filter(row_number()<=20)
#keggrespathways$id <- as.character(keggrespathways$id) 
#keggresids = substr(keggrespathways, start=1, stop=8)

