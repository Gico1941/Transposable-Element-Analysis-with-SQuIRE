#library(dplyr)
rm(list = ls())

library(readr)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library('dplyr')

#library(pathview)
#library(gage)
#library(gageData)
#library("AnnotationDbi")
#library("org.Mm.eg.db")

#library(clusterProfiler)
#library(GSEABase)
#data(kegg.sets.mm)
#data(sigmet.idx.mm)

#columns(org.Mm.eg.db)

#kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]
#library(tidyverse)

#names <- rownames(gene)

#data[grepl('Art2b',rownames(data)),]



data_disassemble <- function(subs = 'Counted/all/GSE211061_gene_subF_counttable.txt',
                             locus = 'Counted/all/GSE211061_gene_TE_counttable.txt'){

  read <- function(addr){
    data <- read.table(addr, sep = "\t", row.names = 1)
    data <-read.table(addr, sep = "\t", row.names = 1)
    colnames(data)<-data[1,]
    data<-data[-1,]
    return(data)
  }
  
    subs_ <- read(subs)
    locus_ <- read(locus)
  
  gene <- subs_[-grep(':',row.names(subs_)),]
  
  sub_TE <- subs_[grep(':',row.names(subs_)),]
  
  locus_TE <- locus_[grep(':',row.names(locus_)),]
  
  
  
  #########
  
  symbol <- lapply(row.names(gene),function(x) strsplit(x,',')[[1]][1])
  
  gene <- aggregate(x=mapply(as.numeric,gene),by=list(unlist(symbol)),FUN = sum)
  colnames(gene)[1]<-'symbol'
  
  locus_TE$symbol <- row.names(locus_TE)
  #locus_TE$symbol <- unlist( lapply(row.names(locus_TE),function(x) strsplit(x,'\\|')[[1]][4]))
  locus_TE <- locus_TE %>% select(symbol,everything())

  
  sub_TE$symbol <- row.names(sub_TE)
  sub_TE <- sub_TE %>% select(symbol,everything())
  
  ####
  
  
  write.table(gene,'Counted/gene_count.txt',sep='\t',quote=F,row.names = F)
  
  
  write.table(sub_TE,'Counted/sub_TE_count.txt',sep='\t',quote=F,row.names = F)
  

  write.table(locus_TE,'Counted/locus_TE_count.txt',sep='\t',quote=F,row.names = F)
  
}



data_disassemble(subs = 'Counted/all/GSE211061_gene_subF_counttable.txt',
                 locus = 'Counted/all/GSE211061_gene_TE_counttable.txt')










PCA_plot <- function(addr='Counted/all/GSE211061_gene_TE_counttable.txt'
                     ,decode='Call/GSE211061_coldata.txt',
                     plot='gene'){
  
  
  samples <- read.csv(decode, sep = "\t", row.names = 1)
  
  data <-read.table(addr, sep = "\t", row.names = 1)
  colnames(data)<-data[1,]
  data<-data[-1,]
  


  TE <- mapply(as.numeric,data[grep(':',row.names(data)),])
  gene <-  mapply(as.numeric,data[-grep(':',row.names(data)),])
  all <- mapply(as.numeric,data)

   
   
   
   
   
  min_read <- 0
  
  PCA_sub_plot <- function(data=gene){
    data <- data[apply(data,1,function(x){max(x)}) >= min_read,]
    groups <- factor(c(rep("Female_LPS",3),rep("Male_LPS",3),rep("Female_CTRL",3),rep("Male_CTRL",3)))
    sampleInfo <- data.frame(groups,row.names=colnames(data))
    dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)
    vsd <- vst(dds, blind = FALSE)
    plotPCA(vsd,intgroup=("groups"))
  }
  
  if(plot=='gene'){
    c<- PCA_sub_plot(gene) #gene,TE,all
  }
  if(plot=='TE'){
    c<- PCA_sub_plot(TE) #gene,TE,all
  }
  if(plot=='all'){
    c<- PCA_sub_plot(all) #gene,TE,all
  }

  
  #c<- lapply(list(data,gene,TE), PCA_sub_plot )
 return (c)
}

PCA_plot('Counted/all/GSE211061_gene_TE_counttable.txt',plot='gene')


PCA_plot('Counted/all/GSE211061_gene_subF_counttable.txt',plot='TE')












pathway_enrich <- function(addr ='Call',data_addr = 'DESeq2_RefSeq_only.txt',p_hold = 10**-6, FC_hold = 2,enrichpaths = 15){
  
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
                      organism = 'mmu',
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
  KEGG_plot(greaters,'greaters')
  KEGG_plot(lowers,'lowers')
}


volcano_plot_p_FC <- function(addr ='Call',data_addr = 'DESeq2_RefSeq_only.txt',p_hold = 10**-4, FC_hold = 2,
                         MAX=5,top_display=15){
  
  ##vlocanno scatter plot and add significance
  data <- read.table(paste0(addr,'/',data_addr), sep = "\t", row.names = 1)
  #categorize dataset
  data$significant <- ifelse(data$padj<p_hold&abs(data$log2FoldChange)>=FC_hold,
                             ifelse(data$log2FoldChange>FC_hold,"Up","Down"),"Stable")
  
  if (data_addr == 'DESeq2_TE_only.txt'){
    if(addr =='Call/sub_level'){
      data $Symbol <- row.names(data)
    }else{
      data $Symbol <- unlist(lapply(row.names(data), function(x) strsplit(x,"\\|")[[1]][4]))
    }
    
   }else{
     data $Symbol <- unlist(lapply(row.names(data), function(x) strsplit(x,",")[[1]][1]))
    }
  
  data <- data[complete.cases(data),]

  plot <- ggplot(data,aes(x=log2FoldChange,y=-log10(pvalue)))+
    theme(plot.margin = unit(c(1,1,1,1),"cm"))+
    geom_point(aes(color = significant),size=2,alpha=0.4)+
    scale_color_manual(values=c('blue','grey','red')[which(c("Down","Stable","Up") %in% unique(data$significant))])+
    geom_text_repel(
      data = rbind(subset(data,padj<p_hold&data$log2FoldChange>=FC_hold) %>%top_n(top_display,log2FoldChange), 
                   subset(data,padj<p_hold&data$log2FoldChange<=-FC_hold)%>%top_n(top_display,log2FoldChange)),
                   
      aes(label = Symbol),size = 3,
      min.segment.length = 0,box.padding = 0.1,max.overlaps =MAX )+
    
    geom_vline(xintercept = c(-FC_hold,FC_hold),lty=4,col='black',lwd=0.8)+
    geom_hline(yintercept = -log10(p_hold),lty=4,col='black',lwd=0.8)+
    theme(legend.position = 'bottom')
  ggsave(paste0(data_addr,'Volcano_plot.pdf'),plot,dpi=1500,width = 10,height = 10)
  return(plot)
}


c <-PCA_plot()

#gene_DE <- pathway_enrich(addr ='Call',data_addr = 'DESeq2_RefSeq_only.txt',p_hold = 10**-6, FC_hold = 2,enrichpaths = 15)

gene_DE <- volcano_plot_p_FC(addr ='Call',data_addr = 'DESeq2_RefSeq_only.txt',p_hold = 10**-4, FC_hold = 2,MAX=10,top_display=50)

TE_DE <- volcano_plot_p_FC(addr ='Call/locus_level',data_addr = 'DESeq2_TE_only.txt',p_hold = 10**-4, FC_hold = 2)



volcano_plot_p_FC(addr ='Call/sub_level',data_addr = 'DESeq2_TE_only.txt',p_hold = 10**-4, FC_hold = 2)







volcano_plot_baseMean_FC <- function(addr ='Call',data_addr = 'DESeq2_RefSeq_only.txt',p_hold = 10**-4, FC_hold = 2,
                              MAX=5,top_display=15){
  
  ##vlocanno scatter plot and add significance
  data <- read.table(paste0(addr,'/',data_addr), sep = "\t", row.names = 1)
  #categorize dataset
  data$significant <- ifelse(data$padj<p_hold&abs(data$log2FoldChange)>=FC_hold,
                             ifelse(data$log2FoldChange>FC_hold,"Up","Down"),"Stable")
  
  if (data_addr == 'DESeq2_TE_only.txt'){
    if(addr =='Call/sub_level'){
      data $Symbol <- row.names(data)
    }else{
      data $Symbol <- unlist(lapply(row.names(data), function(x) strsplit(x,"\\|")[[1]][4]))
    }
    
  }else{
    data $Symbol <- unlist(lapply(row.names(data), function(x) strsplit(x,",")[[1]][1]))
  }
  
  data <- data[complete.cases(data),]
  
  plot <- ggplot(data,aes(x=log2FoldChange,y=log2(baseMean)))+
    theme(plot.margin = unit(c(1,1,1,1),"cm"))+
    geom_point(aes(color = significant),size=2,alpha=0.4)+
    scale_color_manual(values=c('blue','grey','red')[which(c("Down","Stable","Up") %in% unique(data$significant))])+
    geom_text_repel(
      data = rbind(subset(data,padj<p_hold&data$log2FoldChange>=FC_hold) %>%top_n(top_display,log2FoldChange), 
                   subset(data,padj<p_hold&data$log2FoldChange<=-FC_hold)%>%top_n(top_display,log2FoldChange)),
      
      aes(label = Symbol),size = 3,
      min.segment.length = 0,box.padding = 0.1,max.overlaps =MAX )+
    
    #geom_vline(xintercept = c(-FC_hold,FC_hold),lty=4,col='black',lwd=0.8)+
    #geom_hline(yintercept = -log10(p_hold),lty=4,col='black',lwd=0.8)+
    theme(legend.position = 'bottom')
  ggsave(paste0(data_addr,'Volcano_plot.pdf'),plot,dpi=1500,width = 10,height = 10)
  return(plot)
}





volcano_plot_baseMean_FC(addr ='Call',data_addr = 'DESeq2_RefSeq_only.txt',p_hold = 10**-4, FC_hold = 2,MAX=10,top_display=50)


TE_DE <- volcano_plot_baseMean_FC(addr ='Call/locus_level',data_addr = 'DESeq2_TE_only.txt',p_hold = 10**-4, FC_hold = 2)



volcano_plot_baseMean_FC(addr ='Call/sub_level',data_addr = 'DESeq2_TE_only.txt',p_hold = 10**-4, FC_hold = 2)




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

