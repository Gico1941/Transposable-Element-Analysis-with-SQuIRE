library(readr)
#library(circlize)
library(dplyr)
#library(AneuploidyScore) downlaod cytoband of mm39
#ucsc.mm39.cytoband
#library(ComplexHeatmap)



genome_position_file <- function(addr ='Call',
                                 TE_gene_addr = 'DESeq2_TE_only.txt',
                                 level='locus_level',
                                 p_hold = 10**-4, 
                                 FC_hold = 2,
                                 length_restrction=F,
                                 top=F,
                                 percent=F) {

  
  
  
  TE <- read.table(paste0(addr,'/',level,'/',TE_gene_addr), row.names = 1)
  #FILTER_OUT
  TE <- TE[which(TE$padj<p_hold&abs(TE$log2FoldChange)>=FC_hold),]
  TE[,c('chr',"start","end","TE","milliDiv",'strand')] <- t(data.frame(strsplit(rownames(TE),'\\|')))
  
  TE$strand <- unlist(lapply(TE$strand,function(x) strsplit(x,",")[[1]][2]))
  TE <- TE[which(TE$strand !='.'),]
  
  TE_matrix <- mapply(as.numeric,TE[,c('start','end')])
  
  
  
  TE <- TE[order(TE$baseMean,decreasing = TRUE),]
  TE$Fullname <- row.names(TE)
  
  
  if (length_restrction != F){
    TE <- TE[which((TE_matrix[,2] - TE_matrix[,1]) > length_restrction),]
  }
  
  
  if (top != F){
    TE <- TE[top,]
  }
  
  if (percent != F){
    TE $baseMean = TE$baseMean/sum(TE$baseMean)
    p=0
    n=1
    T_<- TE[n,]
    p <- p+TE$baseMean[n]
    
    while(p<percent){
      n <- n + 1
      T_ <- rbind(T_,TE[n,])
      p <- p + TE$baseMean[n]
    }
    TE<-T_
  }
  
  
  
  
  
  
  
  ### baseMean/log2FoldChange
  
  write_tsv(TE[which(TE$log2FoldChange>0),c('chr','start','end','Fullname',"baseMean",'strand')]
            ,"Homer/DE_location/baseMean/up_DE_TE_Homer_locations.txt",col_names = F)
  write_tsv(TE[which(TE$log2FoldChange<0),c('chr','start','end','Fullname',"baseMean",'strand')]
            ,"Homer/DE_location/baseMean/down_DE_TE_Homer_locations.txt",col_names = F)
  write_tsv(TE[,c('chr','start','end','Fullname',"baseMean",'strand')]
            ,"Homer/DE_location/baseMean/all_DE_TE_Homer_locations.txt",col_names = F)
  
  
  write_tsv(TE[which(TE$log2FoldChange>0),c('chr','start','end','Fullname',"log2FoldChange",'strand')]
            ,"Homer/DE_location/log2FC/up_DE_TE_Homer_locations.txt",col_names = F)
  write_tsv(TE[which(TE$log2FoldChange<0),c('chr','start','end','Fullname',"log2FoldChange",'strand')]
            ,"Homer/DE_location/log2FC/down_DE_TE_Homer_locations.txt",col_names = F)
  write_tsv(TE[,c('chr','start','end','Fullname',"log2FoldChange",'strand')]
            ,"Homer/DE_location/log2FC/all_DE_TE_Homer_locations.txt",col_names = F)
  
  
  
  
}




genome_position_file(level = 'locus_level',length_restrction=F,percent=F)     #don't specify top and percent at same time 










############## all expressed TEs

genome_position_file_all <- function(addr ='Counted',
                                 TE_gene_addr = 'locus_TE_count.txt') {
  
  
  
  
  TE <- read_tsv(paste0(addr,'/',TE_gene_addr))
  
  
  
  
  TE$baseMean <- unlist(lapply(1:nrow(TE),function(x) mean(as.numeric(TE[x,c(2:length(colnames(TE)))]))))
  
  #FILTER_OUT

  TE[,c('chr',"start","end","TE","milliDiv",'strand')] <- t(data.frame(strsplit(TE$symbol,'\\|')))
  
  
  
  TE$strand <- unlist(lapply(TE$strand,function(x) strsplit(x,",")[[1]][2]))
  
  TE <- TE[which(TE$strand !='.'),]
  
  TE_matrix <- mapply(as.numeric,TE[,c('start','end')])
  
  
  
  TE <- TE[order(TE$baseMean,decreasing = TRUE),]
  #TE <- TE[which(TE$baseMean>1) ,]
  
  ### baseMean/log2FoldChange
  
  write_tsv(TE[,c('chr','start','end','symbol',"baseMean",'strand')]
            ,"Homer/All_location/all_TE_Homer_locations.txt",col_names = F)

  
  
  
  
}




genome_position_file_all()     #don't specify top and percent at same time 

