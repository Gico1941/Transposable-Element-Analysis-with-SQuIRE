library(readr)
#library(circlize)
library(dplyr)
#library(AneuploidyScore) downlaod cytoband of mm39
#ucsc.mm39.cytoband
#library(ComplexHeatmap)



genome_position_file <- function(addr ='Call',
                                 refseq_gene_addr = 'DESeq2_RefSeq_only.txt',
                                 TE_gene_addr = 'DESeq2_TE_only.txt',
                                 p_hold = 10**-4, 
                                 FC_hold = 2,
                                 location_file='gene_location.csv') {
  gene <- read.table(paste0(addr,'/',refseq_gene_addr),row.names = 1)
  loca <- read.csv(location_file)
  rownames(loca)<-paste0(loca$gene,',',loca$strand)
  gene <- gene[which(gene$padj<p_hold&abs(gene$log2FoldChange)>=FC_hold),]
  gene <- merge(loca,gene,by='row.names')
  
  
  TE <- read.table(paste0(addr,'/',TE_gene_addr), row.names = 1)
  #FILTER_OUT
  TE <- TE[which(TE$padj<p_hold&abs(TE$log2FoldChange)>=FC_hold),]
  TE[,c('chr',"start","end","TE","length",'strand')] <- t(data.frame(strsplit(rownames(TE),'\\|')))
  TE$strand <- unlist(lapply(TE$strand,function(x) gsub(",.",'',x)))
  
  
  write_tsv (gene[which(gene$log2FoldChange>0),c('gene','chr','start','end','strand',"log2FoldChange")]
            ,"co-regulation_map/DE_location/up_DE_gene_co-regulation_map_locations.txt",col_names = F)
  write_tsv (gene[which(gene$log2FoldChange<0),c('gene','chr','start','end','strand',"log2FoldChange")]
             ,"co-regulation_map/DE_location/down_DE_gene_co-regulation_map_locations.txt",col_names = F)
  write_tsv (gene[,c('gene','chr','start','end','strand',"log2FoldChange")]
             ,"co-regulation_map/DE_location/all_DE_gene_co-regulation_map_locations.txt",col_names = F)
  
  write_tsv(TE[which(TE$log2FoldChange>0),c('TE','chr','start','end','strand',"log2FoldChange")]
            ,"co-regulation_map/DE_location/up_DE_TE_co-regulation_map_locations.txt",col_names = F)
  write_tsv(TE[which(TE$log2FoldChange<0),c('TE','chr','start','end','strand',"log2FoldChange")]
            ,"co-regulation_map/DE_location/down_DE_TE_co-regulation_map_locations.txt",col_names = F)
  write_tsv(TE[,c('TE','chr','start','end','strand',"log2FoldChange")]
            ,"co-regulation_map/DE_location/all_DE_TE_co-regulation_map_locations.txt",col_names = F)
  
  
  
  
  
  write_tsv (gene[which(gene$log2FoldChange>0),c('gene','chr','start','end','strand',"baseMean")]
             ,"co-regulation_map/DE_basemean/up_DE_gene_co-regulation_map_locations.txt",col_names = F)
  write_tsv (gene[which(gene$log2FoldChange<0),c('gene','chr','start','end','strand',"baseMean")]
             ,"co-regulation_map/DE_basemean/down_DE_gene_co-regulation_map_locations.txt",col_names = F)
  write_tsv (gene[,c('gene','chr','start','end','strand',"baseMean")]
             ,"co-regulation_map/DE_basemean/all_DE_gene_co-regulation_map_locations.txt",col_names = F)
  
  write_tsv(TE[which(TE$log2FoldChange>0),c('TE','chr','start','end','strand',"baseMean")]
            ,"co-regulation_map/DE_basemean/up_DE_TE_co-regulation_map_locations.txt",col_names = F)
  write_tsv(TE[which(TE$log2FoldChange<0),c('TE','chr','start','end','strand',"baseMean")]
            ,"co-regulation_map/DE_basemean/down_DE_TE_co-regulation_map_locations.txt",col_names = F)
  write_tsv(TE[,c('TE','chr','start','end','strand',"baseMean")]
            ,"co-regulation_map/DE_basemean/all_DE_TE_co-regulation_map_locations.txt",col_names = F)
 
}



genome_position_file()

