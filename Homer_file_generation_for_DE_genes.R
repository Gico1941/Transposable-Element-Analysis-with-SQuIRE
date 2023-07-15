library(readr)
#library(circlize)
library(dplyr)
#library(AneuploidyScore) downlaod cytoband of mm39
#ucsc.mm39.cytoband
#library(ComplexHeatmap)



genome_position_file <- function(addr ='Call',
                                 refseq_gene_addr = 'DESeq2_RefSeq_only.txt',
                                 p_hold = 10**-4, 
                                 FC_hold = 2,
                                 mode='output_names_only',
                                 threshold=T) {
  
  
  
  
  gene <- read.table(paste0(addr,'/',refseq_gene_addr),row.names = 1)
  
  if(threshold==T){
    gene <- gene[which(gene$padj<p_hold&abs(gene$log2FoldChange)>=FC_hold),]
    
  }

  gene <- gene[order(gene$baseMean,decreasing = T),]
  
  
  gene $ Symbol <- lapply(row.names(gene),function(x) strsplit(x,',')[[1]][1]) %>% unlist()
  gene $ Strand <- lapply(row.names(gene),function(x) strsplit(x,',')[[1]][2]) %>% unlist()
  

  DE_genelist<- gene[,c('Symbol','Strand','baseMean','log2FoldChange')]
    


    
    
    
    
    
   
    
    

  
 # write_tsv (gene[which(gene$log2FoldChange>0),c('chr','start','end','gene',"log2FoldChange",'strand')]
            #,"Homer/DE_location/up_DE_gene_Homer_locations.txt",col_names = F)
  #write_tsv (gene[which(gene$log2FoldChange<0),c('chr','start','end','gene',"log2FoldChange",'strand')]
             #,"Homer/DE_location/down_DE_gene_Homer_locations.txt",col_names = F)
  
  #write_tsv (gene[,c('chr','start','end','gene',"log2FoldChange",'strand')]
             #,"Homer/DE_location/all_DE_gene_Homer_locations.txt",col_names = F)
  

  
  
  write_tsv(DE_genelist,'Homer/DE_location/all_DE_gene_Homer_symbol_list.txt')
  
  
  
  
}

genome_position_file(p_hold = 0.05,FC_hold=0)

genome_position_file(threshold=F)

