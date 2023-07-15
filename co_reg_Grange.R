library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(EnrichedHeatmap)

TE_position <- function(DE_addr='co-regulation_map/DE_location',
                        TE_expansion = 'up_DE_TE_Homer_locations.txt',
                        Gene_expansion = 'up_DE_gene_Homer_locations.txt'){
  
  
  
  
  Grange_generator <- function(addr,expansion,type){
    
    file <- read_tsv(paste0(addr,'/',expansion),col_names = F) 
    file <- file[order(file$X6,decreasing = T),]
    file$range <- lapply(c(1:length(row.names(file))),
                       function(x) c(as.numeric(file[x,3]):as.numeric(file[x,4])))
    
    if(type=='Gene'){
      file <- GRanges(seqnames=file$X2,
                      ranges=unlist(file$range),
                      strand=file$X5,
                      FC = file$X6)
    }else{
      file$range <- lapply(c(1:length(row.names(file))),
                           function(x) (as.numeric(file[x,3])+as.numeric(file[x,4]))/2)
      file <- GRanges(seqnames=file$X2,
                      ranges=unlist(file$range),
                      strand=file$X5)
    }
    
    return(file)
  }
  
  
  ########################
  Gene <- Grange_generator(DE_addr,Gene_expansion,'Gene')
  
  TE <- Grange_generator(DE_addr,TE_expansion,'TE')
  
  mat <- normalizeToMatrix(Gene, TE, value_column = "FC", 
                    extend = 3000, mean_mode =  "coverage", w = 50)
  
  
  
  EnrichedHeatmap(mat,col = c("blue","red" ), name = "log2FC")
  
  
}
TE_position ()




###################################






















######################################

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)

TE_position_CS <- function(DE_addr='co-regulation_map/DE_location',
                        TE_expansion = 'up_DE_TE_co-regulation_map_locations.txt',
                        Gene_expansion = 'up_DE_gene_co-regulation_map_locations.txt'){
  
  
  
  
  Grange_generator <- function(addr,expansion){
    
    file <- read_tsv(paste0(addr,'/',expansion),col_names = F) 
    file <- file[order(file$X6,decreasing = T),]
    file$range <- lapply(c(1:length(row.names(file))),
                         function(x) c(as.numeric(file[x,3]):as.numeric(file[x,4])))
      file <- GRanges(seqnames=file$X2,
                      ranges=unlist(file$range),
                      strand=file$X5,
                      FC = file$X6)
    return(file)
  }
  
  
  ########################
  Gene <- Grange_generator(DE_addr,Gene_expansion)
  
  TE <- Grange_generator(DE_addr,TE_expansion)
  
}


Mat<-getTagMatrix(peak = Gene, TxDb = txdb, 
             upstream = 3000, downstream = 3000, 
             type = "start_site", by = "gene", 
             weightCol = "FC")

plotPeakProf(Mat)


Mat<-getTagMatrix(peak = TE, TxDb = txdb, 
                  upstream = 3000, downstream = 5000, 
                  type = "end_site", by = "gene", 
                  weightCol = "FC")


plotPeakProf(Mat)





