library(readr)
library(dplyr)
library(reshape2)           unfinished script
library(ChIPseeker)  


co_reg_map <- function(DE_addr='co-regulation_map/DE_location',
                       flank = 10000,
                       frag_size=501 #(seperate flank region to n frags)
                       ){
  
  up_gene <- read_tsv(paste0(DE_addr,'/','up_DE_gene_Homer_locations.txt'),col_names = F) 
  
  up_TE <- read_tsv(paste0(DE_addr,'/','up_DE_TE_Homer_locations.txt'),col_names = F) 
  
  #down_gene <- read_tsv(paste0(DE_addr,'/','down_DE_gene_Homer_locations.txt'),col_names = F) 

  #down_TE <- read_tsv(paste0(DE_addr,'/','down_DE_TE_Homer_locations.txt'),col_names = F) 
  
  
  
  ###########
  #sort TE by log2(Fc) FC>2,p<10 -4
  
  
  
  #initialize genome map 
  position_array <- function(TE = up_TE, Gene = up_gene,flank=flank){
    
  TE <- TE[order(up_TE$X6,decreasing = T),]
  TE$pos <- (TE$X3+TE$X4)/2

    

  ####
  pos_dt <-c()
  for (n in c(1:length(row.names(TE)))){
    candidate_genes <- left_join(TE[n,],Gene, by = c("X2","X5"),multiple = "all")
    
    #filter by flanks _ distance
    candidate_genes$X3.y <- candidate_genes$X3.y-candidate_genes$pos + flank
    candidate_genes$X4.y <- candidate_genes$X4.y-candidate_genes$pos + flank
    candidate_genes <- candidate_genes[which( (0 < candidate_genes$X3.y&
                                             candidate_genes$X3.y < flank*2+1) | 
                                             (0 < candidate_genes$X4.y &
                                             candidate_genes$X4.y< flank*2+1) |
                                              (candidate_genes$X3.y <0 
                                               & candidate_genes$X4.y >flank*2+1)),]
    if(length(pos_dt)==0){
      pos_dt <-candidate_genes
    }else{
      pos_dt <- rbind(pos_dt,candidate_genes)}}
  
  #### position convert (flank/100 as 1 unit) inclusive mode
  
  pos_dt$X3.y <- unlist(lapply(pos_dt$X3.y, function(x) max(x%/%(flank*2/frag_size),0 )) )
  pos_dt$X4.y <- unlist(lapply(pos_dt$X4.y, function(x) min(x%/%(flank*2/frag_size),frag_size )) ) 
  
  
  #re arrange data for heatmap
  pos_dt$region <- lapply(c(1:length(row.names(pos_dt))),
                          function(x) c(as.numeric(pos_dt[x,9]):as.numeric(pos_dt[x,10])))
  
  # initialize pos_matrix
  
  pos_matrix <- matrix(0,ncol = frag_size,nrow = length(row.names(TE)))
  
  lapply(unique(pos_dt$X1.x), function(x) matrix)
    
  }}


library(EnrichedHeatmap)
set.seed(123)
load(system.file("extdata", "chr21_test_data.RData", package = "EnrichedHeatmap"))
ls()
tss = promoters(genes, upstream = 0, downstream = 1)
mat1 = normalizeToMatrix(H3K4me3, tss, value_column = "coverage", 
                         extend = 5000, mean_mode = "w0", w = 50)


DataFrame(score=11:16, GC=seq(1, 0, length=6))



