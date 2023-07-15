library(readr)







### integration is merged of all samples



all_seq_data_integration <- function(TE_add = "Counted/TE" , 
                                 gene_add = "Counted/Gene_count", 
                                 TE_count ='tot_counts',
                                 gene_count=7,
                                 col_data ='Call/GSE211061_coldata.txt',
                                 genomes_for_mapping = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
                                                         'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
                                                         'chrX','chrY','chrM')
                                 ,gene_down_stream=3000,gene_up_stream=3000
){
  
  ####### creat group/condition levels
  
  col <- read_tsv(col_data)
  keys <- c()
  for(i in unique(col$condition)){
    keys <- c(keys,gsub("\\_.*","",col$sample[which(col$condition == i)]))
  }
  
  
  ############# TE_data_arrange_merge
  TE_data<-c()
  for (file in list.files(TE_add)) {
    if(length(TE_data)==0){
      TE_data <- read_tsv(paste0(TE_add,'/',file))[,c('TE_chr','TE_start','TE_stop','TE_ID',TE_count)] 
      colnames(TE_data)[length(colnames(TE_data))] <-  strsplit(file,"_")[[1]][1]
    }else{
      
      TE_data <- merge(TE_data,read_tsv(paste0(TE_add,'/',file))[,c('TE_ID',TE_count)],
                       by.x ='TE_ID',by.y ='TE_ID',all=T)
      
      colnames(TE_data)[length(colnames(TE_data))] <-  strsplit(file,"_")[[1]][1]}}
  
  #TE_data <- cbind( TE_data[,c(1,2,3,4)],log2( TE_data[-c(1,2,3,4)]+0.01))
  colnames(TE_data)[c(2,3,4)]<- c('chr','start','end')
  TE_data <- TE_data[which(TE_data$chr %in% genomes_for_mapping),]
  ###group samples
  colnames(TE_data)[-c(1:4)] <- colnames(TE_data)[-c(1:4)][order(match(colnames(TE_data)[-c(1:4)], keys))]
  
  TE_data[-c(1:4)][is.na(TE_data[-c(1:4)])] <- 0
  ############## gene_data_arrange_merge
  
  gene_data<-c() 
  for (file in list.files(gene_add)) {
    if(length(gene_data)==0){
      gene_data <- read_tsv(paste0(gene_add,'/',file),col_names = F)[,c(1,2,3,4,6,gene_count)] 
      colnames(gene_data)[length(colnames(gene_data))] <- strsplit(file,"_")[[1]][1]
    }else{
      
      gene_data <- merge(gene_data,read_tsv(paste0(gene_add,'/',file))[,c(1,2,3,4,6,gene_count)],
                         by.x =c(1,2,3,4,5),by.y =c(1,2,3,4,5),all=T)
      
      colnames(gene_data)[length(colnames(gene_data))] <- strsplit(file,"_")[[1]][1]}}
  
  #gene_data <- cbind(gene_data[,c(1,2,3,4)],log2(gene_data[-c(1,2,3,4)]+0.01))
  colnames(gene_data)[c(1,2,3,4,5)]<- c('chr','start','end','Gene_name','strand')
  gene_data <- gene_data[which(gene_data$chr %in% genomes_for_mapping),]
  ###group samples
  colnames(gene_data)[-c(1:4)] <- colnames(gene_data)[-c(1:4)][order(match(colnames(gene_data)[-c(1:4)], keys))]
  
  gene_data[-c(1:4)][is.na(gene_data[-c(1:4)])] <- 0
  
  #gene_data$start <- gene_data$start - gene_up_stream
  #gene_data$end <- gene_data$end + gene_down_stream
  
  
  return(list(gene_data,TE_data))
}




c<-all_seq_data_integration ()

write_csv(c[[1]],'Gene_seq_data_integrated.csv')
write_csv(c[[2]],'TE_seq_data_integrated.csv')



