library(readr)







### integration is merged (shared part of all samples )not all



all_seq_data_integration <- function(TE_add = "Counted/TE" , 
                                 gene_add = "Counted/Gene_count", 
                                 species ='mm39',
                                 gene_count='FPKM', 
                                 TE_count ='tot_counts',
                                 col_data ='Call/GSE211061_coldata.txt',
                                 genomes_for_mapping = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
                                                         'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
                                                         'chrX','chrY')
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
                       by ='TE_ID')
      
      colnames(TE_data)[length(colnames(TE_data))] <-  strsplit(file,"_")[[1]][1]}}
  
  #TE_data <- cbind( TE_data[,c(1,2,3,4)],log2( TE_data[-c(1,2,3,4)]+0.01))
  colnames(TE_data)[c(2,3,4)]<- c('chr','start','end')
  TE_data <- TE_data[which(TE_data$chr %in% genomes_for_mapping),]
  colnames(TE_data)[-c(1:4)] <- factor(colnames(TE_data)[-c(1:4)],level=keys)
  
  
  ############## gene_data_arrange_merge
  
  gene_data<-c() 
  for (file in list.files(gene_add)) {
    if(length(gene_data)==0){
      gene_data <- read_tsv(paste0(gene_add,'/',file),col_names = F)[,c(1,2,3,4,5)] 
      colnames(gene_data)[length(colnames(gene_data))] <- strsplit(file,"_")[[1]][1]
    }else{
      
      gene_data <- merge(gene_data,read_tsv(paste0(gene_add,'/',file))[,c(1,2,3,4,5)],
                         by =c(1,2,3,4))
      
      colnames(gene_data)[length(colnames(gene_data))] <- strsplit(file,"_")[[1]][1]}}
  
  #gene_data <- cbind(gene_data[,c(1,2,3,4)],log2(gene_data[-c(1,2,3,4)]+0.01))
  colnames(gene_data)[c(1,2,3,4)]<- c('chr','start','end','Gene_name')
  gene_data <- gene_data[which(gene_data$chr %in% genomes_for_mapping),]
  
  colnames(gene_data)[-c(1:4)] <- factor(colnames(gene_data)[-c(1:4)],level=keys) 
  #gene_data$start <- gene_data$start - gene_up_stream
  #gene_data$end <- gene_data$end + gene_down_stream
  
  
  return(list(gene_data,TE_data))
}

segement_plot <- function(TE_data='TE_seq_data_integrated.csv',
                          Gene_data='Gene_seq_data_integrated.csv',
                          seg_length = 10000,
                          genomes_for_mapping = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
                                                  'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
                                                  'chrX','chrY'),
                          genome_info ='ref/mm39_chromInfo.txt'){
  
  genome_info <- read_tsv(genome_info,col_names = F)
  genome_info <- genome_info[mapply(function(x) which(x  %in% genomes_for_mapping),genome_info[,1] ),]
  TE_data <- read_csv(TE_data)
  Gene_data <- read_csv(Gene_data)
  
  for (chr in genome_info[,1]){
    position <- 0
    expression <- 0
    while (position <= chr) {
      dataset <- data[which(data$chr == chr),]
      for (var in length(dataset[1,])) {
        
      }
      
      expressed_at_frag <- 
      
      
      
      
    }
    
    
    
    
    
    
  }
  
  
  
  
  
}







c<-all_seq_data_integration ()

write_csv(c[[1]],'Gene_seq_data_integrated.csv')
write_csv(c[[2]],'TE_seq_data_integrated.csv')

segement_plot()







