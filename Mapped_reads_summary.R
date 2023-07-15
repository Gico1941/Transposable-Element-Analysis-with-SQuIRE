library(readr)




#######   remember to set read_length


all_seq_data_integration <- function(TE_add = "Counted/TE" , 
                                     gene_add = "Counted/all", 
                                     species ='mm39',
                                     gene_count='FPKM', 
                                     TE_count ='tot_reads',
                                     col_data ='Call/GSE211061_coldata.txt'
                                     ,gene_down_stream=3000,gene_up_stream=3000,
                                     reads_length=102
){
  
  ####### creat group/condition levels
  
  col <- read_tsv(col_data)
  keys <- c()
  for(i in unique(col$condition)){
    keys <- c(keys,gsub("\\_.*","",col$sample[which(col$condition == i)]))
  }
  
  
  
  ############# TE_reads_arrange_merge
  TE_reads<-c()
  for (file in list.files(TE_add)) {
    if(length(TE_reads)==0){
      TE_reads <- data.frame(TE_reads = sum(read_tsv(paste0(TE_add,'/',file))[,c(TE_count)] )  )
    }else{
      TE_reads <- cbind(data.frame(TE_reads = sum(read_tsv(paste0(TE_add,'/',file))[,c(TE_count)])),TE_reads)}
    colnames(TE_reads)[length(colnames(TE_reads))] <-  strsplit(file,"_")[[1]][1]}
    
  
  
  
  
  ######################
  
  
  TE_add = "Counted/Gene"
  TE_reads<-c()
  for (file in list.files(TE_add)) {
    if(length(TE_reads)==0){
      data <- read_tsv(paste0(TE_add,'/',file))[,c(5,6,7)] 
      exon_reads<- data.frame(exon_reads = sum(data[,3] * abs(data[,2]-data[,1])/reads_length))
      TE_reads <- data.frame(TE_reads = sum(exon_reads)  )
    }else{
      TE_reads <- cbind(data.frame(TE_reads = sum(read_tsv(paste0(TE_add,'/',file))[,c(TE_count)])),TE_reads)}
    colnames(TE_reads)[length(colnames(TE_reads))] <-  strsplit(file,"_")[[1]][1]}
  
  
  ############## gene_reads_arrange_merge\
  #C = (L x N) / G
  
  #where
  
  #C = coverage
  #L = read length (bp)
  #N = number of reads
  #G = haploid genome length (bp)
  
  
  
  #Exon
  
  gene_reads<-c() 
  #recursive files
  files <- list.files(gene_add,recursive=T)
  
  for (file in files) {
    if(length(gene_reads)==0){
      data <- read_tsv(paste0(gene_add,'/',file),skip=2,col_names = F)
      data <- data[which(data$X3=='exon'),c('X4','X5','X9')]   #
      data$X9 <- lapply(data$X9,function(x) strsplit(x,'\"')[[1]][10])
      
      data <- mapply(as.numeric,data)
      
      exon_reads<- data.frame(exon_reads = sum(data[,3] * (abs(data[,2]-data[,1])/reads_length)))
    sum(data[,3])
      
    }else{
      exon_reads <- cbind(data.frame(gene_reads = sum(read_tsv(paste0(gene_add,'/',file))[,c(3)])),gene_reads)}
    colnames(gene_reads)[length(colnames(gene_reads))] <-  strsplit(file,"_")[[1]][1]}
  
###################Transcripts
  gene_reads<-c() 
  #recursive files
  files <- list.files(gene_add,recursive=T)
  
  for (file in files) {
    if(length(gene_reads)==0){
      data <- read_tsv(paste0(gene_add,'/',file),skip=2,col_names = F)
      data <- data[which(data$X3=='transcript'),c('X4','X5','X9')]   #
      data$X9 <- lapply(data$X9,function(x) strsplit(x,'\"')[[1]][8])
      
      data <- mapply(as.numeric,data)
      
      exon_reads<- data.frame(exon_reads = sum(data[,3] * (abs(data[,2]-data[,1])/reads_length)))
      sum(data[,3])
      
    }else{
      exon_reads <- cbind(data.frame(gene_reads = sum(read_tsv(paste0(gene_add,'/',file))[,c(3)])),gene_reads)}
    colnames(gene_reads)[length(colnames(gene_reads))] <-  strsplit(file,"_")[[1]][1]}
  
  return(list(gene_reads,TE_reads))
}


all_seq_data_integration()