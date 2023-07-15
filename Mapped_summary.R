# Load required packages
library(tidyverse)


# Define the path to the directory containing the sample folders


########### this is run on lab PC which has all log files

STAR_logs_integration <- function(path = "Mapped",
                                  decode = 'Call/GSE211061_coldata.txt',
                                  log_extension = "_fwd.log",
                                  addr='Counted/all/GSE211061_gene_subF_counttable.txt',
                                  combined_summary_from_STAR = 'Summary/combined.csv'
){
  
  
  
  
  decode <- read.table(decode,header=T)
  decode$sample <- lapply(decode$sample,function(x) strsplit(x,'_')[[1]][1])
  # Get a list of the sample folders in the directory
  sample_folders <- list.files(path, full.names = TRUE)
  
  # Initialize an empty list to store the data from each log file
  data_list <- list()
  
  # Loop through each sample folder and extract the data from the corresponding log file
  for (folder in sample_folders) {
    # Get the name of the sample from the folder name
    sample_name <- basename(folder)
    
    # Define the path to the log file
    log_file <- file.path(folder, paste0(sample_name,log_extension  ))
    
    
    # Read in the log file as a data frame
    log_data <- read.csv(log_file, sep = "|", strip.white = TRUE, header = FALSE )
    log_data <- log_data[which(log_data$V2!=''),]
    
    
    # Rename the columns of the data frame
    names(log_data) <- c("feature", "value")
    
    log_data$value[grep("%",log_data$value)] <- as.numeric(unlist(strsplit(log_data$value[grep("%",log_data$value)],'%')))
    
    
    # Add a column to the data frame for the sample name
    log_data$sample <- sample_name
    
    # Add the data frame to the list of data frames
    data_list[[sample_name]] <- log_data
  }
  
  # Combine all of the data frames into a single data frame
  combined_data <- bind_rows(data_list)
  
  # Convert the "value" column to numeric
  combined_data$value <- as.numeric(combined_data$value)
  combined_data <- na.omit(combined_data)
  
  #add condition info
  
  
  combined_data <- merge(combined_data,decode,all=T,by.x='sample',by.y='sample')
  
  ##### call TE Gene reads merge
  combined_data <- TE_gene_map_summary(combined_summary_from_STAR=combined_data)
  
  
  write.csv(combined_data,'Summary/Combined.csv')

  
}








######this is run locally to add mapp info from SQUIRE


TE_gene_map_summary <- function(addr='Counted/all/GSE211061_gene_subF_counttable.txt',
                                decode = 'Call/GSE211061_coldata.txt',
                                combined_summary_from_STAR='Summary/combined.csv',
                                save_addr='Summary/combined_.csv'){
  
  
  group <-  read.table(decode,header=T)
  
  reads <- read.table(addr,header=T) 
  
  combined <- read.csv(  combined_summary_from_STAR,header=T)[,-1]
  
  
  TE_reads <- data.frame(sample = gsub('_fwd','',group$sample),
                         feature ='number of reads mapped as TE',
                         value = lapply(group$sample,function(x) sum(reads[grep(':',reads$gene_id ),x])) %>% unlist(),
                         condition=group$condition)
  
  gene_reads <- data.frame(sample = gsub('_fwd','',group$sample),
                           feature ='number of reads mapped as Gene',
                           value = lapply(group$sample,function(x) sum(reads[-grep(':',reads$gene_id ),x])) %>% unlist(),
                           condition=group$condition)
  
  total_reads <- data.frame(sample = gsub('_fwd','',group$sample),
                               feature ='number of total mapped reads',
                            value = lapply(gsub('_fwd','',group$sample),function(x) sum(combined$value[which(combined$sample==x & 
                                                               combined$feature %in% c('Uniquely mapped reads number',
                                                                                       'Number of reads mapped to multiple loci'))])
                            ) %>% unlist(),
                               condition=group$condition)
  

  
  
  TE_reads_percent <- data.frame(sample = gsub('_fwd','',group$sample),
                         feature ='percent of reads mapped as TE %',
                         value = lapply(group$sample,function(x) (sum(reads[grep(':',reads$gene_id ),x])*100)
                                        /total_reads$value[which(total_reads$sample==gsub('_fwd','',x) )]
                                        ) %>% unlist(),
                         condition=group$condition)
  
  
  gene_reads_percent <- data.frame(sample = gsub('_fwd','',group$sample),
                                   feature ='percent of reads mapped as Gene %',
                                   value = lapply(group$sample,function(x) (sum(reads[-grep(':',reads$gene_id ),x])*100)
                                                  /total_reads$value[which(total_reads$sample==gsub('_fwd','',x) )]
                                   ) %>% unlist(),
                                   condition=group$condition)
  


  
  
  
  
  combine <- bind_rows(TE_reads,gene_reads,gene_reads_percent,TE_reads_percent,total_reads,combined )

 write.csv(combine ,save_addr)
}



TE_gene_map_summary()






plot <- function(combined='Summary/Combined_.csv') {
  
  combined_data<-read.csv(combined)
  
  
  for (feature in unique(combined_data$feature)) {
    # Filter the data for the current feature
    feature_data <- combined_data[which(combined_data$feature == feature),]
    
    # Create a ggplot object for the current feature
    plot <- ggplot(feature_data, aes(x = sample, y = value,color=condition)) +
      geom_point(size=5) +
      ggtitle(feature) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    # Save the plot as a pdf file
    ggsave(paste0('Summary/plot/',gsub("[^A-Za-z0-9 ]",' ',feature) , ".pdf"),plot)}
}

plot()








#STAR_logs_integration()











