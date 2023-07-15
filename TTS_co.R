library(readr)
library(dplyr)
library(reshape2)
library(ChIPseeker) 
library(ggplot2)


co_reg_map <- function(DE_addr='co-regulation_map/DE_location'
                       ,save_addr='co-regulation_map/TTS_dis_FC/all.csv'){
  
  all_gene <- read_tsv(paste0(DE_addr,'/','all_DE_gene_co-regulation_map_locations.txt'),col_names = F) 
  
  all_TE <- read_tsv(paste0(DE_addr,'/','all_DE_TE_co-regulation_map_locations.txt'),col_names = F) 
  
  
  
  
  ###########
  #sort TE by log2(Fc) FC>2,p<10 -4
  
  
  
  #initialize genome map 
  position_array <- function(TE = all_TE, Gene = all_gene){
    
    TE <- TE[order(TE$X6,decreasing = T),]
    
    #initiation of co-array
    TE$clos_dis <-NA
    TE$clos_FC <- NA
    #for each DE_TEs

    for(rank in 1:length(TE$X1)){
      
      if(TE$X5[rank]=='+'){             
        corr_genes <- Gene[which(Gene$X2==TE$X2[rank] & Gene$X5== TE$X5[rank]),]  #chr and chain match
        corr_genes <- corr_genes[which(corr_genes$X3 > TE$X3[rank]),] #genes ahead of TE
        
        if(nrow(corr_genes)!=0){
          gene_end_distatnce <- corr_genes$X4 - TE$X3[rank]
          colst_genes <- corr_genes[which( gene_end_distatnce== max(gene_end_distatnce)),] #(find the closest gene)
          
          TE$clos_dis[rank] <- max(0,min(gene_end_distatnce)) #distance =0 if TE is inside gene
          
          TE$clos_FC[rank] <- mean(colst_genes$X6)} #if same distance, FC is the max FC
          
        
      
      }else{
        corr_genes <- Gene[which(Gene$X2==TE$X2[rank] & Gene$X5== TE$X5[rank]),]  #chr and chain match
        corr_genes <- corr_genes[which(corr_genes$X4 > TE$X4[rank]),] #genes ahead of TE
        
        if(nrow(corr_genes)!=0){
          gene_end_distatnce <- corr_genes$X3 - TE$X4[rank]
          colst_genes <- corr_genes[which( gene_end_distatnce== min(gene_end_distatnce)),] #(find the closest gene)
          
          TE$clos_dis[rank] <- min(0,min(gene_end_distatnce)) #distance =0 if TE is inside gene
          
          TE$clos_FC[rank] <- mean(colst_genes$X6)} #if same distance, FC is the max FC
    }}
    
    write_csv(TE,save_addr)
    
  }
  position_array(all_TE,all_gene)
  }
   

co_reg_map('co-regulation_map/DE_basemean','co-regulation_map/TTS_dis_basemean/all.csv')
co_reg_map('co-regulation_map/DE_FC','co-regulation_map/TTS_dis_FC/all.csv')

  
plot <- function(file='co-regulation_map/TTS_dis_FC/all.csv'){
  
  TE <- read_csv(file)
  
  TE_inside <- na.omit(TE[which(TE$clos_dis==0),])
  TE_ouside <- na.omit(TE[which(TE$clos_dis!=0),])
  TE_ouside$clos_dis <- log10(abs(TE_ouside$clos_dis))
  
  
  
  #############
  
  
  
  cor_test <- cor.test(TE_ouside$X6, TE_ouside$clos_FC, method = 'spearman',exact = FALSE)
  
  ann_text<-data.frame(
    x = max(TE_ouside$X6)*5/6, y = max(TE_ouside$clos_FC)*5/6,
    label = paste0('Rho: ',cor_test$estimate,
                   '\nP.value: ',cor_test$p.value,
                   '\nMethod: ',cor_test$method,
                   '\nPercent of outside/all_DE_TE:',nrow(TE_ouside)/nrow(TE))
  )
  
  
  ggplot(TE_ouside,aes(x=X6,y=clos_FC,color=clos_dis))+
    geom_point(alpha=0.5)+
    annotate("text", x=ann_text$x, y=ann_text$y, label=ann_text$label,
      color="black", size=2)+
    scale_colour_gradient(low='#000000',high='#ffffff')+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    annotate("text",x=1,y=1,
             label=paste0(round(nrow(TE_ouside[which(TE_ouside$clos_FC>0&TE_ouside$X6>0),])/nrow(TE_ouside)*100),'%'))+
    annotate("text",x=-1,y=1,
             label=paste0(round(nrow(TE_ouside[which(TE_ouside$clos_FC>0&TE_ouside$X6<0),])/nrow(TE_ouside)*100),'%'))+
    annotate("text",x=1,y=-1,
             label=paste0(round(nrow(TE_ouside[which(TE_ouside$clos_FC<0&TE_ouside$X6>0),])/nrow(TE_ouside)*100),'%'))+
    annotate("text",x=-1,y=-1,
             label=paste0(round(nrow(TE_ouside[which(TE_ouside$clos_FC<0&TE_ouside$X6<0),])/nrow(TE_ouside)*100),'%'))


  
  
  ##############
  
  
  cor_test <- cor.test(TE_inside$X6, TE_inside$clos_FC, method = 'spearman',exact = FALSE)
  
  ann_text<-data.frame(
    x = max(TE_inside$X6)*5/6, y = max(TE_inside$clos_FC)*5/6,
    label = paste0('Rho: ',cor_test$estimate,
                   '\nP.value: ',cor_test$p.value,
                   '\nMethod: ',cor_test$method,
                   '\nPercent of inside/all_DE_TE:',nrow(TE_inside)/nrow(TE))
  )
  
  ggplot(TE_inside,aes(x=X6,y=clos_FC))+
    geom_point(alpha=0.5,color='#000000')+
    annotate("text", x=ann_text$x, y=ann_text$y, label=ann_text$label,
             color="black", size=2)+
    geom_vline(xintercept=0)+
    geom_hline(yintercept=0)+
    
    annotate("text",x=1,y=1,
             label=paste0(round(nrow(TE_inside[which(TE_inside$clos_FC>0&TE_inside$X6>0),])/nrow(TE_inside)*100),'%'))+
    annotate("text",x=-1,y=1,
             label=paste0(round(nrow(TE_inside[which(TE_inside$clos_FC>0&TE_inside$X6<0),])/nrow(TE_inside)*100),'%'))+
    annotate("text",x=1,y=-1,
             label=paste0(round(nrow(TE_inside[which(TE_inside$clos_FC<0&TE_inside$X6>0),])/nrow(TE_inside)*100),'%'))+
    annotate("text",x=-1,y=-1,
             label=paste0(round(nrow(TE_inside[which(TE_inside$clos_FC<0&TE_inside$X6<0),])/nrow(TE_inside)*100),'%'))
  
 
  
}
  ####################







plo_abundance_enrichment<-function(file='co-regulation_map/TTS_dis_basemean/all.csv'){
  TE <- na.omit(read_csv(file))
  colnames(TE)[length(colnames(TE))] <-'clos_basemean'          # LOG2 basemean
  #TE$clos_basemean <- log2(TE$clos_basemean)
  #TE$X6 <- log2(TE$X6)
  
  TE_inside <- TE[which(TE$clos_dis==0),]
  TE_ouside <- TE[which(TE$clos_dis!=0),]
  
  TE_inside_percent <- sum(TE_inside$X6)/sum(TE$X6)
  gene_inside_percent <- sum(TE_inside$clos_basemean)/sum(TE$clos_basemean)
  
  TE_ouside_percent <- sum(TE_ouside$X6)/sum(TE$X6)
  gene_ouside_percent <- sum(TE_ouside$clos_basemean)/sum(TE$clos_basemean)
  
  
  
  
  
  
  
  ###################

  
  
  
  
  
  cor_test <- cor.test(TE_inside$X6, TE_inside$clos_basemean, method = 'spearman',exact = FALSE)
  
  
  TE_inside$X6 <- log2(TE_inside$X6)
  TE_inside$clos_basemean <- log2(TE_inside$clos_basemean)
  ann_text<-data.frame(
    x = max(TE_inside$X6)*5/6, y = max(TE_inside$clos_basemean)*5/6,
    label = paste0('Rho: ',cor_test$estimate,
                   '\nP.value: ',cor_test$p.value,
                   '\nMethod: ',cor_test$method,
                   '\nPercent of readcounts of inside/all_DE_TE:',TE_inside_percent)
  )
  
  ggplot(TE_inside,aes(x=X6,y=clos_basemean))+
    geom_point(alpha=0.5,color='#000000')+
    annotate("text", x=ann_text$x, y=ann_text$y, label=ann_text$label,
             color="black", size=2)
    
  
  
  ###################

  
  # spearmeam perfoemed before log2
  
  
  
  cor_test <- cor.test(TE_ouside$X6, TE_ouside$clos_basemean, method = 'spearman',exact = FALSE)
  
  TE_ouside$X6 <- log2(TE_ouside$X6)
  TE_ouside$clos_basemean <- log2(TE_ouside$clos_basemean)
  
  ann_text<-data.frame(
    x = max(TE_ouside$X6)*5/6, y = max(TE_ouside$clos_basemean)*5/6,
    label = paste0('Rho: ',cor_test$estimate,
                   '\nP.value: ',cor_test$p.value,
                   '\nMethod: ',cor_test$method,
                   '\nPercent of readcounts of ouside/all_DE_TE:',TE_ouside_percent)
  )
  
  ggplot(TE_ouside,aes(x=X6,y=clos_basemean))+
    geom_point(alpha=0.5,color='#000000')+
    annotate("text", x=ann_text$x, y=ann_text$y, label=ann_text$label,
             color="black", size=2)
  
  
  
  
}











    
    
    
    ####
   