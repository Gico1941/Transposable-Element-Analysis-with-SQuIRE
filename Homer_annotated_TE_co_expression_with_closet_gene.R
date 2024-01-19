library(readr)
library(tidyverse)




# significant DE TE (p<10^-4,log2FC>2)   cor with DE Genes (p<0.05)


DE_TE_peak_anno <- 'Homer/Annotatepeaks/annotate with FC/DE_TE_Homer.txt'

#### p_5  p_adj<0.05
co_exp <- function(peak_addr=DE_TE_peak_anno,
                   gene_addr = 'Homer/DE_location/all_DE_gene_Homer_symbol_list_p_5.txt'){
  
  peak_file <- read_delim(peak_addr,skip = 1,col_names = F,delim='\t')
  
  peak_file$X8 <- unlist(lapply(peak_file$X8,function(x) strsplit(x," \\(")[[1]][1]))
  
  gene_data <- read.table(gene_addr,header=T)
  
  
  data_merged <- merge(peak_file,gene_data,by.x='X16',by.y='Symbol',all.x = T)
  

  
  data_merged[,c('baseMean','log2FoldChange')][is.na(data_merged[,c('baseMean','log2FoldChange')]) ] <- 0
  
  
  ##plot 

  
  # create a list with a specific length 
  plot_lst <- vector("list", length = length(unique(data_merged$X8)))
  
  
  
  
  
  
  
  ######### plot by peak type ()UTR,INTRONS....
  
  
  for (i in 1:length(plot_lst)) {
    data_plot <- data_merged[which(data_merged$X8==unique(data_merged$X8)[i]),]
    
    
    cor_test <- cor.test(data_plot$X6, data_plot$log2FoldChange, method = 'spearman',exact = FALSE)
    
    ann_text<-data.frame(
      x = -max(data_merged$X6)*4/6, y = max(data_merged$log2FoldChange)*5/6,
      label = paste0('Rho: ',signif(cor_test$estimate,3 ),
                     '\nP.value: ',signif(cor_test$p.value, 3 ),
                     '\nMethod: ',cor_test$method,
                     '\nPercent of type of TE%:',signif(nrow (data_plot)*100/nrow(data_merged),3 ),
                     '\npercent of TE cor DE gene%:',signif(nrow (data_plot[which(data_plot$log2FoldChange!=0),])*100/nrow (data_plot),3) )
      
      ####### all are number percent
      
      #number of TEs with concordanct DE gene (p<0.05)
      #number of this type of TE/all TE
      
      
      
    )
    
    
    
    g <- ggplot(data_plot,aes(x=X6,y=log2FoldChange))+
      geom_point(alpha=0.3,color='orange')+
      theme_classic()+
      geom_vline(xintercept=0,color='gray',alpha=0.5)+
      geom_hline(yintercept=0,color='gray',alpha=0.5)+
      xlim(min(data_merged$X6)-1,max(data_merged$X6)+1)+
      ylim(min(data_merged$log2FoldChange)-1,max(data_merged$log2FoldChange)+1)+
      annotate("text",x=-1,y=-1,
               label=paste0(round(nrow(data_plot[which(data_plot$X6<0&data_plot$log2FoldChange<0),])/nrow(data_plot)*100),'%'))+
      annotate("text",x=1,y=-1,
               label=paste0(round(nrow(data_plot[which(data_plot$X6>0&data_plot$log2FoldChange<0),])/nrow(data_plot)*100),'%'))+
      annotate("text",x=-1,y=1,
               label=paste0(round(nrow(data_plot[which(data_plot$X6<0&data_plot$log2FoldChange>0),])/nrow(data_plot)*100),'%'))+
      annotate("text",x=1,y=1,
               label=paste0(round(nrow(data_plot[which(data_plot$X6>0&data_plot$log2FoldChange>0),])/nrow(data_plot)*100),'%'))+
      
      annotate("text", x=ann_text$x, y=ann_text$y, label=ann_text$label,
               color="black", size=2)+
      ggtitle(unique(data_plot$X8))+
      
      xlab('Log2FC TE(p<e-4,log2FC>2)')+
      ylab(paste0('Log2FC Cor DE gene(p<0.05)'))+
      theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
    
    
    
    
    
    plot_lst[[i]] <- g
  }
  
  cowplot::plot_grid(plotlist = plot_lst, nrow = 4)
  

  
}


co_exp()