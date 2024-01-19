library(tidyverse)

tops = 20


DE_TE <- read_tsv("E:/bioinfo/TE_differentiated_gene_analysis_ribo_removed/Homer/find_motif/All DEs/DE_TE/knownResults.txt")

DE_Gene <- read_tsv("E:/bioinfo/TE_differentiated_gene_analysis_ribo_removed/Homer/find_motif/All DEs/DE_gene/knownResults.txt")


plot <- function(data=DE_Gene,level='TE',tops=tops ){
  
  data[,7] <- unlist( lapply(unlist(data[,7]), function(x) as.numeric(strsplit(x,'%')[[1]][1])))
  
  data[,9] <- unlist(lapply(unlist(data[,9]), function(x) as.numeric(strsplit(x,'%')[[1]][1])))
  
  
  data[,'S/N ratio'] <- data[,7]/data[,9]
  data <- data[which(data $`q-value (Benjamini)`<0.5),]
  
  data <- data[order(data$`S/N ratio`,decreasing = F),]
  

  
  data$`Motif Name` <- factor(data$`Motif Name`,levels = data$`Motif Name`)
  
  ggplot(data %>% top_n(-tops,`q-value (Benjamini)`),aes(x=`S/N ratio`,y=`Motif Name`,size=`% of Target Sequences with Motif`,color=`q-value (Benjamini)`))+
    geom_point()+
    scale_color_gradient(low='#8db3e6',high='#1c4884')+
  guides(colour = guide_colourbar(order = 1),
         size = guide_legend(order = 2))
         #size = guide_legend(order = 3))
  
  ggsave(paste0("E:/bioinfo/TE_differentiated_gene_analysis_ribo_removed/Homer/find_motif/All DEs/",
                level,".pdf"))
  
}



plot(DE_TE,'TE',tops=20 )

plot(DE_Gene,'Gene',tops=20 )