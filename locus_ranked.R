
library(tidyverse)



follow_same_order = F     ##############one of key parametetr

cumulative = T





all_TE_reads <- read_tsv('Counted/locus_TE_count.txt')

all_TE_reads <- all_TE_reads[-grep('\\.',all_TE_reads$symbol),]
decode <- read_tsv('Call/GSE211061_coldata.txt')

ctrl <- all_TE_reads[,c('symbol',decode$sample[which(decode$condition=='CTRL')])]
others <- all_TE_reads[,c('symbol',decode$sample[which(decode$condition !='CTRL')])]



ctrl$mean <- unlist(lapply(1:nrow(ctrl),function(x) mean(as.numeric(ctrl[x,-1]))))

others$mean <- unlist(lapply(1:nrow(others),function(x) mean(as.numeric(others[x,-1]))))


add_rank <- function(data){
  data <- data[order(data$mean),]
  data$rank <- as.numeric(row.names(data))
return(data)
}
  
if (follow_same_order== T){
  others  <- others [order(ctrl$mean),]
  others $rank <- as.numeric(row.names(others))
}else{
  others <- add_rank(others)
}

ctrl <- add_rank(ctrl)

if(cumulative ==T ){
  ggplot(others ,aes(x=rank,y=cumsum(mean) ,color='b'))+
    
    geom_point(alpha=0.1)+
    
    geom_point(data=ctrl,aes(x=rank,y=cumsum(mean),color='a') ,alpha=0.1)+
    
    scale_color_manual(name = 'group', 
                       values =c('b'='lightgreen','a'='grey'), 
                       labels = c('CTRL','LPS')) 
}else{
  ggplot(others ,aes(x=rank,y=log10(mean+1) ,color='b'))+
    
    geom_point(alpha=0.1)+
    
    geom_point(data=ctrl,aes(x=rank,y=log10(mean+1),color='a' ),alpha=0.1)+
    
    scale_color_manual(name = 'group', 
                       values =c('b'='lightgreen','a'='grey'), 
                       labels = c('CTRL','LPS')) 
}
  





percentage =F

if(percentage ==T ){
  
  
  ggplot(others ,aes(x=rank,y=cumsum(mean/sum(mean)) ,color='b'))+
    
    geom_point(alpha=0.3)+
    
    geom_point(data=ctrl,aes(x=rank,y=cumsum(mean/sum(mean)),color='a') ,alpha=0.3)+
    
    scale_color_manual(name = 'group', 
                       values =c('b'='lightgreen','a'='grey'), 
                       labels = c('CTRL','LPS')) 
  
  
}





  #scale_x_discrete(breaks=pretty(as.numeric(rev(others$rank)), n = 5),labels =pretty((rev(others$rank)), n = 5))
  
  
  #################  cumsum : cumulative plot
  


  
  #pretty(as.numeric(rev(others$rank)), n = 5)






