library(tidyverse)
library(reshape2)
library(ggsignif)

addr = 'Mapped/chr_mapped'

decode=read_tsv('Call/GSE211061_coldata.txt')

mapped <- list.files(addr)


mapped.df <- lapply(mapped,function(x) read_delim(paste0(addr,'/',x),col_names = F,delim=' '))

 k<- mapped.df %>% reduce(full_join, by= 'X1')
 
 k[,-1] <- data.frame(lapply(colnames(k)[-1], function(x) k[,x]*100/sum(k[,x])))
 
 
 
 colnames(k) <- c('chr',unlist(gsub('.txt','',mapped)))
 
 c<-melt(k,id='chr')
 
 c <- merge(c,decode,by.x='variable',by.y='sample',all.x = T)
 
 
 
 
 ggplot(c[which(c$chr=='chrM'),],aes(x=condition,y=value,color=condition))+
   
   geom_boxplot(fill='white',alpha=0.5)+
   geom_jitter()+
   ylab('Mapped_chrM_reads %')+
   geom_signif(
     comparisons = list(c(unique(decode$condition))),
     map_signif_level = F,color='black'
   )
 
 
 
 
 
 
 
 
 
 
 
 
 
 

 
mt <- c[which(c$chr=='chrM'),]


colnames(mt)[c(1,2)] <- c('sample','feature')
mt$feature <- rep('percent of reads mapped to chrM %',nrow(mt))
mt$sample <- unlist(lapply(mt$sample,function(x) strsplit(as.character(x),"_")[[1]][1]))



TE_GENE <- read_csv("E:/bioinfo/TE_differentiated_gene_analysis_ribo_removed/Summary/combined_.csv")
   

TE_GENE <- TE_GENE[grep('percent',TE_GENE$feature),-1]





dt <- rbind(TE_GENE,mt)



ggplot(dt,aes(x=sample,y=value,color=condition,fill=feature))+
   geom_bar(stat='identity',alpha=0.5,size=1)+
   ylab('Mapped_chrM_reads %')+
  theme(axis.text.x = element_text(angle = 90))
 
 
 
 
 