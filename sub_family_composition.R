library(ggplot2)
library(readr)
library(dplyr)
library(sjmisc)
library(reshape2)
library(gplots)
library(gridExtra)
library(ggsignif)
library(grid)
library(cowplot)




data_process_merge <- function(addr ='Counted/Sub',
                     count='tot_reads',
                     normalize_by_library=T,
                     level='Subfamily',
                     save="family_compose/TE_compose_merged.csv"){
  

  data_process_write<-function(file){
    lines <- readLines(paste0(addr,'/',file))
    data <- lapply(lines, function(x) unlist(strsplit(x, ":|\t")))
    data <- as.data.frame(do.call("rbind", data))
    colnames(data)<-data[1,]
    data<-data[-1,]
    if(normalize_by_library==T){
      data[,count]<-as.numeric(data[,count])/as.numeric(data[,"aligned_libsize"])*10000}
    return <- data[,c('Subfamily',"Family","Class",count)]
    colnames(return)[4]<-data$Sample[1]
    return(return)}
 
  merged <- Reduce(function(x, y) merge(x, y, ,by=c('Subfamily',"Family","Class"),all=T), 
             lapply(list.files(addr),data_process_write))
  write_csv(merged,save)
  
}



###############################




plots_overlay_bar <- function(
                  compose_file="family_compose/TE_compose_merged.csv",
                  level='Subfamily',
                  decode='Call/GSE211061_coldata.txt',
                  top=F,
                  overview_level ='Family'
                  ){
  
  decode <- read_table(decode)
  
  
  
  #read_file and process
  data <- read_csv(compose_file)
  data <- aggregate(data[,-c(1:3)],data[,level],sum)
  colnames(data)[-1] <- factor(colnames(data)[-1],levels = decode$sample)
  row.names(data) <- data[,1]
  data<- data[,-1]
  
  # sort by decreasing sum
  data$sum <- unlist(lapply(c(1:length(data$`1`)),function(x) sum(data[x,])))
  data <- data[order(data$sum,decreasing = T),]
  data <- data[,-c(length(colnames(data)))]
  
  
    
  #heatmap(as.matrix(data))
  #legend("topright", legend = seq(min(data), max(data), length.out = 10), 
         #fill = rev(cm.colors(256)), bty = "n")
  
  #select top 200 high TE
  if(top != F){
    data<-data[c(1:top),]
    percent <- lapply(unique(decode$condition),function(x) paste0("top ",top," Rank TEs in average accounts for ",
                                                                  round(sum(data[c(1:top),which(decode$condition==x)])
                                                                        /sum(data[,which(decode$condition==x)])*100,2),
                                                                  "% of total reads in ",x))
  }
  
  
  
  #melt data for overlay bar plot
  data$row <- seq_len(nrow(data))  
  dat <- melt(data, id.vars = "row")
  dat$group <-c(rep("LPS",length(dat$row)/2),rep("Ctrl",length(dat$row)/2))


  ##########   overlay
  ggplot(dat, aes(x =row , y = value, fill = group)) + 
   # scale_fill_gradient(low="darkorange",high="pink")+
    geom_bar(stat = "identity",position = "identity",alpha=0.3) +
    scale_fill_manual(values=c("#eae32e","#2ed2ea")) +
    xlab("TE_rank_by_reads") +
    ylab("TE_reads/tot_reads*10000,calculate per sample") +
    annotate("text",x=top/2,y=150,
             label=paste(unlist(percent),collapse='\n'))+
    theme_bw()
  
  #######   same plot but stack not overlay
  ggplot(dat, aes(x =row , y = value, fill = group)) + 
  # scale_fill_gradient(low="darkorange",high="pink")+
  geom_bar(stat = "identity",alpha=0.3) +
    scale_fill_manual(values=c("#eae32e","#2ed2ea")) +
    xlab("TE_rank_by_reads") +
    ylab("TE_reads/tot_reads*10000,calculate per sample") +
    annotate("text",x=top/2,y=150,
             label=paste(unlist(percent),collapse='\n'))+
    theme_bw()
  
  
}
  





########################################








plots_percentage <- function(
    compose_file="family_compose/TE_compose_merged.csv",
    decode='Call/GSE211061_coldata.txt',
    overview_level ='Family',
    top=F #number of false
){


  ############################################################################################
  #### Family/Class overview plot
  
  
  decode <- read_table(decode)
  
  
  
  data <- read_csv(compose_file)
  
  #sum read_counts based on levels
  data <- aggregate(data[,-c(1:3)],data[,overview_level],sum)
  
  

  
  #remove family with "?"
  undefined_names <- which(unlist(lapply(data[,overview_level], function(x) grepl("\\?",x))))
  if(length(undefined_names) != 0){
    data <- data[-c(undefined_names),]
  }
    
  ####sort data by TE expressiob level
  data$sum <- unlist(lapply(c(1:length(data[,1])),function(x) sum(data[x,-1])))
  data <- data[order(data$sum,decreasing = T),]
  data <- data[,-c(length(colnames(data)))]
  
  #keep already libsize-normalized total expression as total (total of tops)
  
  if(top !=F){
    total <- data.frame(value = unlist(lapply(colnames(data)[-1], function(x) sum(data[1:top,x]))),
                        sample = colnames(data)[-1])
  }else{
    total <- data.frame(value = unlist(lapply(colnames(data)[-1], function(x) sum(data[,x]))),
                        sample = colnames(data)[-1])
  }

  
  total <- merge(total,decode,by='sample')
  total$sample<- factor(total$sample,levels = decode$sample)
  
  
  
  #return percentage
  data[,-1] <-data.frame(lapply(colnames(data)[-1], function(x) data[,x]/sum(data[,x])*100))
  
  #data [level,]
  
 
  colnames(data)[1] <- "level"
  

  

  #  ## calculate percentage of tops 
  
  
  if (top  != F){
    percent <- lapply(unique(decode$condition),function(x) paste0("top ",top," Rank TEs in average accounts for ",
                                                                  round(sum(data[,-1][c(1:top),which(decode$condition==x)])
                                                                        /sum(data[,-1][,which(decode$condition==x)])*100,2),
                                                                  "% of total TE reads in ",x))
  }
  
  
  # take top xx of data              
  if (top  != F){

    data<-data[c(1:top),]
    
  }
  

  ##melt for plot and add annotation
  dt <- melt(data,id.vars='level')
  dt$group <-c(rep(unique(decode$condition)[1],
                   length(dt[,1])/2),rep(unique(decode$condition)[2],length(dt[,1])/2))

  dt$variable <- factor(dt$variable,levels = decode$sample)
  
  

  
  
  
  ############### histogram
  his <-ggplot(dt, aes(x =variable, y =value, fill = level)) + 
  theme_classic()+
    
  geom_bar(stat = "identity")+
  ylab(paste("TE_percentage_by",overview_level))+
  xlab("Sample") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  
  ###############  1. heat map

  hm <- ggplot(data = dt, aes(x = factor(variable), y = reorder(level,value,decreasing = F), fill = value)) + geom_tile() + 
    scale_fill_distiller(name = "TE percentage of all TE", palette = "Reds", direction = 1, na.value = "transparent") +
    scale_x_discrete(breaks = unique(dt$variable), labels = unique(dt$variable)) + theme_classic() +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.title = element_text(size = 15), legend.key.size = unit(1,"cm"),
          legend.text = element_text(size = 7)) +
    guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))+
    
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    #theme(axis.text.y=element_blank())
  

  ###### 1.2.extract map legend
  tmp <- ggplot_gtable(ggplot_build(hm))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  
  ##### 1.3 cleaned heatmap for aggregate
  
  if(overview_level=='Subfamily'){
    hm.clean <- hm +
    theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(), 
          axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position="none")+
      theme(plot.margin = unit(c(0,0,1,1), "cm"))+
      
      geom_vline(xintercept = (2:length(decode$sample)-1)+0.5,color='gray',alpha=0.3)+
      
      guides(x = "none")+
    
    guides(y="none")
    
  }else{
    
    hm.clean <- hm +
      
      theme(axis.title.y = element_blank(),#axis.text.y = element_blank(),axis.ticks.y = element_blank(), 
            axis.title.x = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            legend.position="none")+
      theme(plot.margin = unit(c(0,0,1,1), "cm"))+
      
      
      geom_vline(xintercept = (2:length(decode$sample)-1)+0.5,color='gray',alpha=0.3)+
      
      guides(x = "none")#+
    
    #guides(y="none")
    
    
  }
 

  

   
  ###### 2 heatmap _ y overlay bar (each TE of each individual/sample overlay )
  over_bar <- ggplot(dt, aes(x =reorder(level,value,decreasing = F) , y =value , fill = group)) + 
    theme_classic()+
    # scale_fill_gradient(low="darkorange",high="pink")+
    geom_bar(stat = "identity",position = "identity",alpha=0.3) +
    scale_fill_manual(values=c("#ffff66","#00ccff")) +
    xlab(paste("TE",overview_level)) +
    ylab(paste("TE_percentage_by",overview_level)) +
    coord_flip() +
    
    theme(plot.margin = unit(c(0,0,0.5,0), "cm"))+
    labs(x = "Overlay Percentage of tops for each sample")+
    scale_x_discrete(position = "left")+
    #scale_y_continuous(position = 'left')+
    guides(y="none")
    #guides(x = "none")
  
  
  over_bar.clean <- over_bar +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size=10),
                                   axis.ticks.x = element_blank(), axis.text.y = element_blank(), 
                                   axis.title.y = element_text(size = 10, margin = margin(0,10,0,0), angle = -90),
                                   legend.position="right")
  
  #element_text(size=xx)
  ##### 3 heatmap _ x normal bar ( total expression of all TE)
  
  total_bar <- ggplot(total,aes(x=sample,y=value,fill=condition))+
    geom_bar (stat = "identity") + theme_classic() +
    theme(axis.title.y = element_blank(), axis.text.y = element_text(size=10), 
          axis.ticks.y = element_blank(), axis.text.x = element_blank(), 
          axis.title.x = element_text(size = 10, margin = margin(10,0,0,0)),
          legend.position = "top") +
    #scale_fill_distiller(name = "Value", palette = "Reds", direction = 1) + 
    labs(x = "Total TE reads of tops per sample/aligned_libsize*10^4")+
    guides(x = "none")+
    #guides(y = "none")+
    theme(plot.margin = unit(c(0,0,0.2,ifelse(overview_level=='Subfamily',0.2,
                                       ifelse(overview_level=='Family',2.4,2.0))), "cm"))+
    
    scale_x_discrete(position = "top")
    #scale_y_discrete(position = "right")
    #+
    #guides(x = "none")

  
  ##### merge all three plots
  

  #########

  
  ##########
 if(top !=F){
   grob.title <- textGrob(paste(unlist(percent),collapse='\n'), hjust = 0.5, vjust = 0.5, gp = gpar(fontsize = 10))
 }else{
   grob.title <- textGrob(paste0("TE_",overview_level,'(tops=',top,')' ), hjust = 0.5, vjust = 0.5, gp = gpar(fontsize = 10))
 }
  
  final <- grid.arrange(total_bar, legend, hm.clean, over_bar.clean, nrow = 2, ncol = 2, 
               widths = c(30, 20), heights = c(40, 60), top = grob.title)
  
  ggsave(paste0("family_compose/heatmap_",overview_level,"top",top,'.pdf'),plot=hm,width=10, height=8)
  ggsave(paste0("family_compose/barplot_",overview_level,"top",top,'.pdf'),plot=over_bar,width=10, height=8) 
  return(final)


}

plots_percentage(overview_level ='Subfamily',)   #'Class' 'Subfamily' 'Family'

## to adjust margin of total_bar, search  theme(plot.margin = unit(c(0,0,0.2,ifelse(overview_level=='Subfamily',0.2,2)), "cm"))+


data_process_merge()
plots()






different_level_bars <- function(
    compose_file="family_compose/TE_compose_merged.csv",
    level='Family',
    decode='Call/GSE211061_coldata.txt',
    top=15,
    relative =F,
    percentage_mode=F
  ){
    
    decode <- read_table(decode)
    
    
    
    #read_file and process
    data <- read_csv(compose_file)
    data <- aggregate(data[,-c(1:3)],data[,level],sum)
    colnames(data)[-1] <- factor(colnames(data)[-1],levels = decode$sample)

    
    
    # this return percentage
    if (percentage_mode==T){data[-1] <- data.frame(lapply(colnames(data)[-1], function(x) data[,x]/sum(data[,x])*100))}
    #
    
    
    
    row.names(data) <- data[,1]
    data<- data[,-1]
    
    # sort by decreasing sum
    data$sum <- unlist(lapply(c(1:length(data$`1`)),function(x) sum(data[x,])))
    data <- data[order(data$sum,decreasing = T),]
    data <- data[,-c(length(colnames(data)))]
    
    
    
    
    data<-data[-grep('RNA',rownames(data)),]

    
    #select top 200 high TE
    if(top != F){
      data<-data[c(1:top),]
      percent <- lapply(unique(decode$condition),function(x) paste0("top ",top," Rank TEs in average accounts for ",
                                                                    round(sum(data[c(1:top),which(decode$condition==x)])
                                                                          /sum(data[,which(decode$condition==x)])*100,2),
                                                                    "% of total reads in ",x))
    }else{
      percent ='all'
    }
    
    
    
    #melt data for overlay bar plot
   
    

      
      data$row <-  row.names(data)                  #seq_len(nrow(data))  
      dat <- melt(data, id.vars = "row")
      dat$group <-unlist(lapply(dat$variable,function(x) decode$condition[x] ))
      
      
      #### enable it to create relative expression level
      if (relative ==T){  for(i in unique(dat$row)){
        dat[which(dat$row==i),'value'] <- dat[which(dat$row==i),'value']*length(which(dat$row==i))/sum(dat[which(dat$row==i),'value'])
      }}
    
      
      
      
      
      ggplot(dat, aes(x =group , y = value, color = group)) + 
        # scale_fill_gradient(low="darkorange",high="pink")+
        geom_boxplot(fill='white') +
        geom_jitter()+
        scale_fill_manual(values=c("#eae32e","#2ed2ea")) +
        ylab("") +
        # annotate("text",x=top/2,y=150,
        #label=paste(unlist(percent),collapse='\n'))+
        theme_bw()+
        ggtitle('TE Expression ')+
        #theme(legend.position="none")+
        facet_wrap(~row, scales = "free_x",nrow = 4)+
        stat_summary(fun=mean, colour="black", geom="point", 
                     shape=18, size=2, show.legend=FALSE)
      #geom_signif( comparisons= list(unique(decode$condition)),y_position=300)
    
    

    ##########   overlay
    

  


    
  }


different_level_bars(relative =T)









#colors = colorRampPalette(c("#e6581a", '#37e2d2'), space = "rgb")(2)
#barplot(as.matrix(t(data[c(1:200),c(1:6)])),col=rgb(0,0,1,alpha = 0.5),xlab = 'Samples',ylab = '(TE_reads/libsize)*10000',legend.text = "LPS")
#mtext(paste0("top 200 accounts for ",sum(data[c(1:200),c(1:6)])/sum(data[,c(1:6)]),"of total reads in LPS"))


#barplot(as.matrix(t(data[c(1:200),c(7:12)])),col=rgb(47,203,225,alpha = 0.5),xlab = 'Samples',ylab = '(TE_reads/libsize)*10000')
#mtext(paste0("top 200 accounts for ",sum(data[c(1:200),cc(7:12)])/sum(data[,c(7:12)]),"of total reads in Ctrl"),legend.text = "Ctrl")