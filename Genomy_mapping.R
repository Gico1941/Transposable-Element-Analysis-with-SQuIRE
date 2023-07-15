library(readr)
library(circlize)
library(dplyr)
#library(AneuploidyScore) downlaod cytoband of mm39
#ucsc.mm39.cytoband
library(ComplexHeatmap)



genome_position_file <- function(addr ='Call',
                            refseq_gene_addr = 'DESeq2_RefSeq_only.txt',
                            TE_gene_addr = 'DESeq2_TE_only.txt',
                            p_hold = 10**-4, 
                            FC_hold = 1,
                            location_file='gene_location.csv') {
  gene <- read.table(paste0(addr,'/',refseq_gene_addr),row.names = 1)
  loca <- read.csv(location_file)
  rownames(loca)<-paste0(loca$gene,',',loca$strand)
  gene <- gene[which(gene$padj<p_hold&abs(gene$log2FoldChange)>=FC_hold),]
  gene <- merge(loca,gene,by='row.names')
  
  
  TE <- read.table(paste0(addr,'/',TE_gene_addr), row.names = 1)
  #FILTER_OUT
  TE <- TE[which(TE$padj<p_hold&abs(TE$log2FoldChange)>=FC_hold),]
  TE[,c('chr',"start","end","TE","length",'strand')] <- t(data.frame(strsplit(rownames(TE),'\\|')))
  
  

  write_tsv(gene[,c('chr','start','end',"log2FoldChange",'padj')]
              ,"DE_gene_locations.txt",col_names = F)
  
  
  write_tsv(TE[,c('chr','start','end',"log2FoldChange",'padj')]
              ,"DE_TE_locations.txt",col_names = F)
}

DE_position_map_bar <- function(TE_add = "DE_TE_locations.txt" ,
                                gene_add = "DE_gene_locations.txt",
                                species ='mm39',down_stream=3000,up_stream=3000){
  
  #value1 = log2FC    value2 = padj
  
  TE <- read.table(TE_add,colClasses = c('character','numeric','numeric','numeric','numeric')) 
  colnames(TE) <-c('chr','start','end',"value1",'value2')
  
  gene <- read.table(gene_add,colClasses = c('character','numeric','numeric','numeric','numeric'))
  colnames(gene) <-c('chr','start','end',"value1",'value2')
  
  gene$start <- gene$start - up_stream
  gene$end <- gene$end + down_stream
  
  circos.initializeWithIdeogram(species = species)
  
  max = max(max(TE$value1),max(gene$value1))
  min = min(min(gene$value1),min(TE$value1))
  #col_fun_gene = colorRamp2(c(max(gene$value1), 0, min(gene$value1)), c("#75c8e1", "white", "#e67474"))
  col_fun = colorRamp2(c(max, 0, min), c("#FF000080", "white","#0000FF80" ),transparency = 0.1)
  
  lgd_links = Legend(at = c(round(min,2),round(min/2,2),0,round(max/2,2),round(max,2)), col_fun = col_fun, 
                     title_position = "topleft", direction = "horizontal", title = "Log2FoldChange")
  
  circos.genomicTrack(gene, 
                      panel.fun = function(region, value, ...) {
                        #print(col_fun(value[[1]]))
                        circos.genomicRect(region, value[[1]], ytop.column = 1 , ybottom = 0,col = col_fun(value[[1]]),border=col_fun(value[[1]]),...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  circos.genomicTrack(TE, 
                      panel.fun = function(region, value, ...) {
                        #print(col_fun(value[[1]]))
                        circos.genomicRect(region, value[[1]], ytop.column = 1 , ybottom = 0,col = col_fun(value[[1]]),border=col_fun(value[[1]]),...)
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                      })
  
  draw(lgd_links, x = unit(20, "mm"), y = unit(10, "mm"), just = c("left", "bottom"))
}

DE_position_map_density <- function(TE_add = "DE_TE_locations.txt" ,
                                    gene_add = "DE_gene_locations.txt",
                                    species ='mm39',down_stream=3000,up_stream=3000){
 ###################### 
  #value1 = log2FC    value2 = padj
  
  TE <- read.table(TE_add,colClasses = c('character','numeric','numeric','numeric','numeric')) 
  colnames(TE) <-c('chr','start','end',"value1",'value2')
  
  gene <- read.table(gene_add,colClasses = c('character','numeric','numeric','numeric','numeric'))
  colnames(gene) <-c('chr','start','end',"value1",'value2')
  
  gene$start <- gene$start - up_stream
  gene$end <- gene$end + down_stream
  
  circos.initializeWithIdeogram(species = species)
  
  max = max(max(TE$value1),max(gene$value1))
  min = min(min(gene$value1),min(TE$value1))
  #col_fun_gene = colorRamp2(c(max(gene$value1), 0, min(gene$value1)), c("#75c8e1", "white", "#e67474"))
  col_fun = colorRamp2(c(max, 0, min), c("#e67474", "white","#75c8e1" ),transparency = 0.1)
  
  lgd_links = Legend(at = c(round(min,2),round(min/2,2),0,round(max/2,2),round(max,2)), col_fun = col_fun, 
                     title_position = "topleft", direction = "horizontal", title = "Log2FoldChange")
  
  circos.genomicDensity(gene,col = c("#0000FF80"))
  circos.genomicDensity(TE,col = c("#FF000080"))
  
  
  draw(lgd_links, x = unit(20, "mm"), y = unit(10, "mm"), just = c("left", "bottom"))
  #circos.genomicDensity(TE,col = ifelse(area, "grey", "black"),baseline = 0)
  
  
  #circos.genomicRect(region, value, col, border)
  
  
  #bed = generateRandomBed(nr = 200)
  #circos.genomicTrack(bed, 
                      #panel.fun = function(region, value, ...) {
                       # circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,  ...)
                      #  circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#00000040")
                     # })
}


DE_position_map_Rainfall<- function(TE_add = "DE_TE_locations.txt" ,
                                    gene_add = "DE_gene_locations.txt",
                                    species ='mm39',down_stream=3000,up_stream=3000){
  ###################### 
  #value1 = log2FC    value2 = padj
  
  TE <- read.table(TE_add,colClasses = c('character','numeric','numeric','numeric','numeric')) 
  colnames(TE) <-c('chr','start','end',"value1",'value2')
  
  gene <- read.table(gene_add,colClasses = c('character','numeric','numeric','numeric','numeric'))
  colnames(gene) <-c('chr','start','end',"value1",'value2')
  
  gene$start <- gene$start - up_stream
  gene$end <- gene$end + down_stream
  
  circos.initializeWithIdeogram(species = species)
  
  max = max(max(TE$value1),max(gene$value1))
  min = min(min(gene$value1),min(TE$value1))
  
  circos.genomicRainfall(list(gene,TE), pch = 16, cex = 0.4, col = c("#FF000080", "#0000FF80"))
  
  lgd_links = Legend(at = c("Gene", "TE"), type = "points", 
                      legend_gp = gpar(col = c("#FF000080","#0000FF80")), title_position = "topleft", 
                      title = "Group")
  
  draw(lgd_links, x = unit(20, "mm"), y = unit(10, "mm"), just = c("left", "bottom"))
  
}

 

all_position_heatmap <- function(TE_add = "Counted/TE" , 
                                 gene_add = "Counted/Gene", 
                                 species ='mm39',
                                 gene_count='TPM', 
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
  
   TE_data <- cbind( TE_data[,c(2,3,4)],log2( TE_data[-c(1,2,3,4)]+0.01))
   colnames(TE_data)[c(1,2,3)]<- c('chr','start','end')
   TE_data <- TE_data[which(TE_data$chr %in% genomes_for_mapping),]
   colnames(TE_data)[-c(1:3)] <- factor(colnames(TE_data)[-c(1:3)],level=keys)
   
   
   ############## gene_data_arrange_merge
  
  gene_data<-c() 
  for (file in list.files(gene_add)) {
    if(length(gene_data)==0){
      gene_data <- read_tsv(paste0(gene_add,'/',file))[,c('Reference','Start','End','Gene ID',gene_count)] 
      colnames(gene_data)[length(colnames(gene_data))] <- strsplit(file,"_")[[1]][1]
    }else{
      
      gene_data <- merge(gene_data,read_tsv(paste0(gene_add,'/',file))[,c('Reference','Start','End','Gene ID',gene_count)],
                         by =c('Reference','Start','End','Gene ID'))
      
      colnames(gene_data)[length(colnames(gene_data))] <- strsplit(file,"_")[[1]][1]}}
  
  gene_data <- cbind(gene_data[,c(1,2,3)],log2(gene_data[-c(1,2,3,4)]+0.01))
  colnames(gene_data)[c(1,2,3)]<- c('chr','start','end')
  gene_data <- gene_data[which(gene_data$chr %in% genomes_for_mapping),]
  
  colnames(gene_data)[-c(1:3)] <- factor(colnames(gene_data)[-c(1:3)],level=keys) 
  gene_data$start <- gene_data$start - gene_up_stream
  gene_data$end <- gene_data$end + gene_down_stream

  ###########
  
  heat_plot_by_length <- function(data=gene_data,species=species,count_label=paste0('log2',gene_count),side='outside',col_data=col_data){
    ###
    max = max(data[,-c(1,2,3)])
    min = min(data[,-c(1,2,3)])
    
    col_fun = colorRamp2(c(min, 0, max), c("blue", "white","red" ),transparency = 0.5)
    ###
    lgd_links = Legend(at = c(round(min,2),round(min/2,2),0,round(max/2,2),round(max,2)), col_fun = col_fun, 
                       title_position = "topleft", direction = "horizontal", title = count_label)
    
    #circos.genomicHeatmap(data, col = col_fun, side = side , border = NA,connection_height=NULL)
    
    
                          
                          circos.genomicTrack(data, stack = TRUE, 
                                              panel.fun = function(region, value, ...) {
                                                
                                                print( mean(CELL_META$cell.ylim))
                                                
                                                circos.genomicRect(region, value, col = col_fun(value[[1]]), border = NA, ...)
                                                circos.lines(CELL_META$cell.xlim, c(mean(CELL_META$cell.ylim),
                                                                                    mean(CELL_META$cell.ylim)), lty = 2, col = "#00000040")
                                                
                                              })
                          
                        
    return(lgd_links)
  }
 
  heat_plot_by_regular <-function(data=gene_data,species=species,count_label=paste0('log2',gene_count),side='outside',col_data=col_data){
    ###
    max = max(data[,-c(1,2,3)])
    min = min(data[,-c(1,2,3)])
    
    col_fun = colorRamp2(c(min, 0, max), c("blue", "white","red" ),transparency = 0.5)
    ###
    lgd_links = Legend(at = c(round(min,2),round(min/2,2),0,round(max/2,2),round(max,2)), col_fun = col_fun, 
                       title_position = "topleft", direction = "horizontal", title = count_label)
    
    #circos.genomicHeatmap(data, col = col_fun, side = side , border = NA,connection_height=NULL)
    
    circos.genomicHeatmap(data, col = col_fun, side = side,connection_height=0.1)
    
    return(lgd_links)
  }

  circos.initializeWithIdeogram(species = species)
  ###################
  
  #gene_legend <- heat_plot_by_regular(gene_data[1:50,],species,paste0('Gene ','log2(',gene_count,"+0.01)"),'outside')
  #TE_legend <- heat_plot_by_regular(TE_data[1:50,],species,paste0('TE ','log2(',TE_count,"+0.01)"),'outside')

  
  ###########################
  TE_legend <- heat_plot_by_length(TE_data,species,paste0('TE ','log2(',TE_count,"+0.01)"),'outside')
  gene_legend <- heat_plot_by_length(gene_data,species,paste0('Gene ','log2(',gene_count,"+0.01)"),'outside')
  
  
  ##################
  legend =packLegend(gene_legend,TE_legend)
  
  draw(legend, x = unit(15, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))
  
  
  
  }




genome_position_file()


DE_position_map_Rainfall()
   
DE_position_map_bar()

all_position_heatmap()
  

  






circos.genoimcRainfall(bed)


circos.genomicRainfall(bed, pch = 16, cex = 0.4, col = c("#FF000080", "#0000FF80"))





bed = generateRandomBed(nr = 100, nc = 4)


circos.initializeWithIdeogram(species='mm39')















row.names(TE_PLOT)<- NULL



circos.initializeWithIdeogram(read.table('Book1.txt',sep='\t'))


c <- read.table('Book1.txt',sep='\t')


text(0, 0, "GSE210161_refseq", cex = 1)

circos.genomicDensity(TE_PLOT,col = c("#0000FF80"))
circos.genomicDensity(TE,col = c("#FF000080"))




circos.clear()



