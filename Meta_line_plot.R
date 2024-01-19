library(tidyverse)



meta_rna_peak <- read.table('Homer/Annotatepeaks/Meta annotate/all_TE_meta_profile.txt',skip=1)

  
  
ggplot(meta_rna_peak,aes(x=V1,y=V2))+
  geom_line(color='lightblue',size=1)+
  xlab('Distance from Center')+
  ylab('Peaks per gene per bp')+
  geom_vline(xintercept=0,color='gray',alpha=0.8,size=1)+
  geom_vline(xintercept=10000,color='gray',alpha=0.8,size=1)+
  theme_classic()+
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))+
  annotate("text", x=-1000, y=0, label="TSS", angle=0)+
  annotate("text", x=9000, y=0, label="TTS", angle=0)

ggsave('Homer/Annotatepeaks/Meta annotate/all_metaplot.pdf',width=10,height = 5)

  