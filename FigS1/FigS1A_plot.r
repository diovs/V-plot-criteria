library(ggseqlogo)
library(data.table)
library(ggpubr)

#---------------------------------------DNase
fa<-fread("CTCF_DNase_all_bind.fa",header = F)
fa$V2<-toupper(fa$V2)


ggseqlogo(fa$V2, method="bits",seq_type="DNA" )+#ggtitle(name[i])+
  scale_y_continuous(limits = c(0,2))+
  #scale_x_continuous(expand = c(0.1,0),breaks=seq(1,30,3),limits =c(1,27) )+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_text(size=16,angle =0,hjust = 0.5), 
        axis.title.y=element_text(size = 20), 
        axis.title.x=element_text(size = 20),
        axis.line =  element_line(colour = "white",size = 4),
        panel.border = element_rect(fill=NA,color="white", size=1, linetype="solid"),
        legend.text=element_text(face="italic",  colour="black",  
                                 size=20),
        legend.title=element_text(face="italic",  colour="black", 
                                  size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        plot.margin=unit(c(0.5,1.2,0.5,0.5),'cm'))



ggsave("DNase_bias_CTCF_bound_weblogo.png", width=5, height=3.5, dpi = 300)
