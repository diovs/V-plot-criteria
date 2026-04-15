library(data.table)
library(ggplot2)
library(ggpubr)
#-----------------------------MNase
computematrix<-fread("K562_loMNase_5_3_no_bind.matrix.txt.gz",header = F)

computematrix<-as.data.frame(rbind(cbind(
  signal=apply(computematrix[,7:126],2,function(y){
    mean(y[!is.na(y)])
  }),
  distance=seq(-60,59,1),
  type="5'"),
  cbind(
    signal=apply(computematrix[,127:246],2,function(y){
      mean(y[!is.na(y)])
    }),
    distance=seq(-60,59,1),
    type="3'"))
)

library(ggsci)
computematrix$signal<-as.numeric(computematrix$signal)
computematrix$distance<-as.numeric(computematrix$distance)
computematrix$type<-factor(computematrix$type,levels = c("5'","3'"))


ggplot() + 
  geom_line(data=computematrix,aes(x=distance,y=signal*1000,color=type),size=1,alpha=0.8)+
  scale_y_continuous(limits = c(0,350))+
  scale_color_jco()+
  theme_bw() +
  xlab("Distance from motif (bp)")+
  ylab("Enrichment (10-3)")+
  theme(legend.box = "horizontal",
        legend.title = element_blank(),
        #legend.key.size = unit(18, "pt"),
        legend.position = "none",
        #panel.grid.major = element_blank(),
        legend.text = element_text(colour = "black", size = 15),
        legend.key.height = unit(0.01, 'pt'),
        legend.key.width = unit(20, 'pt'),
        panel.border = element_rect(fill = NA, size = 1, color = 'black'),
        axis.title.x = element_text(colour = "black", size = 20, vjust = 0),
        axis.title.y = element_text(colour = "black", size = 20, vjust = 2),
        axis.title.y.right = element_text(color = "#CD5555", size = 15),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 16, colour = "black"),
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 16, colour = "black"),
        plot.title = element_text(colour = "black", size = 15, vjust = 0.5, hjust = 0.5),
        axis.ticks = element_line(color = "black"),
        text = element_text(size = 15, family = "sans", colour = "black"),
        plot.margin=unit(c(0.5,1.2,0.5,0.5),'cm'))
ggsave("loMNase_fragment_end_CTCF_unbound_lineplot.pdf",width=4.5, height=3.5)

