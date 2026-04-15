library(data.table)
library(ggplot2)
library(ggpointdensity)
library(ggpubr)
CTCF_all_bind_fragL_dist<-fread("CTCF_loMNase_all_bind_fragL_dist.txt")

ggplot(data =CTCF_all_bind_fragL_dist, aes(x=distance, y=fragment)) + 
  
  geom_pointdensity(size=0.2,show.legend = F)+
  scale_color_gradient2(low = "#88a7c5", high = "#104e8b") +
  theme_bw()+xlab("Distance from motif (bp)")+
  ylab("Fragment length (bp)")+
  scale_x_continuous(
    limits = c(-60, 60),
    breaks=seq(-60,60,30),
    expand = c(0,0)) +
  scale_y_continuous(
    limits = c(0, 100),
    # breaks=seq(0,150,50),
    expand = c(0,0)) +
  theme_classic() +
  theme(legend.box = "horizontal",
        legend.title = element_blank(),
        legend.key.size = unit(18, "pt"),
        legend.position = "none",
        #panel.grid.major = element_blank(),
        legend.text = element_text(colour = "black", size = 15),
        legend.key.height = unit(0.7, 'cm'),
        panel.border = element_rect(fill = NA, size = 1, color = 'black'),
        axis.title.x = element_text(colour = "black", size = 20, vjust = 0),
        axis.title.y = element_text(colour = "black", size = 20, vjust = 2),
        axis.title.y.right = element_text(color = "#CD5555", size = 15),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 16, colour = "black"),
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 16, colour = "black"),
        plot.title = element_blank(),
        axis.ticks = element_line(color = "black"),
        text = element_text(size = 15, family = "sans", colour = "black"),
        plot.margin=unit(c(0.5,1.2,0.5,0.5),'cm')) +coord_fixed(ratio = 1)



ggsave("loMNase_CTCF_bound_Vplot.png",width=4.5, height=3.5,dpi=300)
