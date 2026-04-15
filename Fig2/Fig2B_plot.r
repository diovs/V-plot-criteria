library(data.table)
library(ggplot2)
library(ggpointdensity)
library(ggpubr)
CTCF_WT<-fread("CTCF_loMNase_fragL_dist.txt")
colnames(CTCF_WT)<-c("motif","fragment","distance")
ggplot(data =CTCF_WT, aes(x=distance, y=fragment)) + 
  
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
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.y.right = element_text(color = "#CD5555", size = 15),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12, colour = "black"),
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 12, colour = "black"),
        plot.title = element_blank(),
        axis.ticks = element_line(color = "black"),
        text = element_text(size = 15, family = "sans", colour = "black"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm')) +coord_fixed(ratio = 1)



ggsave("loMNase_CTCF_no_bias_Vplot.png",width=3.5, height=3.5,dpi=300)





loMNase_fragment<-fread("K562_loMNase_merge_moredata.bed",header = F)
loMNase_fragment_length_tab<-as.data.frame(table(CTCF_WT$fragment))

loMNase_fragment_length_tab<-loMNase_fragment_length_tab[order(-loMNase_fragment_length_tab$Freq),]
loMNase_fragment_length_tab$Var1<-as.numeric(as.character(loMNase_fragment_length_tab$Var1))
loMNase_fragment_length_tab<-loMNase_fragment_length_tab[loMNase_fragment_length_tab$Var1>=20 &
                                                           loMNase_fragment_length_tab$Var1<=100,]
loMNase_fragment_length_tab$Freq<-loMNase_fragment_length_tab$Freq/sum(loMNase_fragment_length_tab$Freq)

fwrite(loMNase_fragment_length_tab,"loMNase_fragment_length_freq.txt",
       col.names = F,row.names = F,quote = F,eol ="\n",sep = "\t" )


result <- generate_ctcf_fragments(
  ctcf_file = "CTCF_selected_320bp.bed",
  freq_file = "loMNase_fragment_length_freq.txt",
  fragments_per_region = 1000,
  output_prefix = "CTCF_fragments_500per",
  seed = 2024
)



CTCF_random<-fread("CTCF_fragments_random_all_midP_fragL_CTCF_dist.txt")
colnames(CTCF_random)<-c("motif","fragment","distance")


CTCF_random<-CTCF_random[sample(1:nrow(CTCF_random),nrow(CTCF_WT)),]


ggplot(data =CTCF_random, aes(x=distance, y=fragment)) + 
  
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
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.y.right = element_text(color = "#CD5555", size = 15),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12, colour = "black"),
        axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 12, colour = "black"),
        plot.title = element_blank(),
        axis.ticks = element_line(color = "black"),
        text = element_text(size = 15, family = "sans", colour = "black"),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm')) +coord_fixed(ratio = 1)



ggsave("loMNase_CTCF_random_Vplot.png",width=3.5, height=3.5,dpi=300)
