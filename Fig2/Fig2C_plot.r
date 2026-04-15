#------------------------------------------------------------------------------------------------------------------

loMNase_fragL_dist<-fread("MAZ_loMNase_midP_fragL_dist.txt")

colnames(loMNase_fragL_dist)<-c("motif","fragment","distance")
loMNase_fragL_dist_list<-split(loMNase_fragL_dist,f=loMNase_fragL_dist$motif)

a<-80
dist<- 1
frag<- 8
v_outer_list<-lapply(loMNase_fragL_dist_list, function(x){
  temp<-x
  a<-80
  b<-20
  temp<-temp[between(temp$fragment,b,a),]
  
  v_inner_count<-sum(temp$fragment >= -2*temp$distance+frag+2*dist &
                       temp$fragment >= 2*temp$distance+frag-2*dist)

  v_outer_count<-sum((temp$fragment >= -2*(temp$distance-a)+frag+2*dist &
                        temp$distance <=dist+a)|
                       (temp$fragment >= 2*(temp$distance+a)+frag-2*dist &
                          temp$distance >=dist-a))
  
  
  v_outer_count<-ifelse(v_outer_count>0,v_outer_count,1)
  
  as.data.frame(cbind(v_inner_count=v_inner_count,
                      v_outer_count=v_outer_count,
                      FC=v_inner_count/v_outer_count,
                      motif=temp$motif[1]))
})


v_outer_list<-rbindlist(v_outer_list)


MAZ_signal_all<-fread("MAZ_20bp_motif_with_K562_signal.bed")

v_outer_list<-merge(v_outer_list,MAZ_signal_all[,c(4,7,8,9,10,11)],by.x="motif",by.y="V4")
v_outer_list$FC<-as.numeric(v_outer_list$FC)
v_outer_list$v_inner_count<-as.numeric(v_outer_list$v_inner_count)
v_outer_list$v_outer_count<-as.numeric(v_outer_list$v_outer_count)
v_outer_list<-v_outer_list[v_outer_list$FC>1,]

MAZ_1<-v_outer_list[between(v_outer_list$FC,1,3),]
MAZ_2<-v_outer_list[between(v_outer_list$FC,3,9),]
MAZ_3<-v_outer_list[between(v_outer_list$FC,9,27),]
MAZ_4<-v_outer_list[between(v_outer_list$FC,27,81),]
MAZ_5<-v_outer_list[between(v_outer_list$FC,81,571),]
mean(MAZ_1$FC)
mean(MAZ_2$FC)
mean(MAZ_3$FC)
mean(MAZ_4$FC)
mean(MAZ_5$FC)



loMNase_fragL_dist<-fread("MAZ_loMNase_midP_fragL_dist.txt")

colnames(loMNase_fragL_dist)<-c("motif","fragment","distance")
MAZ_FC_quan<-list(MAZ_1,MAZ_2,MAZ_3,MAZ_4,MAZ_5)

segment_table<-as.data.frame(cbind(x=c(81,81,81,37,37,-35,-43,-43,-79),
                                   y=c(80,80,8,80,80,80,80,80,8),
                                   xend=c(81,45,45,1,-35,1,-79,-79,-79),
                                   yend=c(8,80,80,8,80,8,80,8,80)))


for (i in 1:5) {
  temp_bind_fragL_dist<-loMNase_fragL_dist[loMNase_fragL_dist$motif%in%MAZ_FC_quan[[i]]$motif,]
  
  if (nrow(temp_bind_fragL_dist) >200000) {
    temp_bind_fragL_dist<-temp_bind_fragL_dist[sample(1:nrow(temp_bind_fragL_dist),200000),]
    
  }
  
  p<-ggplot(data =temp_bind_fragL_dist, aes(x=distance, y=fragment)) + 
    
    geom_pointdensity(size=0.2,show.legend = F)+
    scale_color_gradient2(low = "#88a7c5", high = "#104e8b") +
    # geom_abline(intercept = 8+2* 1, slope = -2,size = 0.5,linetype="dashed",color="red")+
    # geom_abline(intercept = 8-2* 1, slope = 2,size = 0.5,linetype="dashed",color="red")+
    geom_segment(data =segment_table,
      aes(x = x, y = y, xend =xend , yend = yend),
      size = 0.5,linetype="dashed",color="red"
    )+

    theme_bw()+xlab("Distance from motif (bp)")+
    ylab("Fragment length (bp)")+
    scale_x_continuous(
      limits = c(-85, 85),
      breaks=seq(-80, 80,40),
      expand = c(0,0)) +
    scale_y_continuous(
      limits = c(0, 100),
      # breaks=seq(0,150,50),
      expand = c(0,0)) +
    theme_classic()+
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
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_blank(),
          axis.ticks = element_line(color = "black"),
          text = element_text(size = 15, family = "sans", colour = "black"),
          plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm')) +coord_fixed(ratio = 1)
  
  
  ggsave(paste("MAZ_",i,"_vplot.png",sep = ""),p,width=5, height=3.5,dpi=300)
  print(i)
}


for (i in 1:5) {
  
  temp<-MAZ_signal_all[MAZ_signal_all$V4%in%MAZ_FC_quan[[i]]$motif,1:6]
  temp<-temp[order(temp$V1,temp$V2),]
  fwrite(temp,
         paste("MAZ_",i,".motif",sep = ""),
         col.names = F,row.names = F,quote = F,eol ="\n",sep = "\t" )
}





path<-"/mnt/disk2/2/zqf/private/GB_Vplot_Correspondence/TF_selected"


fileNames <- dir(path)
fileNames<-fileNames[grepl("MAZ_.*_NChIP.matrix.txt.gz",fileNames,perl = T)]
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})
computematrix_name<-gsub("_NChIP.matrix.txt.gz","",fileNames)

for (i in 1:length(filePath)) {
  computematrix<-fread(filePath[i],header = F)
  
  computematrix<-as.data.frame(cbind(
    signal=apply(computematrix[,7:206],2,function(y){
      mean(y[!is.na(y)])
    }),
    distance=seq(-2000,1980,20)))
  
  
  
  
  ggplot() + 
    geom_line(data=computematrix,aes(x=distance,y=signal,color="#104e8b"),size=1,alpha=0.8)+
    scale_y_continuous(limits = c(0,9),breaks = seq(0,8,4))+
    scale_color_jco()+
    theme_bw() +
    xlab("Distance from motif (bp)")+
    ylab("Enrichment")+
    theme(legend.box = "horizontal",
          legend.title = element_blank(),
          #legend.key.size = unit(18, "pt"),
          legend.position = "none",
          #panel.grid.major = element_blank(),
          legend.text = element_text(colour = "black", size = 15),
          legend.key.height = unit(0.01, 'pt'),
          legend.key.width = unit(20, 'pt'),
          panel.border = element_rect(fill = NA, size = 1, color = 'black'),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.y.right = element_text(color = "#CD5555", size = 15),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(colour = "black", size = 15, vjust = 0.5, hjust = 0.5),
          axis.ticks = element_line(size=1,color = "black"),
          text = element_text(size = 15, family = "sans", colour = "black"),
          plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'))
  ggsave(paste("/mnt/disk2/2/zqf/private/GB_Vplot_Correspondence/figure/",computematrix_name[i],"_lineplot.png"),width=5, height=3.5)
  
  
}
