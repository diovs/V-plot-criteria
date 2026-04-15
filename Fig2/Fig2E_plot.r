cl <- makeSOCKcluster(80)
registerDoSNOW(cl)


path<-"/mnt/disk2/2/zqf/private/GB_Vplot_Correspondence/DNase_Bias_3bp_sequence"


fileNames <- dir(path)
fileNames<-fileNames[grepl("_DNase_fragL_dist.txt",fileNames)]
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})


x_crosspoint<-seq(-25,25,1)
pb <- txtProgressBar(max = length(x_crosspoint), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)




bias_name<-gsub("_DNase_fragL_dist.txt","",fileNames)
bias_result_list<-list()
for (i in 1:length(filePath)) {
  temp<-fread(filePath[i])
  colnames(temp)<-c("motif","fragment","distance")
  temp<-temp[between(temp$fragment,30,100),]
  
  result <-foreach(j = 1:length(x_crosspoint),
                   .packages = c("data.table"),
                   .options.snow=opts) %dopar% {
                     options(scipen = 20)
                     dist<-x_crosspoint[j]
                     frag<-20
                     delta_x_intercept<-3
                     
                     Pos_2_inner_points <- sum(temp$fragment >= 2 * temp$distance + frag - 2 * dist &
                                                 temp$fragment <= 2 * (temp$distance + delta_x_intercept) + frag - 2 * dist)
                     
                     Pos_2_outer_points <- sum(temp$fragment <= 2 * temp$distance + frag - 2 * dist &
                                                 temp$fragment >= 2 * (temp$distance - delta_x_intercept) + frag - 2 * dist)
                     
                     Neg_2_inner_points <- sum(temp$fragment >= -2 * temp$distance + frag + 2 * dist &
                                                 temp$fragment <= -2 * (temp$distance - delta_x_intercept) + frag + 2 * dist)
                     
                     Neg_2_outer_points <- sum(temp$fragment <= -2 * temp$distance + frag + 2 * dist &
                                                 temp$fragment >= -2 * (temp$distance + delta_x_intercept) + frag + 2 * dist)
                     
                     Pos_2_inner_points <- ifelse(Pos_2_inner_points > 0, Pos_2_inner_points, 1)
                     Pos_2_outer_points <- ifelse(Pos_2_outer_points > 0, Pos_2_outer_points, 1)
                     Neg_2_inner_points <- ifelse(Neg_2_inner_points > 0, Neg_2_inner_points, 1)
                     Neg_2_outer_points <- ifelse(Neg_2_outer_points > 0, Neg_2_outer_points, 1)
                     
                     Pos_2_inner_outer_ratio <- Pos_2_inner_points / Pos_2_outer_points
                     
                     Neg_2_inner_outer_ratio <- Neg_2_inner_points / Neg_2_outer_points
                     
                     
                     as.data.frame(cbind(dist, frag,
                                         Pos_2_inner_outer_ratio, Neg_2_inner_outer_ratio,
                                         Pos_2_inner_points, Pos_2_outer_points,
                                         Neg_2_inner_points, Neg_2_outer_points
                     ))
                   }
  
  bias_result_list[[i]] <- rbindlist(result)
  bias_result_list[[i]]$TF<-bias_name[i]
  
}
bias_result_list_1<-list()

for (i in 1:length(bias_result_list)) {
  temp<-bias_result_list[[i]]
  
  temp<-temp[order(-temp$Pos_2_inner_outer_ratio),]
  dist1<-mean(temp$dist[1:3])
  #dist1<-temp$dist[1]
  
  temp<-temp[order(-temp$Neg_2_inner_outer_ratio),]
  dist2<-mean(temp$dist[1:3])
  #dist2<-temp$dist[1]
  
  temp<-as.data.frame(cbind(dist1,frag1=20,dist2,frag2=20))
  
  temp$TF<-bias_result_list[[i]]$TF[i]
  bias_result_list_1[[i]]<-temp
}




bias_result_list_1<-rbindlist(bias_result_list_1)




bias_result_list_1$v_smt_dist<-(bias_result_list_1$dist1+bias_result_list_1$dist2)/2

bias_result_list_1$v_smt_frag<- round(2*bias_result_list_1$v_smt_dist+20-2*bias_result_list_1$dist1)

bias_result_list_1<-bias_result_list_1[bias_result_list_1$TF%in%c("TGC","TCC","TGA","CCC","TGG"),]



v_inner_outer_FC<-list()

for (i in 1:nrow(bias_result_list_1)) {
  dist<-bias_result_list_1$v_smt_dist[i]
  frag<-bias_result_list_1$v_smt_frag[i]
  
  temp<-fread(filePath[grepl(bias_result_list_1$TF[i],filePath)])
  
  colnames(temp)<-c("motif","fragment","distance")
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
  
  v_inner_outer_FC[[i]]<-as.data.frame(cbind(v_inner_count=v_inner_count,
                                             v_outer_count=v_outer_count,
                                             FC=v_inner_count/v_outer_count,
                                             TF=bias_result_list_1$TF[i]))
  
  
}





v_inner_outer_FC<-rbindlist(v_inner_outer_FC)          


bias_result_list_2<-merge(bias_result_list_1,v_inner_outer_FC,by="TF",all=T)
bias_result_list_2$FC<-as.numeric(bias_result_list_2$FC)


path<-"/mnt/disk2/2/zqf/private/GB_Vplot_Correspondence/DNase_selected/"


fileNames <- dir(path)
fileNames<-fileNames[grepl("_midP_fragL_dist.txt",fileNames)]
filePath <- sapply(fileNames, function(x){ 
  paste(path,x,sep='/')})




x_crosspoint<-seq(-25,25,1)
pb <- txtProgressBar(max = length(x_crosspoint), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)




TF_name<-gsub("_midP_fragL_dist.txt","",fileNames)
result_list<-list()
for (i in 1:length(filePath)) {
  temp<-fread(filePath[i])
  colnames(temp)<-c("motif","fragment","distance")
  #temp<-temp[between(temp$fragment,20,75),]
  
  result <-foreach(j = 1:length(x_crosspoint),
                   .packages = c("data.table"),
                   .options.snow=opts) %dopar% {
                     options(scipen = 20)
                     dist<-x_crosspoint[j]
                     frag<-20
                     delta_x_intercept<-5
                     
                     Pos_2_inner_points <- sum(temp$fragment >= 2 * temp$distance + frag - 2 * dist &
                                                 temp$fragment <= 2 * (temp$distance + delta_x_intercept) + frag - 2 * dist)
                     
                     Pos_2_outer_points <- sum(temp$fragment <= 2 * temp$distance + frag - 2 * dist &
                                                 temp$fragment >= 2 * (temp$distance - delta_x_intercept) + frag - 2 * dist)
                     
                     Neg_2_inner_points <- sum(temp$fragment >= -2 * temp$distance + frag + 2 * dist &
                                                 temp$fragment <= -2 * (temp$distance - delta_x_intercept) + frag + 2 * dist)
                     
                     Neg_2_outer_points <- sum(temp$fragment <= -2 * temp$distance + frag + 2 * dist &
                                                 temp$fragment >= -2 * (temp$distance + delta_x_intercept) + frag + 2 * dist)
                     
                     Pos_2_inner_points <- ifelse(Pos_2_inner_points > 0, Pos_2_inner_points, 1)
                     Pos_2_outer_points <- ifelse(Pos_2_outer_points > 0, Pos_2_outer_points, 1)
                     Neg_2_inner_points <- ifelse(Neg_2_inner_points > 0, Neg_2_inner_points, 1)
                     Neg_2_outer_points <- ifelse(Neg_2_outer_points > 0, Neg_2_outer_points, 1)
                     
                     Pos_2_inner_outer_ratio <- Pos_2_inner_points / Pos_2_outer_points
                     
                     Neg_2_inner_outer_ratio <- Neg_2_inner_points / Neg_2_outer_points
                     
                     
                     as.data.frame(cbind(dist, frag,
                                         Pos_2_inner_outer_ratio, Neg_2_inner_outer_ratio,
                                         Pos_2_inner_points, Pos_2_outer_points,
                                         Neg_2_inner_points, Neg_2_outer_points
                     ))
                   }
  
  result_list[[i]] <- rbindlist(result)
  result_list[[i]]$TF<-TF_name[i]
  
}
result_list_1<-list()

for (i in 1:length(result_list)) {
  temp<-result_list[[i]]
  
  temp<-temp[order(-temp$Pos_2_inner_outer_ratio),]
  dist1<-mean(temp$dist[1:5])
    #dist1<-temp$dist[1]

  temp<-temp[order(-temp$Neg_2_inner_outer_ratio),]
  dist2<-mean(temp$dist[1:5])
    #dist2<-temp$dist[1]
  
  temp<-as.data.frame(cbind(dist1,frag1=20,dist2,frag2=20))
  
  temp$TF<-result_list[[i]]$TF[i]
  result_list_1[[i]]<-temp
}








result_list_1<-rbindlist(result_list_1)


result_list_1$v_smt_dist<- (result_list_1$dist1+result_list_1$dist2)/2

result_list_1$v_smt_frag<- round(2*result_list_1$v_smt_dist+20-2*result_list_1$dist1)


v_inner_outer_FC<-list()

for (i in 1:nrow(result_list_1)) {
  dist<-result_list_1$v_smt_dist[i]
  frag<-result_list_1$v_smt_frag[i]
  
  temp<-fread(filePath[grepl(result_list_1$TF[i],filePath)])
  
  colnames(temp)<-c("motif","fragment","distance")
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
  
  v_inner_outer_FC[[i]]<-as.data.frame(cbind(v_inner_count=v_inner_count,
                                             v_outer_count=v_outer_count,
                                             FC=v_inner_count/v_outer_count,
                                             TF=result_list_1$TF[i]))
  
  
}

v_inner_outer_FC<-rbindlist(v_inner_outer_FC)          


result_list_2<-merge(result_list_1,v_inner_outer_FC,by="TF",all=T)
result_list_2$FC<-as.numeric(result_list_2$FC)
result_list_2<-result_list_2[result_list_2$TF%in%c("CTCF","ZNF148","ZNF281","PATZ1","EGR1","KLF1",
                                                   "MAZ","KLF6","SP1","SP2"),]



plot_data_DNase<-rbind(cbind(bias_result_list_2[,c("TF","v_smt_frag","FC")],
                       type="bias"),
                 cbind(result_list_2[,c("TF","v_smt_frag","FC")],
                       type="TF"))




library(ggrepel)
ggplot() + 
  
  
  
  # stat_density_2d(geom = "tile", aes(fill = ..density..), contour = FALSE,n=200) +
  # scale_fill_gradient2(low = "white", high = "#0c3e6f") +
  geom_point(data =plot_data_DNase[plot_data_DNase$type=="bias",], aes(x=v_smt_frag, y=log2(FC)),
             alpha = 1, size = 4, color = "#FF6A6A",shape=1) +
  geom_point(data =plot_data_DNase[plot_data_DNase$type=="TF",], aes(x=v_smt_frag, y=log2(FC)),
             alpha = 0.8, size = 4, color = "#104e8b",shape=16) +

  geom_text_repel(data =plot_data_DNase[plot_data_DNase$type=="TF",], 
                  aes(x=v_smt_frag, y=log2(FC),label=TF),min.segment.length = 0,max.overlaps = 20,segment.colour = NA)+
  
  scale_y_continuous(expand = c(0,0),limits = c(-0.5,3))+
  scale_x_continuous(expand = c(0,0),limits = c(-1,26),breaks = seq(0,26,2))+
  
  theme_bw()+xlab("V-channel width (bp)")+
  ylab("log2(V-in/V-out)")+
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
        plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'))



ggsave("DNase_valley_enrichment_scatterplot.pdf",width=7.5, height=3.5)