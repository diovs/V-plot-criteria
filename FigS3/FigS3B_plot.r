path<-"/mnt/disk2/2/zqf/private/GB_Vplot_Correspondence/TF_selected/fragL_dist"


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
  temp<-temp[between(temp$fragment,20,75),]
  
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
  
  temp$Pos_2_inner_outer_wt<-(temp$Pos_2_inner_points+temp$Pos_2_outer_points)/sum(temp$Pos_2_inner_points+temp$Pos_2_outer_points)
  temp$Pos_2_inner_outer_ratio_wt<-temp$Pos_2_inner_outer_ratio*temp$Pos_2_inner_outer_wt
  
  temp$Neg_2_inner_outer_wt<-(temp$Neg_2_inner_points+temp$Neg_2_outer_points)/sum(temp$Neg_2_inner_points+temp$Neg_2_outer_points)
  temp$Neg_2_inner_outer_ratio_wt<-temp$Neg_2_inner_outer_ratio*temp$Neg_2_inner_outer_wt
  
  temp<-temp[order(-temp$Pos_2_inner_outer_ratio_wt),]
  dist1<-mean(temp$dist[1:5])
  
  temp<-temp[order(-temp$Neg_2_inner_outer_ratio_wt),]
  dist2<-mean(temp$dist[1:5])
  
  temp<-as.data.frame(cbind(dist1,frag1=20,dist2,frag2=20))
  
  temp$TF<-result_list[[i]]$TF[i]
  result_list_1[[i]]<-temp
}










result_list_1<-rbindlist(result_list_1)



result_list_1$v_smt_dist<- (result_list_1$dist1+result_list_1$dist2)/2

result_list_1$v_smt_frag<- round(2*result_list_1$v_smt_dist+20-2*result_list_1$dist1)


p<-list()
for (i in 1:nrow(result_list_1)) {
  temp<-fread(filePath[i])
  
  colnames(temp)<-c("motif","fragment","distance")
  dist<-result_list_1$v_smt_dist[i]
  frag<-result_list_1$v_smt_frag[i]
  segment_table<-as.data.frame(cbind(x=c(dist,dist,dist-(100-frag)/2,dist-frag,dist+frag),
                                     y=c(frag,frag,100,frag,frag),
                                     xend=c(dist-(100-frag)/2,dist+(100-frag)/2,dist+(100-frag)/2,dist-frag-(100-frag)/2,dist+frag+(100-frag)/2),
                                     yend=c(100,100,100,100,100)))
  if (nrow(temp)>200000) {
    temp<-temp[sample(1:nrow(temp),200000),]
    p[[i]]<-ggplot(data =temp, aes(x=distance, y=fragment)) + 
      
      geom_pointdensity(size=0.2,show.legend = F)+
      scale_color_gradient2(low = "#88a7c5", high = "#104e8b") +
      geom_segment(data =segment_table,
                   aes(x = x, y = y, xend =xend , yend = yend),
                   size = 0.5,linetype="dashed",color="red"
      )+

    theme_bw()+xlab("Distance from motif (bp)")+
      ylab("Fragment length (bp)")+
      ggtitle(result_list_1$TF[i])+
      scale_x_continuous(
        limits = c(-60, 60),
        breaks=seq(-60,60,30),
        expand = c(0,0)) +
      scale_y_continuous(
        limits = c(0, 100),
        # breaks=seq(0,150,50),
        expand = c(0,0)) +
      theme_bw() +
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
            axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 16, colour = "black"),
            axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 16, colour = "black"),
            plot.title = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 20, colour = "black"),
            axis.ticks = element_line(color = "black"),
            text = element_text(size = 15, family = "sans", colour = "black"),
            plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm')) +coord_fixed(ratio = 1)
    print(i)
  }else{
    p[[i]]<-ggplot(data =temp, aes(x=distance, y=fragment)) + 
      
      geom_pointdensity(size=0.2,show.legend = F)+
      scale_color_gradient2(low = "#88a7c5", high = "#104e8b") +
      geom_segment(data =segment_table,
                   aes(x = x, y = y, xend =xend , yend = yend),
                   size = 0.5,linetype="dashed",color="red"
      )+

      
      theme_bw()+xlab("Distance from motif (bp)")+
      ylab("Fragment length (bp)")+
      ggtitle(result_list_1$TF[i])+
      scale_x_continuous(
        limits = c(-60, 60),
        breaks=seq(-60,60,30),
        expand = c(0,0)) +
      scale_y_continuous(
        limits = c(0, 100),
        # breaks=seq(0,150,50),
        expand = c(0,0)) +
      theme_bw() +
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
            axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 16, colour = "black"),
            axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 16, colour = "black"),
            plot.title = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 20, colour = "black"),
            axis.ticks = element_line(color = "black"),
            text = element_text(size = 15, family = "sans", colour = "black"),
            plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm')) +coord_fixed(ratio = 1)
    print(i)
  }
  
  
  
  ggsave(paste("/mnt/disk2/2/zqf/private/GB_Vplot_Correspondence/TF_selected/",result_list_1$TF[i],"_V_edge.png",sep = ""),p[[i]],width=4.5, height=3.5,dpi=300)
  
}