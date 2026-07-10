#-------------------------------------generate Vplot fragL_dist file
#---------CTCF
nohup closestBed -a K562_DNase_merge_midP_fragL.bed \
 -b CTCF_DNase_all_bind.motif -d -t first | \
awk '{if($13>=0 && $13<=300){if(($9<=$3 && $12=="+")||($9>$3 && $12=="-")) print $10,$5,$13; else print $10,$5,$13*(-1)}}' OFS='\t' \
  > CTCF_DNase_all_bind_fragL_dist.txt &

nohup closestBed -a K562_DNase_merge_midP_fragL.bed \
 -b CTCF_DNase_no_bind.motif -d -t first | \
awk '{if($13>=0 && $13<=300){if(($9<=$3 && $12=="+")||($9>$3 && $12=="-")) print $10,$5,$13; else print $10,$5,$13*(-1)}}' OFS='\t' \
  > CTCF_DNase_all_bind_fragL_dist.txt &


nohup closestBed -a K562_DNase_merge_midP_fragL.bed \
 -b DNase_Bias_top5_shuffled_no_CTCF.bed -d -t first | \
awk '{if($13>=0 && $13<=300){if(($9<=$3 && $12=="+")||($9>$3 && $12=="-")) print NR,$5,$13; else print NR,$5,$13*(-1)}}' OFS='\t' \
  > DNase_Bias_top5_DNase_fragL_dist.txt &



nohup closestBed -a K562_NChIP_CTCF_merge_midP_fragL.bed \
 -b CTCF_DNase_all_bind.motif -d -t first | \
awk '{if($13>=0 && $13<=300){if(($9<=$3 && $12=="+")||($9>$3 && $12=="-")) print $10,$5,$13; else print $10,$5,$13*(-1)}}' OFS='\t' \
  > CTCF_NChIP_all_bind_fragL_dist.txt &

nohup closestBed -a K562_NChIP_CTCF_merge_midP_fragL.bed \
 -b CTCF_DNase_no_bind.motif -d -t first | \
awk '{if($13>=0 && $13<=300){if(($9<=$3 && $12=="+")||($9>$3 && $12=="-")) print $10,$5,$13; else print $10,$5,$13*(-1)}}' OFS='\t' \
  > CTCF_NChIP_no_bind_fragL_dist.txt &

nohup closestBed -a K562_NChIP_CTCF_merge_midP_fragL.bed \
 -b DNase_Bias_top5_shuffled_no_CTCF.bed -d -t first | \
awk '{if($13>=0 && $13<=300){if(($9<=$3 && $12=="+")||($9>$3 && $12=="-")) print NR,$5,$13; else print NR,$5,$13*(-1)}}' OFS='\t' \
  > DNase_Bias_top5_NChIP_fragL_dist.txt &