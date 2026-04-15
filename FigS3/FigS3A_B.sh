
ls /mnt/disk2/2/zqf/private/GB_Vplot_Correspondence/TF_selected/*.motif |while read id;
do
filename=`basename $id .motif`
 nohup closestBed -a K562_loMNase_merge_midP_fragL.bed -b ${filename}.motif -d -t first | \
 awk '{if($13>=0 && $13<=300){if(($9<=$3 && $12=="+")||($9>$3 && $12=="-")) print $10,$5,$13; else print $10,$5,$13*(-1)}}' OFS='\t' \
  > ${filename}_midP_fragL_dist.txt &
done



for motif in AGG AGC TGC TGG AAG; do


 awk '{print $1,$2-12,$3+12,$4,$5,$6}' OFS='\t'  Hs_CTCF.motif |\
 intersectBed -sorted -v -s -f 1 -a ${motif}_shuffled.bed -b - -wa \
  > ${motif}_shuffled_no_CTCF.bed

 nohup closestBed -a K562_loMNase_merge_midP_fragL.bed \
 -b ${motif}_shuffled_no_CTCF.bed -d -t first | \
 awk '{if($13>=0 && $13<=300){if(($9<=$3 && $12=="+")||($9>$3 && $12=="-")) print $10,$5,$13; else print $10,$5,$13*(-1)}}' OFS='\t' \
  > ${motif}_loMNase_fragL_dist.txt &

done