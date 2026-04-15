nohup closestBed -a K562_loMNase_merge_midP_fragL.bed \
 -b CTCF.motif -d -t first | \
awk '{if($13>=0 && $13<=300){if(($9<=$3 && $12=="+")||($9>$3 && $12=="-")) print $10,$5,$13; else print $10,$5,$13*(-1)}}' OFS='\t' \
  > CTCF_loMNase_fragL_dist.txt &



awk '{if($2>160){print $1,$2-160,$3+160,$4,$5,$6}}' OFS='\t' CTCF.motif > CTCF_selected_320bp.bed


awk '{print $1,$2,$3}' OFS='\t' CTCF_fragments_500per.bed |sort -k1,1 -k2,2n -S 60% > CTCF_fragments_random.bed

awk '{if($2>100){print $1,$2-10,$3+10,$4,$5,$6}}' OFS='\t' CTCF.motif |\
intersectBed -sorted -v -a CTCF_fragments_random.bed -b - -wa  > CTCF_fragments_random_no_CTCF_overlap.bed
awk '{if($2>100){print $1,$2-10,$3+10,$4,$5,$6}}' OFS='\t' CTCF.motif |\
intersectBed -sorted -F 1 -a CTCF_fragments_random.bed -b - -wa  > CTCF_fragments_random_full_CTCF_overlap.bed

count=`cat CTCF_fragments_random_full_CTCF_overlap.bed | wc -l`

shuf -n $count/2 CTCF_fragments_random_no_CTCF_overlap.bed |\
cat - CTCF_fragments_random_full_CTCF_overlap.bed |sort -k1,1 -k2,2n -S 60% > CTCF_fragments_random_all.bed
rm CTCF_fragments_random_no_CTCF_overlap.bed CTCF_fragments_random_full_CTCF_overlap.bed

awk '{print $1,int(($2+$3)/2),int(($2+$3)/2)+1,NR,$3-$2,NR}' OFS='\t' CTCF_fragments_random_all.bed |\
sort -k1,1 -k2,2n -S 100% > CTCF_fragments_random_all_midP_fragL.bed
closestBed -a CTCF_fragments_random_all_midP_fragL.bed -b CTCF.motif -d -t first | \
awk '{if($13>=0 && $13<=300){if(($9<=$3 && $12=="+")||($9>$3 && $12=="-")) print $10,$5,$13; else print $10,$5,$13*(-1)}}' OFS='\t' \
  > CTCF_fragments_random_all_midP_fragL_CTCF_dist.txt