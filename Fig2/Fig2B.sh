nohup closestBed -a K562_loMNase_merge_midP_fragL.bed -b Hs_MAZ.motif -d -t first | \
 awk '{if($13>=0 && $13<=300){if(($9<=$3 && $12=="+")||($9>$3 && $12=="-")) print $10,$5,$13; else print $10,$5,$13*(-1)}}' OFS='\t' \
  > MAZ_loMNase_midP_fragL_dist.txt &

awk '{if($2>20){print $1,$2-10,$3+10,$4,0,$6}}' OFS='\t' Hs_MAZ.motif  > Hs_MAZ_20bp.bed



nohup bigWigAverageOverBed K562_loMNase_merge.bw Hs_MAZ_20bp.bed MAZ_K562_loMNase_20.tab &

nohup bigWigAverageOverBed K562_MAZ_ab.bw Hs_MAZ_20bp.bed MAZ_K562_MAZ_XChIP_20.tab &
nohup bigWigAverageOverBed K562_N-ChIP_MAZ_50mM_merge.bw Hs_MAZ_20bp.bed MAZ_K562_MAZ_NChIP_20.tab &

nohup bigWigAverageOverBed K562_DNase_R1_frag25_100.bw Hs_MAZ_20bp.bed MAZ_K562_DNase_20.tab &
nohup bigWigAverageOverBed K562_ATAC_R1_frag25_100.bw Hs_MAZ_20bp.bed MAZ_K562_ATAC_20.tab &




join -1 1 -2 1 MAZ_K562_MAZ_XChIP_20.tab MAZ_K562_MAZ_NChIP_20.tab |\
join -1 1 -2 1 - MAZ_K562_loMNase_20.tab |\
join -1 1 -2 1 - MAZ_K562_DNase_20.tab |\
join -1 1 -2 1 - MAZ_K562_ATAC_20.tab |\
 sort -k1,1n -S 100% |join -1 4 -2 1 Hs_MAZ_20bp.bed - | awk '{print $2,$3+10,$4-10,$1,$5,$6,$10,$15,$20,$25,$30}' OFS='\t'|\
 sort -k1,1 -k2,2n -S 100% > MAZ_20bp_motif_with_K562_signal.bed

rm *.tab MAZ_20bp.bed



ls /mnt/disk2/2/zqf/private/GB_Vplot_Correspondence/TF_selected/MAZ_*.motif |while read id;
do
filename=`basename $id .motif`

nohup computeMatrix reference-point -S /mnt/disk2/2/zqf/private/GB_Vplot_Correspondence/other_data/K562_N-ChIP_MAZ_50mM_merge.bw \
-R $id -a 2000 -b 2000 --referencePoint center -bs 20 \
-o ${filename}_NChIP.matrix.txt.gz -p max --missingDataAsZero &

done