
#-------------------------------------MNase fragment end bigwig
#------CTCF
awk '{print $1,$2,$2+1}' OFS='\t' K562_ATAC_merge.bed > K562_ATAC_end_5.bed
awk '{print $1,$3-1,$3}' OFS='\t' K562_ATAC_merge.bed > K562_ATAC_end_3.bed

cat K562_ATAC_end_5.bed K562_ATAC_end_3.bed | sort -k1,1 -k2,2n -S 100% > K562_ATAC_end.bed

 
awk '{print $1,$2-60,$3+60,$4,$5,$6}' OFS='\t'  /mnt/disk4/public/RefBed/CTCF/Hs_CTCF.motif |\
intersectBed -sorted -f 1 -a K562_ATAC_end.bed -b - -wa -wb | \
awk '{print $1,$2,$3,$4,$5+60,$6-60,$7,$8,$9}' OFS='\t' |\
awk '{if($9=="+" && $2<=$5){print $1,$2,$3}else if($9=="-" && $2>$5){print $1,$2,$3}}' OFS='\t' \
> K562_ATAC_end_up.bed

awk '{print $1,$2-60,$3+60,$4,$5,$6}' OFS='\t'  /mnt/disk4/public/RefBed/CTCF/Hs_CTCF.motif |\
intersectBed -sorted -f 1 -a K562_ATAC_end.bed -b - -wa -wb | \
awk '{print $1,$2,$3,$4,$5+60,$6-60,$7,$8,$9}' OFS='\t' |\
awk '{if($9=="+" && $2>=$5){print $1,$2,$3}else if($9=="-" && $2<$5){print $1,$2,$3}}' OFS='\t' \
> K562_ATAC_end_down.bed

count=`cat K562_ATAC_end_up.bed | wc -l`

genomeCoverageBed -bg -i K562_ATAC_end_up.bed -g /ssd/genome/hg38_chromsize.txt | \
awk -v a=${count} '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4*1000000/a}' | \
awk '{$4/=1;print}' OFS='\t' > K562_ATAC_end_up.bdg
bedGraphToBigWig K562_ATAC_end_up.bdg /ssd/genome/hg38_chromsize.txt K562_ATAC_end_up.bw
rm K562_ATAC_end_up.bdg

count=`cat K562_ATAC_end_down.bed | wc -l`

genomeCoverageBed -bg -i K562_ATAC_end_down.bed -g /ssd/genome/hg38_chromsize.txt | \
awk -v a=${count} '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4*1000000/a}' | \
awk '{$4/=1;print}' OFS='\t' > K562_ATAC_end_down.bdg
bedGraphToBigWig K562_ATAC_end_down.bdg /ssd/genome/hg38_chromsize.txt K562_ATAC_end_down.bw
rm K562_ATAC_end_down.bdg


computeMatrix reference-point -S K562_ATAC_end_up.bw K562_ATAC_end_down.bw \
-R /mnt/disk2/2/zqf/private/GB_Vplot_Correspondence/motif_intersect_bed/CTCF_ATAC_no_bind.motif -a 60 -b 60 --referencePoint center -bs 1 \
-o K562_ATAC_5_3_no_bind.matrix.txt.gz -p max --missingDataAsZero

computeMatrix reference-point -S K562_ATAC_end_up.bw K562_ATAC_end_down.bw \
-R /mnt/disk2/2/zqf/private/GB_Vplot_Correspondence/motif_intersect_bed/CTCF_ATAC_all_bind.motif -a 60 -b 60 --referencePoint center -bs 1 \
-o K562_ATAC_5_3_all_bind.matrix.txt.gz -p max --missingDataAsZero

shuf -n 100000 /mnt/disk2/2/zqf/private/GB_Vplot_Correspondence/MNase_Bias_top5_shuffled_no_CTCF.bed|\
awk '{print $1,$2,$3,NR,$5,$6}' OFS='\t'|\
sort -k1,1 -k2,2n -S 100% > /mnt/disk2/2/zqf/private/GB_Vplot_Correspondence/motif_intersect_bed/MNase_Bias_top5_shuffled_no_CTCF_100k.bed

computeMatrix reference-point -S K562_ATAC_end_up_MNase_Bias.bw K562_ATAC_end_down_MNase_Bias.bw \
-R /mnt/disk2/2/zqf/private/GB_Vplot_Correspondence/motif_intersect_bed/MNase_Bias_top5_shuffled_no_CTCF_100k.bed -a 60 -b 60 --referencePoint center -bs 1 \
-o K562_ATAC_5_3_MNase_Bias_top5.matrix.txt.gz -p max --missingDataAsZero