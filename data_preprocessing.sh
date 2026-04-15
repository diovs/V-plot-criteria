#--------------------------------------------loMNase
samtools merge -@ 80 -o K562_loMNase.bam \
 /mnt/disk4/public/publication/2407_MSB/loMNase/bam/K562_loMNase.bam /mnt/disk4/public/publication/2407_MSB/loMNase/bam/K562_loMNase_DMSO_merge.bam

samtools sort -n -@ 80 K562_loMNase.bam -o K562_loMNase_sorted.bam
bamToBed -bedpe -mate1 -i K562_loMNase_sorted.bam |\
awk '{if($9=="+"){print $1,$2,$6,$9} else{print $1,$5,$3,$9}}' OFS='\t' |\
sort -k1,1 -k2,2n -k3,3n -k4,4 -S 100% | cut -f 1-3 |\
awk '$1!~/chr[CLMT]/' OFS='\t' |awk '{if($3-$2<=100){print $0}}' OFS='\t'|\
awk '$1!~/^Mm/' OFS='\t' > K562_loMNase_merge.bed

nohup cat K562_loMNase_merge.bed 241107_loMNase_N10_q10_cutoff_100bp.bed |\
 awk '{print $1,int(($2+$3)/2),int(($2+$3)/2)+1,NR,$3-$2,NR}' OFS='\t'|\
sort -k1,1 -k2,2n -S 100% > K562_loMNase_merge_midP_fragL.bed &


nohup awk '{print $1,int(($2+$3)/2),int(($2+$3)/2)+1,NR,$3-$2,NR}' OFS='\t' K562_loMNase_merge.bed|\
sort -k1,1 -k2,2n -S 100% > K562_loMNase_merge_midP_fragL.bed &

count=`cat K562_loMNase_merge.bed | wc -l`

genomeCoverageBed -bg -i K562_loMNase_merge.bed -g /ssd/genome/hg38_chromsize.txt | \
awk -v a=${count} '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4*1000000/a}' | \
awk '{$4/=1;print}' OFS='\t' > K562_loMNase_merge.bdg
bedGraphToBigWig K562_loMNase_merge.bdg /ssd/genome/hg38_chromsize.txt K562_loMNase_merge.bw
rm K562_loMNase_merge.bdg

###----------------------------------DNase-seq---------------------------------------------
samtools sort -n -@ 80 K562_DNase_R1.bam -o K562_DNase_R1_sorted.bam
bamToBed -bedpe -mate1 -i K562_DNase_R1_sorted.bam |\
awk '{if($9=="+"){print $1,$2,$6,$9} else{print $1,$5,$3,$9}}' OFS='\t' |\
sort -k1,1 -k2,2n -k3,3n -k4,4 -S 100% | cut -f 1-3 |\
awk '$1!~/chr[CLMT]/' OFS='\t' |awk '{if($3-$2<=100){print $0}}' OFS='\t'|\
awk '$1!~/^Mm/' OFS='\t' > K562_DNase_R1_merge.bed

awk '$3-$2>25 && $3-$2<100' OFS='\t' K562_DNase_R1_merge.bed > K562_DNase_R1_frag25_100.bed

nohup awk '{print $1,int(($2+$3)/2),int(($2+$3)/2)+1,NR,$3-$2,NR}' OFS='\t' K562_DNase_R1_frag25_100.bed|\
sort -k1,1 -k2,2n -S 100% > K562_DNase_R1_midP_fragL.bed &

count=`cat K562_DNase_R1_frag25_100.bed | wc -l`

genomeCoverageBed -bg -i K562_DNase_R1_frag25_100.bed -g /ssd/genome/hg38_chromsize.txt | \
awk -v a=${count} '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4*1000000/a}' | \
awk '{$4/=1;print}' OFS='\t' > K562_DNase_R1_frag25_100.bdg
bedGraphToBigWig K562_DNase_R1_frag25_100.bdg /ssd/genome/hg38_chromsize.txt K562_DNase_R1_frag25_100.bw

###----------------------------------ATAC-seq---------------------------------------------
samtools sort -n -@ 80 K562_ATAC_R1.bam -o K562_ATAC_R1_sorted.bam
bamToBed -bedpe -mate1 -i K562_ATAC_R1_sorted.bam |\
awk '{if($9=="+"){print $1,$2,$6,$9} else{print $1,$5,$3,$9}}' OFS='\t' |\
sort -k1,1 -k2,2n -k3,3n -k4,4 -S 100% | cut -f 1-3 |\
awk '$1!~/chr[CLMT]/' OFS='\t' |awk '{if($3-$2<=100){print $0}}' OFS='\t'|\
awk '$1!~/^Mm/' OFS='\t' > K562_ATAC_R1_merge.bed


awk '$3-$2>25 && $3-$2<100' OFS='\t' K562_ATAC_R1_merge.bed |sort -k1,1 -k2,2n -S 100% > K562_ATAC_R1_frag25_100.bed

nohup awk '{print $1,int(($2+$3)/2),int(($2+$3)/2)+1,NR,$3-$2,NR}' OFS='\t' K562_ATAC_R1_frag25_100.bed |\
sort -k1,1 -k2,2n -S 100% > K562_ATAC_R1_midP_fragL.bed &

count=`cat K562_ATAC_R1_frag25_100.bed | wc -l`

genomeCoverageBed -bg -i K562_ATAC_R1_frag25_100.bed -g /ssd/genome/hg38_chromsize.txt | \
awk -v a=${count} '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4*1000000/a}' | \
awk '{$4/=1;print}' OFS='\t' > K562_ATAC_R1_frag25_100.bdg
bedGraphToBigWig K562_ATAC_R1_frag25_100.bdg /ssd/genome/hg38_chromsize.txt K562_ATAC_R1_frag25_100.bw







#---------------X-ChIP
#-----------------------------------------CTCF

samtools sort -n -@ 80 /mnt/disk4/public/publication/2407_MSB/X-ChIP/bam/K562_X-ChIP_CTCF_merge.bam -o K562_X-ChIP_CTCF_merge_sorted.bam
bamToBed -bedpe -mate1 -i K562_X-ChIP_CTCF_merge_sorted.bam |\
awk '{if($9=="+"){print $1,$2,$6,$9} else{print $1,$5,$3,$9}}' OFS='\t' |\
sort -k1,1 -k2,2n -k3,3n -k4,4 -S 100% | cut -f 1-3 |\
awk '$1!~/chr[CLMT]/' OFS='\t' |\
awk '$1!~/^Mm/' OFS='\t' > K562_X-ChIP_CTCF_merge.bed

nohup awk '{print $1,int(($2+$3)/2),int(($2+$3)/2)+1,NR,$3-$2,NR}' OFS='\t' K562_X-ChIP_CTCF_merge.bed|\
sort -k1,1 -k2,2n -S 100% > K562_X-ChIP_CTCF_merge_midP_fragL.bed &

count=`cat K562_X-ChIP_CTCF_merge.bed | wc -l`

genomeCoverageBed -bg -i K562_X-ChIP_CTCF_merge.bed -g /ssd/genome/hg38_chromsize.txt | \
awk -v a=${count} '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4*1000000/a}' | \
awk '{$4/=1;print}' OFS='\t' > K562_X-ChIP_CTCF_merge.bdg
bedGraphToBigWig K562_X-ChIP_CTCF_merge.bdg /ssd/genome/hg38_chromsize.txt K562_X-ChIP_CTCF.bw
rm K562_X-ChIP_CTCF_merge.bdg

#-----------------------------------------NChIP
samtools sort -n -@ 80 K562_CTCF_NChIP.bam -o K562_CTCF_NChIP_sorted.bam
bamToBed -bedpe -mate1 -i K562_CTCF_NChIP_sorted.bam |\
awk '{if($9=="+"){print $1,$2,$6,$9} else{print $1,$5,$3,$9}}' OFS='\t' |\
sort -k1,1 -k2,2n -k3,3n -k4,4 -S 100% | cut -f 1-3 |\
awk '$1!~/chr[CLMT]/' OFS='\t' |awk '{if($3-$2<=100){print $0}}' OFS='\t'|\
awk '$1!~/^Mm/' OFS='\t' > K562_CTCF_NChIP_merge.bed

nohup cat K562_CTCF_NChIP_merge.bed 241107_CTCF_N10_q10_cutoff_100bp.bed |\
 awk '{print $1,int(($2+$3)/2),int(($2+$3)/2)+1,NR,$3-$2,NR}' OFS='\t'|\
sort -k1,1 -k2,2n -S 100% > K562_CTCF_NChIP_merge_midP_fragL.bed &


nohup awk '{print $1,int(($2+$3)/2),int(($2+$3)/2)+1,NR,$3-$2,NR}' OFS='\t' K562_CTCF_NChIP_merge.bed|\
sort -k1,1 -k2,2n -S 100% > K562_CTCF_NChIP_merge_midP_fragL.bed &

count=`cat K562_CTCF_NChIP_merge.bed | wc -l`

genomeCoverageBed -bg -i K562_CTCF_NChIP_merge.bed -g /ssd/genome/hg38_chromsize.txt | \
awk -v a=${count} '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4*1000000/a}' | \
awk '{$4/=1;print}' OFS='\t' > K562_CTCF_NChIP_merge.bdg
bedGraphToBigWig K562_CTCF_NChIP_merge.bdg /ssd/genome/hg38_chromsize.txt K562_CTCF_NChIP_merge.bw
rm K562_CTCF_NChIP_merge.bdg


#-----------------------------------------NChIP
samtools sort -n -@ 80 K562_MAZ_NChIP.bam -o K562_MAZ_NChIP_sorted.bam
bamToBed -bedpe -mate1 -i K562_MAZ_NChIP_sorted.bam |\
awk '{if($9=="+"){print $1,$2,$6,$9} else{print $1,$5,$3,$9}}' OFS='\t' |\
sort -k1,1 -k2,2n -k3,3n -k4,4 -S 100% | cut -f 1-3 |\
awk '$1!~/chr[CLMT]/' OFS='\t' |awk '{if($3-$2<=100){print $0}}' OFS='\t'|\
awk '$1!~/^Mm/' OFS='\t' > K562_MAZ_NChIP_merge.bed

nohup cat K562_MAZ_NChIP_merge.bed 241107_MAZ_N10_q10_cutoff_100bp.bed |\
 awk '{print $1,int(($2+$3)/2),int(($2+$3)/2)+1,NR,$3-$2,NR}' OFS='\t'|\
sort -k1,1 -k2,2n -S 100% > K562_MAZ_NChIP_merge_midP_fragL.bed &


nohup awk '{print $1,int(($2+$3)/2),int(($2+$3)/2)+1,NR,$3-$2,NR}' OFS='\t' K562_MAZ_NChIP_merge.bed|\
sort -k1,1 -k2,2n -S 100% > K562_MAZ_NChIP_merge_midP_fragL.bed &

count=`cat K562_MAZ_NChIP_merge.bed | wc -l`

genomeCoverageBed -bg -i K562_MAZ_NChIP_merge.bed -g /ssd/genome/hg38_chromsize.txt | \
awk -v a=${count} '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4*1000000/a}' | \
awk '{$4/=1;print}' OFS='\t' > K562_MAZ_NChIP_merge.bdg
bedGraphToBigWig K562_MAZ_NChIP_merge.bdg /ssd/genome/hg38_chromsize.txt K562_MAZ_NChIP_merge.bw
rm K562_MAZ_NChIP_merge.bdg