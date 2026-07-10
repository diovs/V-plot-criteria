#-------------------------------------generate bias sequence fasta file for CTCF motif
awk '{print $1,$2-12,$3+12,$4,$5,$6}' OFS='\t' CTCF_ATAC_no_bind.motif | sort -k1,1 -k2,2n -S 100%| \
 bedtools getfasta -fi /ssd/index/bismark/hg38XY/hg38XY.fa -bed - -s -name -tab \
  > CTCF_ATAC_no_bind.fa

awk '{print $1,$2-12,$3+12,$4,$5,$6}' OFS='\t' CTCF_ATAC_all_bind.motif | sort -k1,1 -k2,2n -S 100%| \
 bedtools getfasta -fi /ssd/index/bismark/hg38XY/hg38XY.fa -bed - -s -name -tab \
  > CTCF_ATAC_all_bind.fa

awk '{if($6=="+") {print $1,$2,$3+24,$4,$5,$6}else{print $1,$2-24,$3,$4,$5,$6}}' OFS='\t' MNase_Bias_top5_shuffled_no_CTCF.bed | sort -k1,1 -k2,2n -S 100%| \
 bedtools getfasta -fi /ssd/index/bismark/hg38XY/hg38XY.fa -bed - -s -name -tab \
  > MNase_Bias_top5_shuffled_no_CTCF.fa